# library(devtools)
# devtools::install("/juno/work/iacobuzc/haochen/bulk_analysis/tools/signature.tools.lib")
# BiocManager::install("BSgenome")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("VariantAnnotation")

# ‘VariantAnnotation’, 
# ‘BSgenome.Hsapiens.UCSC.hg38’, 
# ‘BSgenome.Hsapiens.1000genomes.hs37d5’, 
# ‘BSgenome’, 
# ‘limSolve’ are not available for package ‘signature.tools.lib’

library(signature.tools.lib)
library(data.table)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(glue)
library(ggplot2)
# library(ggpubr)
library(RColorBrewer)
library(stringr)
ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5
genome_version = "hg19"

setwd("/juno/work/iacobuzc/haochen/bulk_analysis/signatures_analysis/signatures")

# ===== do the below in bash to construct input TSV =====
# echo -e "sample\tsv\tmutations\tcnv" > s_RA17_13_67_1__s_RA17_13_2_1.tsv
# echo -e "s_RA17_13_67_1__s_RA17_13_2_1\ts_RA17_13_67_1__s_RA17_13_2_1.final.bedpe\ts_RA17_13_67_1__s_RA17_13_2_1.somatic.maf\ts_RA17_13_67_1__s_RA17_13_2_1_hisens.facets.filtered.copynumber.csv" >> s_RA17_13_67_1__s_RA17_13_2_1.tsv
# Rscript HRDetect_wrapper.R s_RA17_13_67_1__s_RA17_13_2_1.tsv hg19 2
input_tsv <- "WGS_tempo_pipeline_output_fs.tsv"
input.files<-fread(input_tsv, header = T, data.table = T)
# fill NA with empty string
input.files[is.na(input.files)] <- ""
sample_names<-input.files$sample

dir.create("intermediate")

# ----- no need to rerun now that we have "intermediate/" created -----
correctMutations <- function(this_sample){
  print(paste0("processing ", this_sample))
  if (file.exists(paste0("intermediate/", this_sample, ".indels"))){
    # print that this sample's intermediates are already generated
    print("sample {this_sample} already processed! Skipping")
    return()
  }
  # if the input file is empty, skip
  if (input.files[input.files$sample==this_sample,]$mutations == ""){
    print("[WARNING] empty input file for sample {this_sample}! Skipping...")
    return()
  }
  this_mutations<-fread(input.files[input.files$sample==this_sample,]$mutations, header = T, data.table = F)
  setnames(this_mutations, "Chromosome", "chr")
  setnames(this_mutations, "vcf_pos", "position")
  setnames(this_mutations, "vcf_id", "ID")
  setnames(this_mutations, "Reference_Allele", "REF")
  setnames(this_mutations, "Allele", "ALT")
  setnames(this_mutations, "vcf_qual", "QUAL")
  this_indels<-this_mutations[this_mutations$Variant_Type %in% c("DEL", "DNP", "INS", "TNP"),]
  this_indels <- this_indels %>%
      dplyr::rowwise() %>%
      mutate(left_b = as.character(ref.genome[[chr]][position])) %>%
      mutate(REF = ifelse(ALT == "-",paste0(left_b,REF),REF),
	     ALT = ifelse(ALT == "-",left_b,ALT)) %>%
      mutate(ALT = ifelse(REF == "-",paste0(left_b,ALT),ALT),
	     REF = ifelse(REF == "-",left_b,REF)) %>%
      select(-c(left_b))

  this_snv<-this_mutations[!this_mutations$Variant_Type %in% c("DEL", "DNP", "INS", "TNP"),]
  print(table(rbind(this_indels, this_snv)$Variant_Type))
  write.table(file = paste0("intermediate/", this_sample, ".indels"), x = this_indels, quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(file = paste0("intermediate/", this_sample, ".snv"), x = this_snv, quote = F, row.names = F, col.names = T, sep = "\t")
}

lapply(sample_names,correctMutations)
# -----

# ----- load the intermediates -----
Indels_tab_files <- paste0("intermediate/",sample_names,".indels")
SNV_tab_files <- paste0("intermediate/",sample_names,".snv")
names(Indels_tab_files) <- sample_names
names(SNV_tab_files) <- sample_names

# #load SNV data and convert to SNV mutational catalogues
# SNVcat_list <- list()
# message("")
# message("========= Converting SNV mutational catalogues =========")
# for (i in 1:length(SNV_tab_files)){
#   message(sample_names[i])
#   tmpSNVtab <- read.table(SNV_tab_files[i],sep = "\t", fill = T, quote="", header = TRUE,check.names = FALSE, stringsAsFactors = FALSE)
#   res <- tabToSNVcatalogue(subs = tmpSNVtab,genome.v = genome_version)
#   colnames(res$catalogue) <- sample_names[i]
#   SNVcat_list[[i]] <- res$catalogue
# }
# SNV_catalogues <- do.call(cbind,SNVcat_list)

#Initialize feature matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = length(sample_names), ncol = length(col_hrdetect), dimnames = list(sample_names,col_hrdetect))

# Compute the proportion of indels at micro-homology 
Indel.del.mh.prop_list <- list()
Indel.del.mh.prop_list.split_clonality <- list()
case_names_clonality <- list()
sample_names_clonality <- list()
message("")
message("========= Computing the proportion of indels at micro-homology =========")
for (i in 1:length(Indels_tab_files)){
  message(sample_names[i])
  tmpIndeltab <- read.table(Indels_tab_files[i],sep = "\t", fill = T, quote="", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # filter out sites where REF or ALT is empty
  tmpIndeltab <- tmpIndeltab[!(tmpIndeltab$REF == "") & !(tmpIndeltab$ALT == ""),]
  res <- tabToIndelsClassification(tmpIndeltab, sample_names[i], genome_version)
  Indel.del.mh.prop_list[[i]] <- res$count_proportion

  # # for each sample, split the mutations classified in column "clonality" into "CLONAL" and "SUBCLONAL"
  # # and compute the proportion of clonal mutations with microhomology
  # # and the proportion of subclonal mutations with microhomology
  # res_clonal <- tabToIndelsClassification(tmpIndeltab[tmpIndeltab$clonality == "CLONAL",], paste0(sample_names[i], "_clonal"), genome_version)
  # res_subclonal <- tabToIndelsClassification(tmpIndeltab[tmpIndeltab$clonality == "SUBCLONAL",], paste0(sample_names[i], "_subclonal"), genome_version)
  # Indel.del.mh.prop_list.split_clonality[[2*i-1]] <- res_clonal$count_proportion
  # Indel.del.mh.prop_list.split_clonality[[2*i]] <- res_subclonal$count_proportion

  # append the sample name to the list
  # duplicate case names
  case_names_clonality <- c(case_names_clonality, sample_names[i], sample_names[i])

  sample_names_clonality[[2*i-1]] <- paste0(sample_names[i], ".clonal")
  sample_names_clonality[[2*i]] <- paste0(sample_names[i], ".subclonal")
}
names(Indel.del.mh.prop_list)<-sample_names
# names(Indel.del.mh.prop_list.split_clonality) <- sample_names_clonality
Indel.del.mh.prop <- do.call(rbind,Indel.del.mh.prop_list)
# Indel.del.mh.prop.split_clonality <- do.call(rbind, Indel.del.mh.prop_list.split_clonality)

# write output TSV
# fname would be "WGS_indel_results.tsv" + today's date
fname <- paste0("WGS_indel_results.", format(Sys.Date(), "%Y-%m-%d"), ".tsv")
write.table(file = fname, x = Indel.del.mh.prop, quote = F, row.names = T, col.names = T, sep = "\t")
# write.table(file = "WGS_indel_results_split_by_clonality.tsv", x = Indel.del.mh.prop.split_clonality, quote = F, row.names = T, col.names = T, sep = "\t")

# ===== Plot indel proportions =====
# hrdetect_output <- "tempo_pipeline/analysis/HRD/HRDetect_output.tsv"
# hrdetect_output <- read.table(hrdetect_output, header=T, sep="\t")

# plot barplot of sample vs del.mh.prop
# Indel.del.mh.prop$sample <- factor(Indel.del.mh.prop$sample, levels=c("M04-3", "RSX582", "RSX640", "RA15_6", "RA16_8", "RA17_13", "RA17_22", "BPA-1-RSX", "BPA-1-IR", "BPA-2-RSX"))

# BRCA2 status:
# M04-3 -> BRCA2 biallelic inactivated
# RSX582 -> BRCA2 intact
# RSX640 -> BRCA2 biallelic inactivated
# RA17_22 -> BRCA2 intact
# RA15_6 -> BRCA2 intact
# RA17_13 -> BRCA2 intact
# RA16_8 -> BRCA2 intact
# M13-1 --> BRCA2 intact
# BPA-1-RSX -> BRCA2 biallelic inactivated
# BPA-1-IR -> BRCA2 biallelic inactivated
# BPA-2-RSX -> BRCA2 biallelic inactivated
# BPA-2-IR -> BRCA2 biallelic inactivated
# BPA-4-RSX -> BRCA2 biallelic inactivated
# add to the table
Indel.del.mh.prop %>% 
  mutate(BRCA2_status = factor(c("biallelic inactivated", "intact", "biallelic inactivated", "intact", "intact", "intact", "intact", "intact","biallelic inactivated", "biallelic inactivated", "biallelic inactivated", "biallelic inactivated", "biallelic inactivated"))) %>%
  mutate(BRCA2_status = factor(BRCA2_status, levels=c("biallelic inactivated", "intact"))) %>% 
  # sort by BRCA2 status, then sample_name
  arrange(BRCA2_status) -> wgs_sample_indel_props

# sort the sample column by del.mh.prop value and make into factor
wgs_sample_indel_props$sample <- factor(wgs_sample_indel_props$sample, levels=wgs_sample_indel_props$sample[order(wgs_sample_indel_props$del.mh.prop)])

hrdetect_plot <- ggplot(wgs_sample_indel_props) + 
  aes(x=sample, y=del.mh.prop, fill=BRCA2_status) +
  geom_bar(stat="identity") + 
  theme(aspect.ratio = 1) + 
  theme_minimal() +
  geom_text(aes(label=round(del.mh.prop, 2)), vjust=-0.5, size=3) +
  labs(title="Proportion of deletions with microhomology",y="Proportion") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    # disable x axis title
    axis.title.x = element_blank(),
    # change all fonts to 12, Arial
    text = element_text(size=24),
    )

ggsave("indel_at_mh_site_props.pdf", hrdetect_plot, width=10, height=10, units="in")
ggsave("indel_at_mh_site_props.png", hrdetect_plot, width=10, height=10, units="in", dpi=300)

# ========== indel props split by clonality ============
# sort the sample column by del.mh.prop value and make into factor
Indel.del.mh.prop.split_clonality$sample <- factor(Indel.del.mh.prop.split_clonality$sample, levels=Indel.del.mh.prop.split_clonality$sample[order(Indel.del.mh.prop.split_clonality$del.mh.prop)])
# add case_name, which is to split the last part of the sample name by "_"
# Indel.del.mh.prop.split_clonality$sample <- sample_names_clonality
Indel.del.mh.prop.split_clonality$case_name <- as.character(case_names_clonality)
# as factor

# clonal status
Indel.del.mh.prop.split_clonality$clonality <- ifelse(grepl("subclonal", Indel.del.mh.prop.split_clonality$sample), "subclonal", "clonal")

hrdetect_plot_split_clonality <- ggplot(Indel.del.mh.prop.split_clonality) + 
  aes(y=del.mh.prop, x = case_name, fill=clonality) +
  geom_bar(position="dodge", stat="identity") + 
  theme(aspect.ratio = 1) + 
  theme_minimal() +
  geom_text(aes(label=round(del.mh.prop, 2)), vjust=-0.5, size=3) +
  labs(title="Proportion of deletions with microhomology",y="Proportion") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    # disable x axis title
    axis.title.x = element_blank(),
    # change all fonts to 12, Arial
    text = element_text(size=24),
    )

ggsave("indel_at_mh_site_props.split_by_clonality.pdf", hrdetect_plot_split_clonality, width=10, height=10, units="in")
ggsave("indel_at_mh_site_props.split_by_clonality.png", hrdetect_plot_split_clonality, width=10, height=10, units="in", dpi=300)

# ----- quick test -----
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

ggplot(Indel.del.mh.prop.split_clonality, aes(fill=clonality, y=del.mh.prop, x=case_name)) + 
    geom_bar(position="dodge", stat="identity")

# ========== SBS signatures ==========
# fetch all sample_data in the wd
tempo_pipeline_dir <- ("/juno/work/iacobuzc/haochen/bulk_analysis/tempo_pipeline_cmo/results/somatic")
sample_metadata_fs <- Sys.glob(paste0(tempo_pipeline_dir, "/*/meta_data/*sample_data.txt"))
# read all sample_data, rbind
sample_data <- do.call(rbind, lapply(sample_metadata_fs, read.table, header=T, sep="\t"))

# ----- also concatenate the controls -----
# @HZ TODO

# make stacked barcharts, one bar per sample
# filter out signatures with too small props, and normalize the others' proportions to 1

cols_of_interest <- c("sample")
# select only columns with ".observed"
cols_of_interest <- c(cols_of_interest, grep(".observed", colnames(sample_data), value=T))
sample_data <- sample_data[,cols_of_interest]
sample_data_long <- sample_data %>% pivot_longer(
  !sample, names_to = "signature", values_to = "proportion"
  )
# get rid of ".observed" in signature names
sample_data_long$signature <- gsub(".observed", "", sample_data_long$signature)
# filter out signatures with proportions smaller than 0.05
sample_data_long <- sample_data_long[sample_data_long$proportion > 0.05,]
# normalize the rest, within each sample
sample_data_long <- sample_data_long %>% group_by(sample) %>% mutate(proportion = proportion/sum(proportion))

# order the SBS by string sort after stripping the "SBS"
sbs_order <- c("SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS7a", "SBS8", "SBS12", "SBS13","SBS35","SBS39", "SBS40", "SBS44", "SBS50", "SBS51")

sample_data_long$signature <- factor(sample_data_long$signature, levels=sbs_order)
# rename the SBS
recode(
  sample_data_long$signature ,
  "SBS1"= "SBS1 (clock-like)",
  "SBS2"= "SBS2 (APOBEC cytidine deaminases)",
  "SBS3"= "SBS3 (defective HRD)",
  "SBS4"= "SBS4 (tobacco smoking)",
  "SBS5"= "SBS5 (clock-like)",
  "SBS7a"= "SBS7a (UV)",
  "SBS8"= "SBS8 (defective HRD)",
  "SBS12"= "SBS12 (clock-like)",
  "SBS13"= "SBS13 (APOBEC cytidine deaminases)",
  "SBS35"= "SBS35 (platinum chemo treatment)",
  "SBS39"= "SBS39 (unknown)",
  "SBS40"= "SBS40 (clock-like)",
  "SBS44"= "SBS44 (defective MMR)",
  "SBS50"= "SBS50 (possible sequencing artifact)",
  "SBS51"= "SBS51 (possible sequencing artifact)"
) -> sample_data_long$signature

# Add more colors to this palette :
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(length(unique(sample_data_long$signature)))
# create map
sig_colors <- setNames(mycolors, unique(sample_data_long$signature))
# highlight the HRD signatures
sig_colors["SBS3 (defective HRD)"] <- "#FF0000"
sig_colors["SBS8 (defective HRD)"] <- "#FF0000"

sig_compo_plot <- ggplot(sample_data_long) + 
  aes(x=sample, y=proportion, fill=signature) +
  geom_bar(stat="identity") + 
  # highlight the HRD signatures
  scale_fill_manual(values=sig_colors) +
  labs(title="SBS signatures composition", y="Proportion") +
  # width and height ratio
  theme(aspect.ratio = 2) + 
  theme_minimal() +
  # vertical x axis labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  # reduce legend font
  theme(legend.key.size = unit(0.4, "cm"), legend.text=element_text(size=8))
# save
ggsave("snv_signatures.pdf", sig_compo_plot, width=5, height=5, units="in")


# ========== indel enrichment at microhomo ==========
hrdetect_output <- "tempo_pipeline/analysis/HRD/HRDetect_output.tsv"
hrdetect_output <- read.table(hrdetect_output, header=T, sep="\t")

# plot barplot of sample vs del.mh.prop
hrdetect_output$sample <- factor(hrdetect_output$sample, levels=c("M04-3", "RSX582", "RSX640", "RA15_6", "RA16_8", "RA17_13", "RA17_22"))

hrdetect_plot <- ggplot(hrdetect_output) + 
  aes(x=sample, y=del.mh.prop, fill="#c73232") +
  geom_bar(stat="identity",  color="black") + 
  theme(aspect.ratio = 1) + 
  guides(fill=FALSE)


ggsave("tempo_pipeline/analysis/HRD/del.mh.prop.pdf", hrdetect_plot, width=5, height=5, units="in")
ggsave("tempo_pipeline/analysis/HRD/del.mh.prop.png", hrdetect_plot, width=5, height=5, units="in", dpi=300)
