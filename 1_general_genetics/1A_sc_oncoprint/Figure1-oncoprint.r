library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(dplyr)
library(openxlsx)
library(tidyr)
library(readr)
library(stringr)
library(ggnewscale)
library(patchwork)

# change directory into the same folder as the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ----- 0. per-case clinical info table -----
patient_info_table <- read.xlsx("../../Tapestri_batch2_samples_MASTER.xlsx", sheet = "all_case_genetics", startRow = 2, colNames = TRUE)
# preprocess 
patient_info_table <- patient_info_table %>% filter(censor == 0) %>% rename(patient_name = case_ID)
# str_trim everything
patient_info_table <- patient_info_table %>% mutate_all(str_trim)

# # adjust survival column
# patient_info_table$os <- as.numeric(patient_info_table$Overall.survivial)

composite_maf %>% mutate(patient_name = factor(as.character(patient_name), levels = patient_info_table$patient_name)) -> composite_maf

# ----- 1. per-case SNV MAF -----

combined_sc_maf <- read.table("pan_cohort.scMAF.csv", header = TRUE, sep = ",")
# only keep SNVs with > 0.01 CCF; scale mean AF to 0-1
combined_sc_maf <- combined_sc_maf %>%
  filter(CCF > 0.01) %>% mutate(mean_sc_AF_in_mut_filtered = (mean_sc_AF_in_mut_filtered /100) )
genes_of_interest <- c("KRAS", "CDKN2A", "TP53", "SMAD", "TGFB", "ACVR", "ARID", "ATM", "POL", "CREBBP", "IRF6", "PIK3CA","STK11", "FGFR1", "MYC","FBXW7", "RNF43", "BRCA1", "BRCA2", "PALB2")

# +++++ blacklist +++++
blacklist_df <- read.csv("tumor_pon/all_tumor_cases.pon_blacklist.csv", header=TRUE)

# save into an array
blacklist_snvs <- as.list(blacklist_df$X)

snv_maf_for_plotting <- NULL
for (gene in genes_of_interest) {
  snv_maf_for_plotting <- rbind(snv_maf_for_plotting, combined_sc_maf %>% filter(grepl(gene, gene_name, fixed=TRUE)) )
}
# remove blacklist
snv_maf_for_plotting <- snv_maf_for_plotting %>% filter(!condensed_format %in% blacklist_snvs)
# rename RA21_17hpo to RA21_17
# snv_maf_for_plotting$patient_name <- gsub("RA21_17hpo", "RA21_17", snv_maf_for_plotting$patient_name)

# rename columns alteration_class -> snv_class
snv_maf_for_plotting <- snv_maf_for_plotting %>% rename(alteration_class = snv_class)
# rename alteration classes
snv_maf_for_plotting$alteration_class <- str_replace_all(snv_maf_for_plotting$alteration_class, c("stop_gained" = "nonsense", "missense_variant" = "missense", "frameshift_truncation" = "fs_indel", "frameshift_elongation" = "fs_indel"))

# ----- 2. per-case CNV MAF -----
combined_cnv_maf <- read.xlsx("manual_stats/manual_scMAF_homdel_amps.xlsx",sheet = 1, colNames = TRUE)
# remove certain genes
blacklist_genes <- c("KRAS")

combined_cnv_maf <- combined_cnv_maf %>% filter(!gene_name %in% blacklist_genes) %>% rename("alteration_class" = "cnv_class", "fraction_altered"= "fraction_of_gene_affected") %>% filter(patient_name %in% unique(patient_info_table$patient_name)) %>% rename(CCF = CCF_manual)

# ----- 3. combined SNV and CNV MAF -----
columns_of_interest = c("patient_name", "gene_name", "alteration_class", "CCF")
composite_maf <- rbind(
  subset(snv_maf_for_plotting, select = columns_of_interest),
  subset(combined_cnv_maf, select = columns_of_interest)
)

# ----- 5. plot -----
# sort the mutations by CCF
composite_maf <- composite_maf[order(composite_maf$CCF),]

genes_to_remove <- c("ARID1B", "PALB2", "PIK3CA")
composite_maf <- composite_maf %>% filter(!gene_name %in% genes_to_remove)

# manually order genes
manual_gene_order <- rev(c("KRAS", "TP53","CDKN2A", "SMAD4", "SMAD2", "SMAD3", "TGFBR1", "TGFBR2", "ACVR1B", "BMPR1A", "ARID1A", "ARID2", "BRCA2", "ATM", "BAP1", "FGFR1","FBXW7", "RNF43", "POLD1", "IRF6", "GATA6", "MYC", "MTOR"))

# which of unique gene_name is not in manual_gene_order?
# setdiff(unique(composite_maf$gene_name), manual_gene_order)
composite_maf <- composite_maf %>% filter(gene_name %in% manual_gene_order)
composite_maf$gene_name <- factor(composite_maf$gene_name, levels = manual_gene_order)

# create grid of patient (x-axis) by gene names (y-axis)
alteration_order <- c("fs_indel", "missense", "nonsense", "splice_site_variant", "amp", "homdel")

manual_scale_colors <- setNames(
  c("#b55e02", "#119e02", "#99082e", "#00ebd7", "#f542e9", "#000000"),
  c("fs_indel", "missense", "nonsense", "splice_site_variant", "amp", "homdel")
)

# from the composite_maf, order the patient names by number of unqiue genes altered
ordered_case_names <- composite_maf %>% group_by(patient_name) %>% summarise(num_genes_altered = n_distinct(gene_name)) %>% arrange(desc(num_genes_altered))
ordered_case_names <- unique(ordered_case_names$patient_name)

# order the patient names by number of unqiue genes altered
composite_maf$patient_name <- factor(composite_maf$patient_name, levels = ordered_case_names)
patient_info_table$patient_name <- factor(patient_info_table$patient_name, levels = ordered_case_names)

fig_a <- ggplot(composite_maf, aes(
  x=patient_name, y=gene_name, 
  # size=as.factor(cut(CCF, breaks=c(0, 0.25, 0.5, 0.75, 1))),
  )) + 
  geom_point(
    shape = 22,
    data=composite_maf, 
    aes(
      fill=alteration_class, # height=CCF, width=CCF, 
      size = CCF, stroke = 0.5,
      # bin the CCF into 4 groups
      # size=cut(CCF, breaks=c(0, 0.25, 0.5, 0.75, 1)),
      # alpha=fraction_altered,
      ), 
    ) + 
  # per: https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
  # scale_size_manual(
  #   name = "CCF", 
  #   values = c(0.25, 0.5, 0.75, 1),
  #   guide = guide_legend(
  #     title.position = "top", byrow = TRUE, order=0,
  #     override.aes = list(size = c(3, 6, 9, 12),alpha = 1)
  #     )
  #   ) + 
  scale_fill_manual(
    aesthetics = "fill", 
    values = manual_scale_colors, 
    labels = alteration_order[1:4], 
    breaks = names(manual_scale_colors)[1:4], 
    name = "SNV", 
    guide = guide_legend(
      override.aes = list(size = 10),
      title.position = "top", byrow = TRUE, order=1)
    ) +
  new_scale_fill() +
  geom_point(
    shape = 22,
    data=composite_maf, 
    aes(
      fill=alteration_class, # height=CCF, width=CCF, 
      size = CCF, stroke = 0.5,
      # alpha=fraction_altered,
      ), 
    # colour="grey"
    ) + 
  scale_fill_manual(
    aesthetics = "fill", 
    values = manual_scale_colors, 
    labels = alteration_order[5:6], 
    breaks = names(manual_scale_colors)[5:6], 
    name = "CNV", 
    guide = guide_legend(
      override.aes = list(size = 10),
      title.position = "top", byrow = TRUE, order=2)
    ) +
  theme(panel.background = element_blank()) +
  geom_hline(yintercept = seq(1, nlevels(composite_maf$patient_name)) - 0.5, alpha = 0.3) +
  geom_vline(xintercept = seq(1, nlevels(composite_maf$gene_name)) - 0.5, alpha = 0.3) + 
  # xlab("Case_ID") + 
  ylab("Gene") +
  theme(
    # aspect.ratio = 1,
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold.italic"),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(), axis.line = element_blank(),
    axis.ticks = element_blank(), axis.title.x = element_blank(),
    legend.position = "right", 
    legend.direction = "vertical",
    legend.key = element_rect(fill = "white"),
    legend.spacing.y = unit(2, "mm")
  ) + 
  scale_size(range = c(0,10))
  # discrete x scale, with ticks at each patient - 0.5
  # scale_x_discrete(breaks = seq(1, nlevels(composite_maf$patient_name), by = 1) - 0.5) +
  # scale_y_discrete(breaks = seq(1, nlevels(composite_maf$gene_name), by = 1) - 0.5)

# # plot diag.stage
# fig_b <- ggplot(data = patient_info_table) + 
#   geom_tile(aes(x=patient_name, y=1, fill=diag.stage, width=0.9)) +
#   theme(panel.background = element_blank()) + 
#   scale_fill_manual(values=c("#C7E9B4","#41B6C4","#225EA8")) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
#     # legend.position = "top", 
#     legend.justification='top',
#     legend.direction = "vertical",
#     legend.key = element_rect(fill = "white"),
#   )

# plot collection.method
fig_c <- ggplot(data = patient_info_table) + 
  geom_tile(aes(x=patient_name, y=1, fill=collection_method, width=0.9)) +
  theme(panel.background = element_blank()) + 
  scale_fill_manual(values=c("#34ebdb","#3455eb","#eb6434", "#ffe75d", "#c800ff")) +
  # xlab("Case_ID") + 
  # add x axis ticks (patient_name)
  # scale_x_discrete(breaks = patient_info_table$patient_name) +
  theme(
    axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    legend.position = "right", 
    legend.justification='top',
    legend.direction = "vertical",
    legend.key = element_rect(fill = "white"),
  )

# # plot overall survival
# fig_d <- ggplot(data = patient_info_table, aes(x=patient_name, y=os)) + 
#   geom_bar(stat = "identity", width=0.9) +
#   # scale_fill_gradient(low="white", high="blue") +
#   # xlab("Case_ID") + 
#   ylab("Overall survival\n(months)") +
#   theme(
#     axis.title.y = element_text(size = 8),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(), 
#     axis.title.x = element_blank(), axis.text.x = element_blank(),
#     legend.position = "right", 
#     legend.justification='top',
#     legend.direction = "vertical",
#     legend.key = element_rect(fill = "white"),
#   )

# # for each gene in snv_maf_for_plotting, what's the mean_sc_AF_in_mut_filtered? produce the result and sort
# snv_maf_for_plotting %>% group_by(gene_name) %>% summarise(mean_sc_AF_in_mut_filtered = mean(mean_sc_AF_in_mut_filtered)) %>% arrange(desc(mean_sc_AF_in_mut_filtered))

# # save the final composite_maf
# write.csv(composite_maf, "B-oncoprint/composite_maf.FINAL.csv", row.names = FALSE)


# # ========== patchwork method ==========
combined.plot <- fig_a / fig_c + 
  plot_layout(ncol = 1, heights = c(20, 1), guides = "collect") + 
  theme(
    legend.justification = c(0.5, 0),
    )

pdf("Figure_1B_sc_oncoprint.pdf", height=9, width=11, useDingbats = FALSE)
grid.draw(combined.plot)
dev.off()

# also save to PNG with dpi=400
png("Figure_1B_sc_oncoprint.png", height=9, width=11, units = "in", res = 400)
grid.draw(combined.plot)
dev.off()

# # ----- get stats on SMAD4 and CDKN2A genes -----
# cdkn2a_stats <- composite_maf %>% filter(gene_name %in% c("CDKN2A")) %>% group_by(gene_name, alteration_class, CCF, fraction_altered) %>% summarise(freq = n())

# smad4_stats <- composite_maf %>% filter(gene_name %in% c("SMAD4")) %>% group_by(gene_name, alteration_class, CCF, fraction_altered) %>% summarise(freq = n())

# # group smad4 homdel by subclonal/sub-gene
# smad4_stats <- smad4_stats %>% mutate(alteration_class = ifelse(alteration_class == "homdel" & fraction_altered<1, "subclonal/subgene-deletion", alteration_class))
# smad4_stats_count <- smad4_stats %>% group_by(alteration_class) %>% summarise(freq = n())
# # fill in WT
# smad4_stats_count <- rbind(smad4_stats_count, data.frame(alteration_class = "WT", freq = 5))

# # group cdkn2a homdel by subclonal/sub-gene
# cdkn2a_stats <- cdkn2a_stats %>% mutate(alteration_class = ifelse(alteration_class == "homdel" & fraction_altered<1, "subclonal/subgene-deletion", alteration_class))
# cdkn2a_stats_count <- cdkn2a_stats %>% group_by(alteration_class) %>% summarise(freq = n())
# # fill in WT
# cdkn2a_stats_count <- rbind(cdkn2a_stats_count, data.frame(alteration_class = "WT", freq = 5))



