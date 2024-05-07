library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
# library(maftools)

# excel
library(readxl)

# ----- read and process snDNA-seq MAF -----
maf_f <- "1B_tree_analysis/1B_pan_cohort_phylogeny_scMAF.csv"
maf <- read_excel(maf_f)

# ----- plot oncoprint -----
genes_of_interest <- c("SMAD4", "SMAD2", "SMAD3", "TGFBR1", "TGFBR2", "ACVR1B", "BMPR1A", "ARID1A", "ARID2")
maf <- maf %>% filter(gene_name %in% genes_of_interest) %>% filter(case_type != "organoid")
maf$gene_name <- factor(maf$gene_name, levels = rev(genes_of_interest))
# create column "truncal_or_mrca": it would be
# if case type is RSX, look at "truncal" column: "truncal" if yes, "subtruncal" if no.
# if case_type is autopsy, look at "MRCA_of_met_clones" column; :Pre-MRCA_of_mets" if yes, "Post-MRCA_of_mets" if no.
maf <- maf %>% mutate(truncal_or_mrca = if_else(
	case_type == "RSX",
	if_else(truncal == "yes", "Truncal", "Subtruncal"),
	if_else(MRCA_of_met_clones == "yes", "Pre-MRCA_of_mets", "Post-MRCA_of_mets"),
	))


oncoplot <- ggplot(maf, aes(x = case_name, y = gene_name)) +
	geom_tile(
		aes(fill = truncal_or_mrca), 
		color = "black", lwd = 1, linetype = 1
	) + 
	scale_fill_manual(
	name = "MRCA_of_met_clones?",
	values = c(
		"Truncal" = "#1a35e3", "Subtruncal" = "#7cc4ff",
		"Pre-MRCA_of_mets" = "#E31A1C", "Post-MRCA_of_mets" = "#ffd0d1"
		)
	) +
	theme_classic() + 
	theme(
		axis.title.x = element_blank(), axis.title.y = element_blank(), 
		axis.text.x = element_text(angle = 60, vjust = 0.5, size=15),
		axis.text.y = element_text(face="bold.italic", size=15),
	) + 
	coord_fixed(ratio = 1.2) + 
	guides(fill = guide_legend(reverse = TRUE))
# adjust order of legend to be c("truncal", "subtruncal","Pre-MRCA_of_mets""Post-MRCA_of_mets")

# save 
ggsave("Publication/Tapestri_main_manuscript/Figure3/tgfb_mrca_oncoplot.pdf", oncoplot, width = 10, height = 5, units = "in", useDingbats = FALSE)
ggsave("Publication/Tapestri_main_manuscript/Figure3/tgfb_mrca_oncoplot.png", oncoplot, width = 10, height = 5, units = "in", dpi=300)

