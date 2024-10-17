library(ggplot2)
library(dplyr)
library(ggpubr)
# library(plotly)

output_dir <- "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/N_BRCA_analysis/hrdetect"
# ===== Plot indel proportions =====
indel_result <- "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/N_BRCA_analysis/hrdetect/WGS_indel_results.2024-07-19.tsv"
indel_props <- read.table(indel_result, header=T, sep="\t")

# Rename the samples, change into factor
indel_props$sample <- factor(recode(
  indel_props$sample,
  "M04_normal" = "PC01-normal",
  "M04-3" = "PC01-RSX",
  "M13-1" = "PC05-RSX",
  "RSX425" = "PC23-RSX",
  "RSX476" = "PC24-RSX",
  "RSX476_normal" = "PC22-normal",
  "RSX582" = "PDAC_wo_gBRCA",
  "RSX640" = "PC22-RSX",
  "RSX640_normal" = "PC22-normal",
  "RA15_6" = "PC06-autopsy",
  "RA16_8" = "PC08-autopsy",
  "RA17_13" = "PC10-autopsy",
  "RA17_22" = "PC11-autopsy",
  "BPA-1-RSX" = "PC20-RSX",
  "BPA-1-IR" = "PC20-IR",
  "BPA-2-RSX" = "PC21-RSX",
  "BPA-2-IR" = "PC21-IR",
))

# filter out the samples with "autopsy"
indel_props <- indel_props[!grepl("autopsy", indel_props$sample),]

# Add BRCA2 status column
indel_props$BRCA2_status <- factor(recode(
  indel_props$sample,
  "PC01-normal" = "intact",
  "PC01-RSX" = "biallelic inactivated",
  "PC05-RSX" = "intact",
  "PC23-RSX" = "biallelic inactivated",
  "PC22-normal" = "intact",
  "PC24-RSX" = "sBRCA2 not detected by Tapestri",
  "PDAC_wo_gBRCA" = "intact",
  "PC22-RSX" = "biallelic inactivated",
  "PC06-autopsy" = "intact",
  "PC08-autopsy" = "intact",
  "PC10-autopsy" = "intact",
  "PC11-autopsy" = "intact",
  "PC20-RSX" = "biallelic inactivated",
  "PC20-IR" = "biallelic inactivated",
  "PC21-RSX" = "biallelic inactivated",
  "PC21-IR" = "biallelic inactivated",
  "PC23-RSX" = "biallelic inactivated"
))

indel_props %>% 
  # sort by BRCA2 status, then sample_name
  arrange(BRCA2_status) -> wgs_sample_indel_props

hrdetect_plot <- ggbarplot(wgs_sample_indel_props, x="sample", y="del.mh.prop", fill="BRCA2_status", border="black", rotate = TRUE, ggtheme = theme_minimal()) +
  scale_fill_manual(values=c("intact" = "lightgreen", "biallelic inactivated" = "#ff1b3d", "sBRCA2 not detected by Tapestri" = "#ff9f39")) +
  geom_text(aes(label=round(del.mh.prop, 2)), vjust=-0.5, size=5) +
  labs(title="Proportion of indels at microhomology sites",y="Proportion") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size=24),
    )

ggsave(paste0(output_dir, "/indel_at_mh_site_props.pdf"), hrdetect_plot, width=15, height=10, units="in")
ggsave(paste0(output_dir, "/indel_at_mh_site_props.png"), hrdetect_plot, width=15, height=10, units="in", dpi=300)

# HRdetect results
hr_results_dir <- "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/N_BRCA_analysis/hrdetect/hrdetect_results"
hr_results_files <- list.files(hr_results_dir, pattern = "*hrdetect.tsv", full.names = TRUE)
hr_results <- lapply(hr_results_files, read.delim)
hr_results <- bind_rows(hr_results)

# sample would be sample split by "__" and take the first item
hr_results$sample_renamed <- sapply(strsplit(hr_results$sample, "__"), function(x) x[1])
hr_results$sample_renamed <- sapply(strsplit(hr_results$sample_renamed, "_tumor"), function(x) x[1])
# recode samples
hr_results$sample_renamed <- factor(recode(
  hr_results$sample_renamed,
  "M04-3" = "PC01-RSX",
  "M13-1" = "PC05-RSX",
  "RSX425" = "PC23-RSX",
  "RSX435_normal" = "PC01-normal",
  "RSX476" = "PC24-RSX",
  "RSX582" = "PDAC_wo_gBRCA",
  "RSX640" = "PC22-RSX",
  "RSX640_normal" = "PC22-normal",
  "BPA-1-RSX" = "PC20-RSX",
  "BPA-1-IR" = "PC20-IR",
  "BPA-2-RSX" = "PC21-RSX",
  "BPA-2-IR" = "PC21-IR",
  "BPA-4-RSX" = "PC23-RSX"
))

# Add BRCA2 status column
hr_results$BRCA2_status <- factor(recode(
  hr_results$sample_renamed,
  "PC01-normal" = "intact",
  "PC01-RSX" = "biallelic inactivated",
  "PC05-RSX" = "intact",
  "PC23-RSX" = "biallelic inactivated",
  "PC24-normal" = "intact",
  "PC24-RSX" = "sBRCA2 not detected by Tapestri",
  "PDAC_wo_gBRCA" = "intact",
  "PC22-RSX" = "biallelic inactivated",
  "PC20-RSX" = "biallelic inactivated",
  "PC20-IR" = "biallelic inactivated",
  "PC21-RSX" = "biallelic inactivated",
  "PC21-IR" = "biallelic inactivated",
  "PC23-RSX" = "biallelic inactivated"
))
# remove duplicated samples, keeping the higher HRdetect score
hr_results <- hr_results %>%
  group_by(sample_renamed) %>%
  filter(Probability == max(Probability))

# Use ggpubr to make a barplot of the HRdetect score per sample
hrdect_score_plt <- ggbarplot(hr_results, x="sample_renamed", y="Probability", fill="BRCA2_status", border="black", rotate = TRUE, ggtheme = theme_minimal()) +
  scale_fill_manual(values=c("intact" = "lightgreen", "biallelic inactivated" = "#ff1b3d", "sBRCA2 not detected by Tapestri" = "#ff9f39")) +
  geom_text(aes(label=round(Probability, 2)), vjust=-0.5, size=5) +
  labs(title="HRDetect score") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size=24),
    )
# save
ggsave(paste0(output_dir, "/hrdetect_score_per_sample.pdf"), hrdect_score_plt, width=15, height=10, units="in")
ggsave(paste0(output_dir, "/hrdetect_score_per_sample.png"), hrdect_score_plt, width=15, height=10, units="in", dpi=300)
