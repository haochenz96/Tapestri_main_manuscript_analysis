library(dplyr)
library(ggplot2)
hr_results_dir <- "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/bulk/brca/hrdetect_results"
hr_results_files <- list.files(hr_results_dir, pattern = "*hrdetect.tsv", full.names = TRUE)
hr_results <- lapply(hr_results_files, read.delim)
hr_results <- bind_rows(hr_results)

# sample would be sample split by "__" and take the first item
hr_results$sample <- sapply(strsplit(hr_results$sample, "__"), function(x) x[1])



# %%
# plot the HRdetect score per sample
# Plot the HRdetect score per sample
ggplot(hr_results, aes(x = sample, y = Probability)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "HRdetect Score per Sample", x = "Sample", y = "HRdetect Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

