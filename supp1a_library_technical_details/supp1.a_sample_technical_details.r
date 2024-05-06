library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)

# configure default font size, font family, theme, figure size
theme_set(theme_pubr())
options(repr.plot.width = 10, repr.plot.height = 5)

# read in master sample sheet
master_sample_sheet <- read_excel(
  "Tapestri_batch2_samples_MASTER.xlsx", 
  sheet = "all_sample_clinical_info"
)
master_sample_sheet %>% 
  filter(censor == 0) -> master_sample_sheet
# the collection_method should be ordered:
master_sample_sheet$collection_method <- factor(
  master_sample_sheet$collection_method, 
  levels = c("derived_organoid", "autopsy", "resection", "biopsy")
)
# sort
master_sample_sheet <- master_sample_sheet %>% arrange(collection_method)

##############################################################################
# ===== plot per-sample library nuclei yield, split by collection method =====
##############################################################################

# palette for collection_method:
# derived_organoid: "#eb6434"
# autopsy: "#3455eb"
# resection: "#34ebdb"
# biopsy: "#ffe75d"
collection_palette <- c(
  "derived_organoid" = "#eb6434",
  "autopsy" = "#3455eb",
  "resection" = "#34ebdb",
  "biopsy" = "#ffe75d"
)

# plot bar plot of each sample's "library nuclei yield", color by collection method
ggbarplot(
  master_sample_sheet,
  x = "sample", y = "library nuclei yield",
  fill = "collection_method",
) +
  scale_fill_manual(values = collection_palette) +
  theme_classic() +
  coord_fixed(ratio = 0.005) + 
  # disable x ticks/tick marks, set y range to be 0 through 8000
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylim(0, 8000)

# save
ggsave(
  "supp1a_library_technical_details/supp1.1_sample_technical_details.pdf", 
  width = 10, height = 5, units = "in", useDingbats = FALSE
)

#######################
# ===== get stats =====
#######################
# overall median library nuclei yield
# overall standard deviation of library nuclei yield
# overall total number of single nuclei libraries

master_sample_sheet %>% 
  summarise(
    mean_library_nuclei_yield = mean(`library nuclei yield`),
    sd_library_nuclei_yield = sd(`library nuclei yield`),
    total_single_nuclei_libraries = sum(`library nuclei yield`)
  ) -> stats
# also get total number of single nuclei per case, not sample, and take median of that
master_sample_sheet %>% 
  group_by(case) %>% 
  summarise(
    total_single_nuclei_libraries = sum(`library nuclei yield`)
  ) -> stats_by_case
mean_by_case <- mean(stats_by_case$total_single_nuclei_libraries)
mean_by_case
# median library nuclei yield for each collection method
# standard deviation of library nuclei yield for each collection method

master_sample_sheet %>% 
  group_by(collection_method) %>% 
  summarise(
    mean_library_nuclei_yield = mean(`library nuclei yield`),
    sd_library_nuclei_yield = sd(`library nuclei yield`)
  ) -> stats_by_collection_method
