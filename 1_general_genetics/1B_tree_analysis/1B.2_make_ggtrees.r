# BiocManager::install("YuLab-SMU/treedataverse")
# BiocManager::install("phytools")
library(treedataverse)
set.seed(42)
# Load the necessary libraries
library(ggtree)
library(phytools)
library(dplyr)
library(ggrepel)
library(stringr)
library(patchwork)

# change directory into the same folder as the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in NHX tree
NHX_tree_dir <- "_NHX_trees"
OUTPUT_DIR <- "."

NONFUNC_SO <- c('2kb_upstream_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'intron_variant', 'synonymous_variant')

# Create a function to read in the NHX tree file and process it for plotting
process_nhx_tree <- function(nhx_tree_f) {
  nhx_tree <- read.nhx(nhx_tree_f)
  tree_tbl <- as_tibble(nhx_tree)

  # ========== deal with edge cases ==========
  # @TEMPORARY FIX FOR TP12:
  # if any line has no parent and no node, then add them
  if (any(is.na(tree_tbl$parent)) & any(is.na(tree_tbl$node))) {
    pseudo_root <- tree_tbl$node[tree_tbl$parent == tree_tbl$node]
    pseudo_root <- pseudo_root[!is.na(pseudo_root)]
    new_node_name = max(tree_tbl$node, na.rm=TRUE) + 1
    tree_tbl$parent[is.na(tree_tbl$parent)] <- new_node_name
    tree_tbl$node[is.na(tree_tbl$node)] <- new_node_name
    tree_tbl$parent[tree_tbl$node == pseudo_root] <- new_node_name
  }

  # see if the tree is rooted at diploid: that is the row whose `parent` and `node` are equal, but with name as NA 
  is_root <- tree_tbl$parent == tree_tbl$node & is.na(tree_tbl$name)
  # when the tree is not really rooted, we will add a root node with the name "Diploid"
  if (!any(is_root)) {
    # find the pseudo-root, the node whose parent is itself
    pseudo_root <- tree_tbl$node[tree_tbl$parent == tree_tbl$node]
    if (length(pseudo_root) > 1) {
      # raise error
      stop(">2 pseudo-roots found in the tree. Please check the NHX tree file.")
    }
    num_nodes <- nrow(tree_tbl)
    # add a root node. @HZ: for some reason, the new node should be numbered the biggest, instead of 0. Otherwise will be reorganized at as.treedata
    tree_tbl <- bind_rows(tree_tbl, tibble(parent=num_nodes+1, node=num_nodes+1, name="Diploid", label="Diploid"))
    # make the pseudo-root the child of the root
    tree_tbl$parent[tree_tbl$node == pseudo_root] <- num_nodes+1
  }
  # ========== end of edge cases ==========

  # normalize clone_size by total number of cells
  tree_tbl$clone_size <- tree_tbl$clone_size / sum(tree_tbl$clone_size, na.rm=TRUE)

  # fill in the internal nodes' (the nodes with label as `NA`) names with IN-x, where x is the number of the internal node
  tree_tbl$label[is.na(tree_tbl$label)] <- paste0("IN-", 1:sum(is.na(tree_tbl$label)))

  # filter out all the somatic SNV gain/LOH elements that have any of the NONFUNC_SO as a substring, and filter out any NA or empty string
  tree_tbl$somatic_snv_gains <- sapply(
    str_split(tree_tbl$somatic_snv_gains, "\\|"), 
    function(x) paste(x[!str_detect(x, paste(NONFUNC_SO, collapse = "|"))], collapse = "|")
  )
  # convert empty strings to NA
  tree_tbl[tree_tbl == ""] <- NA
  
  # count the number of somatic SNV gains, ignore NAs
  tree_tbl$n_somatic_snv_gains <- sapply(str_split(tree_tbl$somatic_snv_gains, "\\|"), function(x) sum(x != "NA"))
  tree_tbl$n_somatic_snv_gains[is.na(tree_tbl$n_somatic_snv_gains)] <- 0

  # count the number of somatic SNV LOHs, assign 0 to NA values
  tree_tbl$n_somatic_snv_lohs <- sapply(str_split(tree_tbl$somatic_snv_lohs, "\\|"), function(x) sum(x != "NA"))
  tree_tbl$n_somatic_snv_lohs[is.na(tree_tbl$n_somatic_snv_lohs)] <- 0

  # create triangles to represent the number of somatic SNV gains and LOHs, ignore NAs # type out a triangle here --> ▲
  tree_tbl$n_somatic_snv_gains_geom <- sapply(tree_tbl$n_somatic_snv_gains, function(x) paste0(rep("▲", x), collapse = ""))
  tree_tbl$n_somatic_snv_lohs_geom <- sapply(tree_tbl$n_somatic_snv_lohs, function(x) paste0(rep("▲", x), collapse = ""))

  # Split the labels in somatic_snv_gains by the "|" character and concatenate with a newline character
  tree_tbl$somatic_snv_gains <- sapply(str_split(tree_tbl$somatic_snv_gains, "\\|"), function(x) paste(x, collapse = "\n"))
  # substitute the NA values in somatic_snv_gains with an empty string
  tree_tbl$somatic_snv_gains[tree_tbl$somatic_snv_gains == "NA"]  <- ""
  tree_tbl$somatic_snv_lohs <- sapply(str_split(tree_tbl$somatic_snv_lohs, "\\|"), function(x) paste(x, collapse = "\n"))
  tree_tbl$somatic_snv_lohs[tree_tbl$somatic_snv_lohs == "NA"]  <- ""

  tree_tbl$n_germline_snp_lohs[is.na(tree_tbl$n_germline_snp_lohs)] <- 0

  # Create n squares to represent the number of germline_snp_lohs, and store in the n_germline_snp_lohs_geom column. Create a new line every 3 squares
  tree_tbl$n_germline_snp_lohs_geom <- sapply(tree_tbl$n_germline_snp_lohs, function(x) paste0(rep("■", x), collapse = ""))
  nhx_tree <- as.treedata(tree_tbl)

  return(nhx_tree)
}

# +++++ DEBUG M12/TP12 +++++
nhx_tree_f <- file.path(NHX_tree_dir, "TP12_HZ_ETE_tree.nhx")
nhx_tree <- read.nhx(nhx_tree_f)
tree_tbl <- as_tibble(nhx_tree)

tree <- process_nhx_tree(nhx_tree_f)
d <- as_tibble(tree)

# plot a simple tree
ggtree(tree, layout = "rectangular", ladderize = FALSE) +
  geom_label(aes(label=label)) + 
  geom_point(aes(size=clone_size), color="#7e7e7e") + guides(colour=FALSE, size=FALSE)

# Plot the tree with the added tags and annotations
ggtree(tree, ladderize = FALSE, layout = "roundrect", branch.length = "dist") +
  geom_tippoint(aes(size=clone_size), color="#7e7e7e") + 
  geom_nodepoint(aes(size=clone_size), color="#7e7e7e") + guides(colour=FALSE, size=FALSE) + 
  # geom_tiplab(size=5, color="black", hjust = -2) + 
  # add in leaf size: if not NA, size = clone_size; if NA, size = 5 (default)
  # only color the tips 
  # geom_point(aes(size = clone_size, color=isTip)) +
  theme_tree2() + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_gains_geom, paste0("(.{8})"), "\\1\n")), 
    color="#ff0000", vjust=0, size=10, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#9500ffff", vjust=0.6, size=10, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_germline_snp_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#008000", vjust=1.2, size=8, label.size=NA, fill=NA) + 
  theme_tree2() +
  scale_size_area(max_size=20) + 
  coord_fixed(ratio=50)
# +++++ END OF DEBUG +++++

# ===== Plot all trees =====
# read in all the tree files
tree_files <- list.files(NHX_tree_dir, pattern = "*.nhx", full.names = TRUE)
# exclude TP6
tree_files <- tree_files[!str_detect(tree_files, "TP6")]

# create a list to store the plots
tree_objs <- list()
# plot each patient's tree, store in the list
# also calculate the max depth of eahc tree
tree_max_depths <- list()
for (tree_f in tree_files) { 
  sample_name <- strsplit(basename(tree_f), "_HZ_ETE")[[1]][1]
  tree <- process_nhx_tree(tree_f)
  tree_objs[[sample_name]] <- tree
  # depth first search to calculate the max depth for each tree
  tree_tbl = as_tibble(tree)
  max_depth <- 0
  for (node in tree_tbl$node) {
    depth <- 0
    # note that the root is defined as the node whose parent is itself
    while (node != tree_tbl$parent[tree_tbl$node == node]) {
      depth <- depth + tree_tbl$dist[tree_tbl$node == node]
      node <- tree_tbl$parent[tree_tbl$node == node]
    }
    if (depth > max_depth) {
      max_depth <- depth
    }
  }
  tree_max_depths[[sample_name]] <- max_depth
}
# order the trees by the total tree length
tree_objs <- tree_objs[order(unlist(tree_max_depths), decreasing=FALSE)]
# 
class(tree_objs) = "multiPhylo"
# add yscale which is the number of clones

ggtree(tree_objs, layout = "roundrect", branch.length = "dist", ladderize = FALSE, ) +
  facet_wrap(~.id, scales="free_y", nrow=24, strip.position = "right") + 
  geom_tippoint(aes(size=clone_size*100), color="#7e7e7e") + 
  geom_nodepoint(aes(size=clone_size*100), color="#7e7e7e") + guides(colour=FALSE, size=FALSE) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_gains_geom, paste0("(.{8})"), "\\1\n")), 
    color="#ff0000", vjust=0, size=5, label.size=NA, fill=NA, lineheight = .5) + # need to reduce the space between lines
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#9500ffff", vjust=0.5, size=5, label.size=NA, fill=NA, lineheight = .5) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_germline_snp_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#008000", vjust=1.2, size=5, label.size=NA, fill=NA, lineheight = .5) + 
  theme_tree2() + 
  scale_size_area(max_size=10) + 
  scale_x_reverse() + 
  # scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) + 
  # layout_dendrogram() + # to make the trees vertical
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    # disable x axis lines
    axis.line.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0, "lines")
    ) -> all_trees_plot


# # save to PDF # @HZ: have issue! might have to use CairoPDF
# ggsave(
#   paste0(OUTPUT_DIR, "/all_trees_plot.pdf"),
#   all_trees_plot, 
#   device=pdf, 
#   # family="Arial Unicode MS",
#   width = 15, height = 20)

# save a PNG with dpi=400
ggsave(
  paste0(OUTPUT_DIR, "/all_trees_plot.png"),
  all_trees_plot, 
  width = 15, height = 24, dpi=400)

# ================= save a table 


