suppressPackageStartupMessages({
  library(tibble)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(harmony)
})
options(future.globals.maxSize = 8000  * 1024^3)

source("~/Scripts/functions/processing.r")
setwd("/vast/scratch/users/zhang.ji/data/proseg_segmentation/Xenium_seurat/multimodal")

main_dir <- "/vast/projects/SpatialBench/segmentation/Proseg/updated_seurat_objects/Xenium"
batch34_dir <- list.files(main_dir, pattern = "^Batch34__0032118")
sample_info <- tibble(file = batch34_dir) %>% 
  mutate(
    sample_id = case_when(
      grepl("Region_1", file) ~ "ctrl_173",
      grepl("Region_2", file) ~ "wt_710",
      grepl("Region_3", file) ~ "wt_709",
      grepl("Region_4", file) ~ "ctrl_174",
      grepl("Region_5", file) ~ "ctrl_172",
      grepl("Region_6", file) ~ "embryo"
    ),
    batch = "batch34",
    slide = "slide0032118"
  ) %>% 
  filter(sample_id != "embryo") %>% 
  tidyr::separate(col = sample_id, into = c("stim", "sample"), sep = "_", remove = FALSE)

# create a list of seurat objects
seurat_list <- create_seurat_list(sample_info, main_dir)

# per-sample quality control (based on median 10th percentile of gene counts across 9 samples), normalisation, and merging
merged_seurat <- preprocessing(seurat_list, assay = "Xenium", cutoff = 10)
saveRDS(merged_seurat, file = "xenium_ctrl_wt_ko.rds")

# harmony integration to remove inter-sample variances, dimensional reduction, and clustering
merged_seurat <- run_harmony_clustering(merged_seurat)

# save harmonised seurat object
saveRDS(merged_seurat, file = "xenium_harmony_ctrl_wt_ko.rds")
