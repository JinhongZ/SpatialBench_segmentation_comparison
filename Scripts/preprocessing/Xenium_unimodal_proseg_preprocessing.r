suppressPackageStartupMessages({
  library(tibble)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(harmony)
})
options(future.globals.maxSize = 8000  * 1024^3)

source("~/Scripts/functions/processing.r")
setwd("/vast/scratch/users/zhang.ji/data/proseg_segmentation/Xenium_seurat/unimodal")
main_dir <- "/vast/projects/SpatialBench/segmentation/Proseg/updated_seurat_objects/Xenium"


# list file directories and their metadata
sample_info <- tribble(
  ~file,                            ~sample, ~stim,  ~batch,    ~slide,
  "Batch24__0011456__Region_1.rds", "166",   "ko",   "batch24", "slide0011456",
  "Batch24__0011456__Region_2.rds", "167",   "ko",   "batch24", "slide0011456",
  "Batch27__0017329__Region_3.rds", "168",   "ko",   "batch27", "slide0017329",
  "Batch27__0017329__Region_2.rds", "172",   "ctrl", "batch27", "slide0017329",
  "Batch24__0011456__Region_4.rds", "173",   "ctrl", "batch24", "slide0011456",
  "Batch27__0017329__Region_4.rds", "174",   "ctrl", "batch27", "slide0017329",
  "Batch27__0017329__Region_5.rds", "709",   "wt",   "batch27", "slide0017329",
  "Batch27__0017329__Region_1.rds", "710",   "wt",   "batch27", "slide0017329",
  "Batch24__0011456__Region_3.rds", "713",   "wt",   "batch24", "slide0011456"
) %>%
  mutate(
    sample_id = paste(stim, sample, sep = "_")
  )


# create a list of seurat objects
seurat_list <- create_seurat_list(sample_info, main_dir)

# per-sample quality control (based on median 10th percentile of gene counts across 9 samples), normalisation, and merging
merged_seurat <- preprocessing(seurat_list, assay = "Xenium", cutoff = 10)
saveRDS(merged_seurat, file = "xenium_ctrl_wt_ko.rds")

# harmony integration to remove inter-sample variances, dimensional reduction, and clustering
merged_seurat <- run_harmony_clustering(merged_seurat)

# save harmonised seurat object
saveRDS(merged_seurat, file = "xenium_harmony_ctrl_wt_ko.rds")
