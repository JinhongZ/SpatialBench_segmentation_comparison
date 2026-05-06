suppressPackageStartupMessages({
  library(tibble)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(harmony)
})
options(future.globals.maxSize = 8000  * 1024^3)

source("~/Scripts/functions/processing.r")
setwd("/vast/scratch/users/zhang.ji/data/default_segmentation/Xenium_seurat/unimodal")
main_dir <- "/vast/projects/SpatialBench/analysis/seurat_xenium/xenium_v2010"

# define files and their metadata
sample_info <- tribble(
  ~file,        ~sample, ~stim,  ~batch,    ~slide,
  "ko166.rds",   "166",   "ko",   "batch24", "slide0011456",
  "ko167.rds",   "167",   "ko",   "batch24", "slide0011456",
  "ko168b.rds",  "168",   "ko",   "batch27", "slide0017329",
  "ctrl172b.rds","172",   "ctrl", "batch27", "slide0017329",
  "ctrl173.rds", "173",   "ctrl", "batch24", "slide0011456",
  "ctrl174b.rds","174",   "ctrl", "batch27", "slide0017329",
  "wt709b.rds",  "709",   "wt",   "batch27", "slide0017329",
  "wt710.rds",   "710",   "wt",   "batch27", "slide0017329",
  "wt713.rds",   "713",   "wt",   "batch24", "slide0011456"
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
