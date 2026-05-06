suppressPackageStartupMessages({
  library(tibble)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(harmony)
})
options(future.globals.maxSize = 8000  * 1024^3)

source("~/Scripts/functions/processing.r")

main_dir <- "/vast/projects/SpatialBench/segmentation/Proseg/updated_seurat_objects/Merscope"
obj_dir <- list.files(main_dir, pattern = "rds$")
obj_name <- c("ctrl_172c", "ctrl_174", "ko_166",
              "ctrl_173", "wt_710b", "ko_168",
              "wt_713b", "ctrl_172", "ko_167",
              "wt_709", "wt_710", "wt_713")
names(obj_dir) <- obj_name

# define files and their metadata
sample_info <- tribble(
  ~sample, ~stim,  ~batch,    ~slide,
  "166",   "ko",   "batch22", "slide1",
  "167",   "ko",   "batch27", "slide1_2",
  "168",   "ko",   "batch26", "slide3_4",
  "172",  "ctrl", "batch27", "slide1_2",
  "173",   "ctrl", "batch26", "slide3_4",
  "174",  "ctrl", "batch22", "slide1",
  "709",   "wt",   "batch33", "slide1",
  "710",   "wt",   "batch33", "slide1",
  "713",   "wt",   "batch33", "slide1"
) %>%
  mutate(
    sample_id = paste(stim, sample, sep = "_")
  )
sample_info$file <- obj_dir[sample_info$sample_id]

# create a list of seurat objects
seurat_list <- create_seurat_list(sample_info, main_dir)

# per‐sample QC
processed_list <- map(seurat_list, function(so){
  DefaultAssay(so) <- "Vizgen"
  
  # QC filter based on median 10th percentile of gene counts across 9 samples
  if (names(so@images) == "wt709") {
    so <- subset(so, subset = nCount_Vizgen > 3)
  } else {
    so <- subset(so, subset = nCount_Vizgen > 10)
  }
  
  # SCTransform normalisation
  so <- SCTransform(so, assay = "Vizgen", clip.range = c(-10, 10))
  
  so
})

# name your list by sample_id so merge() can use those prefixes
names(processed_list) <- sample_info$sample_id

# now merge using those names
merged_seurat <- merge(
  x           = processed_list[[1]],
  y           = processed_list[-1],
  add.cell.ids= names(processed_list),
  project     = "vizgen_ctrl_wt_ko"
)

# choose your output directory
out_dir <- "/vast/scratch/users/zhang.ji/data/proseg_segmentation/Merscope_seurat"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# save as an RDS
saveRDS(merged_seurat,
        file = file.path(out_dir, "vizgen_ctrl_wt_ko.rds"))

gc() %>% invisible()

# harmony integration to remove inter-sample variances, dimensional reduction, and clustering
merged_seurat <- run_harmony_clustering(merged_seurat)

# save harmonised seurat object
saveRDS(merged_seurat, file = file.path(out_dir, "vizgen_harmony_ctrl_wt_ko.rds"))