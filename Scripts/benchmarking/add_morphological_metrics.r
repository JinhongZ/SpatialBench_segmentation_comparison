suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(future.apply)
  library(sf)
  library(sp)
})

source("~/Scripts/functions/compute_morphological_metrics.r")

# laod MERSCOPE data
setwd("/vast/scratch/users/zhang.ji/data")
merscope_obj_dir <- c(
  "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",
  "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",
  "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds"
)

# load Xenium unimodal data
xenium_uni_obj_dir <- c(
  "default_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",
  "proseg_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds"
)

# load Xenium multimodal data
xenium_multi_obj_dir <- c(
  "default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",
  "proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",
  "cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds"
)

for (dir in xenium_multi_obj_dir[2]) {
  message("Processing directory ", dir)
  obj <- readRDS(dir)
  
  cell_area <- extract_cell_area(obj)
  new_metadata <- data.frame(
    cell_area = cell_area,
    aspect_ratio = compute_aspect_ratio(obj),
    log10_signal_density = log10(obj$nCount_Xenium / cell_area)
  )
  obj <- AddMetaData(obj, new_metadata)
  obj <- computeSpatialOutlier(obj, computeBy = "log10_signal_density", method = "both")
  
  message("Saving updated ", dir)
  saveRDS(obj, file = dir)
  gc() %>% invisible()
}