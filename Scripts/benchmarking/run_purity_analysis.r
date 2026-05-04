suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(dplyr)
  library(tibble)
})

setwd("/vast/scratch/users/zhang.ji/data")
source("~/SpatialBench_segmentation_comparison/Scripts/functions/purity_analysis.r")

data_info <- tribble(
  ~file,                                                                           ~segmentation, ~platform,           ~assay,
  "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",          "Default",     "MERSCOPE",          "Vizgen",
  "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",           "Proseg",      "MERSCOPE",          "Vizgen",
  "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",         "Cellpose2",   "MERSCOPE",          "Vizgen",
  "default_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",   "Default",     "Xenium unimodal",   "Xenium",
  "proseg_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",    "Proseg",      "Xenium unimodal",   "Xenium",
  "default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds", "Default",     "Xenium multimodal", "Xenium", 
  "proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",  "Proseg",      "Xenium multimodal", "Xenium",
  "cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds",           "Cellpose2",   "Xenium multimodal", "Xenium"
)

# run purity analysis for spatial data
MECR_sp <- list()
for (i in 1:nrow(data_info)) {
  file_dir <- data_info$file[i]
  seg <- data_info$segmentation[i]
  plat <- data_info$platform[i]
  assay <- data_info$assay[i]
  
  obj <- readRDS(file_dir)
  MECR_sp[[paste(plat, seg, sep = "_")]] <- getMECR_panel(obj, assay_use = assay, layer_use = "counts")
}

# run purity analysis for single cell reference
sc <- readRDS("~/data/sce.labelled.fine.RDS")
MECR_sc <- getMECR_panel(obj = sc, layer_use = "counts")

saveRDS(MECR_sp, file = "~/SpatialBench_segmentation_comparison/Data/MECR_sp_customised.rds")
saveRDS(MECR_sc, file = "~/SpatialBench_segmentation_comparison/Data/MECR_sc.rds")
