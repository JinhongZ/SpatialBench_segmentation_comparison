suppressPackageStartupMessages({
  library(tidyverse)
  library(tidyr)
  library(purrr)
  library(data.table)
  library(dplyr)
  library(Seurat)
})

setwd("~/SpatialBench_segmentation_comparison/Data")
source("../Scripts/functions/benchmarking_functions.r")

main_dir <- "/vast/scratch/users/zhang.ji/data"
sample_info <- tribble(
  ~platform,  ~segmentation,          ~sub_dir,                                                                            ~count_col,          ~feature_col,          ~model,          
  "MERSCOPE", "Cellpose1",            "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",              "nCount_Vizgen",     "nFeature_Vizgen",     NA_character_,
  "MERSCOPE", "Proseg",               "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",               "nCount_Vizgen",     "nFeature_Vizgen",     NA_character_,
  "MERSCOPE", "Cellpose2",            "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",             "nCount_Vizgen",     "nFeature_Vizgen",     NA_character_,
  
  "Xenium",   "Nuclear expansion",    "default_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",       "nCount_Xenium",     "nFeature_Xenium",     "unimodal",
  "Xenium",   "Proseg",               "proseg_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",        "nCount_Xenium",     "nFeature_Xenium",     "unimodal",
  
  "Xenium",   "Nuclear expansion",    "default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",     "nCount_Xenium",     "nFeature_Xenium",     "multimodal",
  "Xenium",   "Proseg",               "proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",      "nCount_Xenium",     "nFeature_Xenium",     "multimodal",
  "Xenium",   "Cellpose2",            "cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds",               "nCount_Xenium",     "nFeature_Xenium",     "multimodal"
) %>% 
  mutate(
    file_path = file.path(main_dir, sub_dir)
  )

out_path = "~/SpatialBench_segmentation_comparison/Data"
merscope_obj <- readRDS(sample_info$file_path[1])
common_genes <- rownames(merscope_obj)

save_sample_and_cell_df(sample_info, common_genes = common_genes, out_path = out_path)
  
