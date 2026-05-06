suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(purrr)
  library(future)
  library(future.apply)
  library(limma)
  library(edgeR)
})
plan("multisession", workers = 4)
options(future.globals.maxSize = 80 * 1024^3)  # 80 GB

main_dir <- "/vast/scratch/users/zhang.ji/data"
source("~/Scripts/functions/pseudo_bulk_DE.r")

data_info <- tribble(
  ~file,                                                                           ~segmentation, ~platform, 
  "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",          "Default",     "MERSCOPE",
  "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",           "Proseg",      "MERSCOPE",
  "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",         "Cellpose2",   "MERSCOPE",
  "default_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",   "Default",     "Xenium unimodal",
  "proseg_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",    "Proseg",      "Xenium unimodal"
)

groups <- c("KOvsWT", "KOvsCTRL", "WTvsCTRL")

de_all <- future_lapply(1:nrow(data_info), function(i) {
  obj <- readRDS(file.path(main_dir, data_info$file[i]))
  cell_types <- levels(obj$ScType_label_res.0.6)
  
  purrr::map_dfr(cell_types, function(cell_type) {
    purrr::map_dfr(groups, function(grp) {
      run_pseudo_bulk_de(obj, cell_label = "ScType_label_res.0.6", cell_type = cell_type, group = grp) %>% 
        mutate(cell_type = cell_type, group = grp, segmentation = data_info$segmentation[i], platform = data_info$platform[i])
    })
  })
}, future.seed = TRUE)

de_all <- bind_rows(de_all)

write.csv(de_all, file = "~/SpatialBench_segmentation_comparison/Data/pseudo-bulk_DE_results.csv", row.names = FALSE)

