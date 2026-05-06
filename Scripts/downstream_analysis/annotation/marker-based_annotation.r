suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  # library(HGNChelper)
  library(Seurat)
})

setwd("/vast/scratch/users/zhang.ji/data")
source("~/Scripts/functions/ScType_annotation.r")

# # Load merged Merscope object
# merged_vizgen <- readRDS("default_segmentation/Merscope_seurat/vizgen_harmony_ctrl_wt_ko.rds")
# merged_proseg <- readRDS("proseg_segmentation/Merscope_seurat/vizgen_harmony_ctrl_wt_ko.rds")
# merged_cellpose <- readRDS("cellpose_segmentation/Merscope_seurat/vizgen_harmony_ctrl_wt_ko.rds")

# # Load merged Xenium unimodal object
# merged_default <- readRDS("default_segmentation/Xenium_seurat/unimodal/xenium_harmony_ctrl_wt_ko.rds")
# merged_proseg <- readRDS("proseg_segmentation/Xenium_seurat/unimodal/xenium_harmony_ctrl_wt_ko.rds")

# load merged Xenium multimodal object
# merged_default <- readRDS("default_segmentation/Xenium_seurat/multimodal/xenium_harmony_ctrl_wt_ko.rds")
merged_proseg <- readRDS("proseg_segmentation/Xenium_seurat/multimodal/xenium_harmony_ctrl_wt_ko.rds")
# merged_cellpose <- readRDS("cellpose_segmentation/Xenium_seurat/xenium_harmony_ctrl_wt_ko.rds")

# Set custom markers for annotation
marker_list <- list(
  `Macrophages` = c("Adgre1","Cd209b","Cd274","Cd68","Cd80","Csf1r"),
  `B cells` = c("Cd19","Cd22","Ighd"),
  `GC B cells` = c("Aicda","Bcl6","Rgs13"),
  `Neutrophils` = c("Ngp","S100a9"),
  `Plasma cells` = c("Cd38","Jchain"),
  `T cells` = c("Cd3d","Cd3e","Cd4","Cd8a","Trac"),
  `NK cells` = c("Ncr1","Gzma"),
  `Dendritic cells` = c("Xcr1","Siglech","Spib", "Ffar2", "Cox6a2"),
  `Endothelial cells` = c("Egfl7","Madcam1"),
  `Fibroblastic reticular cells` = c("Ccl19","Dpt"),
  `Erythrocyte-like` = c("Tpx2", "Rrm1", "Ezh2", "Cdca8", "Pola1", "Ccnd3") # NB: the gene panel does not have strong markers for erythrocytes, those are found from single cell reference subset to MERSCOPE gene panel
  # `GC B cells (DZ)` = c("Rrm1","Cdca8","Gpsm2","Nek2","Pfn2",
  #                       "Lmo4","Akap12","Otub2","Scn8a","Tifa"),
  # `GC B cells (LZ)` = c("Actb","Cd38","Stx11","Cd83","Cd86",
  #                       "Cd40","B3gnt5","Fcer2a","Ankrd33b")
)

# Define the cluster resolution level for annotation
res <- 0.6

# # Assign cluster label to Merscope data
# merged_vizgen <- AssignCluster(merged_vizgen, res, marker_list)
# merged_proseg <- AssignCluster(merged_proseg, res, marker_list)
# merged_cellpose <- AssignCluster(merged_cellpose, res, marker_list)

# # Assign cluster label to Xenium unimodal data
# merged_default <- AssignCluster(merged_default, res, marker_list)
# merged_proseg <- AssignCluster(merged_proseg, res, marker_list)

# # Assign cluster label to Xenium multimodal data
# merged_default <- AssignCluster(merged_default, res, marker_list)
merged_proseg <- AssignCluster(merged_proseg, res, marker_list)
# merged_cellpose <- AssignCluster(merged_cellpose, res, marker_list)

# Save annotated results
# saveRDS(merged_vizgen, file = "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds")
# saveRDS(merged_proseg, file = "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds")
# saveRDS(merged_cellpose, file = "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds")
# saveRDS(merged_default, file = "default_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds")
# saveRDS(merged_proseg, file = "proseg_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds")
# saveRDS(merged_default, file = "default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds")
saveRDS(merged_proseg, file = "proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds")
# saveRDS(merged_cellpose, file = "cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds")

