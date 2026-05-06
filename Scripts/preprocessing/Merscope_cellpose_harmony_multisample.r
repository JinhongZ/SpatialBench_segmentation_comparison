suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})
options(future.globals.maxSize = 8000  * 1024^3)

# load the merged object
main_dir     <- "/vast/scratch/users/zhang.ji/data/cellpose_segmentation/Merscope_seurat"
merged_path  <- file.path(main_dir, "vizgen_ctrl_wt_ko.rds")
merged_seurat <- readRDS(merged_path)

# PCA 
merged_seurat <- RunPCA(
  merged_seurat,
  assay          = "SCT",
  features       = rownames(merged_seurat),
  npcs           = 30,
  verbose        = FALSE
)

# Harmony integration to correct for sample‐to‐sample variation
set.seed(42)
merged_seurat <- harmony::RunHarmony(
  merged_seurat,
  group.by.vars = "sample_id",    # or "batch"
  reduction.use     = "pca",
  reduction.save = "harmony",
  assay.use     = "SCT",
  project.dim   = FALSE
)

# Build SNN graph & cluster on the Harmony embedding
merged_seurat <- FindNeighbors(
  merged_seurat,
  reduction  = "harmony",
  dims       = 1:30,
  verbose    = FALSE
)
merged_seurat <- FindClusters(
  merged_seurat,
  resolution = 0.6,
  verbose    = FALSE
)

# UMAP for visualization
merged_seurat <- RunUMAP(
  merged_seurat,
  reduction   = "harmony",
  dims        = 1:30,
  seed.use    = 42,
  verbose     = FALSE
)

# check with quick plot
# p1 <- DimPlot(
#  merged_seurat,
#  group.by = "sample_id",
#  pt.size  = 0.3
# ) + ggtitle("UMAP: by Sample ID")
# 
# p2 <- DimPlot(
#  merged_seurat,
#  group.by = "SCT_snn_res.0.4",
#  label   = TRUE,
#  pt.size = 0.3
# ) + ggtitle("UMAP: Harmony Clusters")
# 
# print(p1)
# print(p2)

# prep for differential expression
merged_seurat <- PrepSCTFindMarkers(merged_seurat, assay = "SCT")

# save the integrated object
out_file <- file.path(main_dir, "vizgen_harmony_ctrl_wt_ko.rds")
saveRDS(merged_seurat, file = out_file)