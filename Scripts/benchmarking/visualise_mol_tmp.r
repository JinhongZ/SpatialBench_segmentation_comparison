suppressPackageStartupMessages({
  library(Seurat)
})

options(future.globals.maxSize = 8000  * 1024^3)

setwd("/vast/scratch/users/zhang.ji/data")
xenium_multi_default <- readRDS("default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds")
xenium_multi_proseg <- readRDS("proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds")
xenium_multi_cellpose <- readRDS("cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds")

roi <- Crop(xenium_multi_proseg[["wt710"]], x = c(2150, 2950), y = c(2350, 3150), coords = "tissue")
xenium_multi_proseg[["roi"]] <- roi

roi <- Crop(xenium_multi_cellpose[["wt710"]], x = c(2150, 2950), y = c(2350, 3150), coords = "tissue")
xenium_multi_cellpose[["roi"]] <- roi

roi <- Crop(xenium_multi_proseg[["wt709"]], x = c(4350, 5150), y = c(2250, 3050))
xenium_multi_proseg[["roi"]] <- roi

roi <- Crop(xenium_multi_cellpose[["wt709"]], x = c(4350, 5150), y = c(2250, 3050))
xenium_multi_cellpose[["roi"]] <- roi

# visualise neutrophil marker molecules
cols <- c("T cells" = "yellowgreen", 
          "Erythrocyte-like" = "violet", 
          "Neutrophils" = "mediumblue", 
          "B cells" = "tomato", 
          "Plasma cells" = "thistle1", 
          "GC B cells" = "yellow", 
          "Macrophages" = "skyblue2",
          "Endothelial cells" = "orange",
          "Muscle cells" = "deeppink",
          "Unknown" = "grey",
          "Dendritic cells" = "darkgreen",
          "NK cells" = "dodgerblue3",
          "Fibroblastic reticular cells" = "tan3")

Idents(xenium_multi_proseg) <- "orig.ident"
Idents(xenium_multi_cellpose) <- "orig.ident"

# Neutrophils & T cells
ImageDimPlot(xenium_multi_proseg, fov = "roi", axes = T, 
             boundaries = "segmentation", 
             cols = "polychrome",
             # group.by = "ScType_label_res.0.6",
             border.color = "white", border.size = 0.1, 
             molecules = c("Ngp", "Cd3d", "Cd3e","Cd4","Cd8a","Trac"), 
             mols.size = 0.5, alpha = 0.5)

ImageFeaturePlot(xenium_multi_proseg, fov = "roi", axes = T, 
                 boundaries = "segmentation", border.color = "white", border.size = 0.1,
                 features = c("Ngp", "Cd3d", "Cd3e","Cd4","Cd8a","Trac"))

ImageDimPlot(xenium_multi_cellpose, fov = "roi", axes = T, 
             boundaries = "segmentation", cols = "polychrome",
             border.color = "white", border.size = 0.1, 
             molecules = c("Ngp", "Cd3d", "Cd3e","Cd4","Cd8a","Trac"), 
             mols.size = 0.5, alpha = 0.5)


# Fibroblastic reticular cells & GC B cells
ImageDimPlot(xenium_multi_proseg, fov = "roi", axes = T, 
             boundaries = "segmentation", cols = "polychrome",
             # group.by = "ScType_label_res.0.6",
             border.color = "white", border.size = 0.1, 
             molecules = c("Ccl19", "Dpt", "Aicda","Bcl6","Rgs13"), 
             mols.size = 0.5, alpha = 0.5)

ImageFeaturePlot(xenium_multi_proseg, fov = "roi", axes = T, 
                 boundaries = "segmentation", border.color = "white", border.size = 0.1,
                 features = c("Ccl19", "Dpt", "Aicda","Bcl6","Rgs13"))
            