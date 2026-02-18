# Required packages are Seurat, dplyr

# marker list used for cell type annotation
marker_list_ann <- list(
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
)

# set up markers for mutually exclusive correlation rate
marker_df_spleen_panel <- data.frame(
  gene = unlist(marker_list_ann),
  cell_type = rep(names(marker_list_ann), times = lengths(marker_list_ann)),
  row.names = unlist(marker_list_ann)
)

# function slightly modified from (https://github.com/Center-for-Spatial-OMICs/SpatialQM/blob/60b9217cdab5f95c8f34a7261fe76d1bbbcefebf/R/utils_update_final.R#L1729)
getMECR_panel <- function(seu_obj,
                          assay_use = NULL,
                          layer_use = c("counts","data"),
                          marker_df = marker_df_spleen_panel,
                          max_genes = 25) {
  
  stopifnot(!is.null(seu_obj))
  
  layer_use <- match.arg(layer_use)
  
  # pick assay
  if (is.null(assay_use)) assay_use <- Seurat::DefaultAssay(seu_obj)
  
  exp <- Seurat::GetAssayData(seu_obj, assay = assay_use, layer = layer_use)
  
  # intersect markers with data
  genes <- intersect(rownames(exp), marker_df$gene)
  
  if (length(genes) < 6) {
    return(list(
      MECR = NA_real_,
      n_genes_used = length(genes),
      genes_used = genes,
      note = "Too few marker genes present in this assay/panel to estimate MECR robustly."
    ))
  }
  
  # subset and coerce
  mtx <- as.matrix(exp[genes, , drop = FALSE])
  
  # binarise the matrix
  mtx_b <- mtx > 0
  
  # compute all pairwise intersection at once, returning a gene x gene matrix and each cell indicating the number of cells where both genes are detected
  inter <- mtx_b %*% t(mtx_b) 
  
  # compute unions
  gene_sums <- rowSums(mtx_b)
  # outer function computes all pairwise operations between two vectors
  union <- outer(gene_sums, gene_sums, "+") - inter
  
  # jaccard matrix
  jaccard <- inter / union
  jaccard[union == 0] <- NA # avoid cases when neither gene is detected in any cell
  
  # build mask
  genes <- rownames(mtx)
  ct <- marker_df[genes, "cell_type"]
  
  # only interest in marker expression in different cell types
  diff_ct <- outer(ct, ct, "!=")
  
  # avoid duplicated pairs (e.g., {Adgre1, Aicda} is the same as {Aicda, Adgre1})
  lower_tri <- lower.tri(jaccard)
  
  valid_pairs <- diff_ct & lower_tri
  
  # final MECR
  idx <- which(valid_pairs, arr.ind = TRUE)
  g1 <- rownames(jaccard)[idx[, 1]]
  g2 <- colnames(jaccard)[idx[, 2]]
  MECR_df <- data.frame(
    gene1 = g1,
    cell_type1 = marker_df[g1, "cell_type"],
    gene2 = g2,
    cell_type2 = marker_df[g2, "cell_type"],
    MECR = round(jaccard[valid_pairs], 4)
  )
  
  list(
    MECR_df = MECR_df,
    n_genes_available = length(intersect(rownames(exp), marker_df$gene)),
    n_genes_used = length(genes),
    genes_used = genes,
    n_pairs = nrow(MECR_df),
    assay_use = assay_use,
    layer_use = layer_use
  )
}

# # examples
# setwd("/vast/scratch/users/zhang.ji/data")
# merscope_default <- readRDS("default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds")
# merscope_default_mecr <- getMECR_panel(merscope_default, assay_use = "Vizgen", layer_use = "counts")
# merscope_default_mecr$MECR_df %>% arrange(desc(MECR))



