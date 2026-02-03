# Required packages are Seurat, tidyverse, dplyr and harmony

# Create a list of all seurat objects and add relevant metadata
create_seurat_list <- function(sample_info, main_dir) {
  sample_info %>%
    # each row stores a Seurat object
    mutate(
      so = map(file, ~ readRDS(file.path(main_dir, .x))) # map each row of file into .x
    ) %>%
    mutate(
      # perform parallel map to iterate row-wise over columns in list 
      so = pmap(
        list(so, sample, sample_id, stim, batch, slide),
        function(so, sample, sample_id, stim, batch, slide) {
          
          # check Seurat object version and update if needed
          current_obj_ver <- tryCatch({
            SeuratObject::Version(so)
          }, error = function(e) NA)
          installed_ver <- packageVersion("Seurat")
          if (is.na(current_obj_ver) || current_obj_ver < installed_ver) {
            message("Updating Seurat object to match Seurat package version ", installed_ver)
            so <- UpdateSeuratObject(so)
          }
          
          # update FOV's name if it's not equivalent to sample name
          sample_name <- paste0(stim, sample)
          fov_name <- names(so@images)
          if (length(fov_name) == 1 & fov_name != sample_name) {
            names(so@images) <- sample_name
          }
          
          # add metadata
          so <- AddMetaData(so, metadata = sample,    col.name = "sample")
          so <- AddMetaData(so, metadata = sample_id, col.name = "sample_id")
          so <- AddMetaData(so, metadata = stim,      col.name = "stim")
          so <- AddMetaData(so, metadata = batch,     col.name = "batch")
          so <- AddMetaData(so, metadata = slide,     col.name = "slide")
          so
        }
      )
    ) %>%
    # pull out the list of Seurat objects
    pull(so)
}

# Per-sample preprocessing and merging
preprocessing <- function(seurat_list, assay, cutoff) {
  # per‐sample QC → Normalize → Scale → Banksy
  processed_list <- map(seurat_list, function(so){
    DefaultAssay(so) <- assay
    
    # QC filter based on median 10th percentile of gene counts across 9 samples
    if (assay == "Vizgen") {
      so <- subset(so, subset = nCount_Vizgen > cutoff)
    } else if (assay == "Xenium") {
      so <- subset(so, subset = nCount_Xenium > cutoff)
    } else {
      return("Assay type must be either Vizgen or Xenium")
    }
    
    # SCTransform normalisation
    so <- SCTransform(so, assay = assay, clip.range = c(-10, 10))
    
    so
  })
  
  # name your list by sample_id so merge() can use those prefixes
  names(processed_list) <- sample_info$sample_id
  
  # now merge using those names
  project_name = ifelse(assay == "Vizgen", "vizgen_ctrl_wt_ko", "xenium_ctrl_wt_ko")
  merged_seurat <- merge(
    x           = processed_list[[1]],
    y           = processed_list[-1],
    add.cell.ids= names(processed_list),
    project     = project_name
  )
  
  return(merged_seurat)
}

# Use harmony to remove inter-sample variances and perform dimensional reduction and clustering
run_harmony_clustering <- function(merged_seurat) {
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
  
  # Prepare for differential expression
  merged_seurat <- PrepSCTFindMarkers(merged_seurat, assay = "SCT")
  
  return(merged_seurat)
}



