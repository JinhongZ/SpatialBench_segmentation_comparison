# Required packages are Seurat, Matrix, data.table, dplyr, ggplot2, scales

# find the transcripts counts for genes overlapping MERSCOPE and Xenium panels
add_common_gene_counts <- function(obj, common_genes) {
  # use dgCMatrix for faster processing
  counts <- Matrix::Matrix(obj[["Xenium"]]@layers$counts, sparse = TRUE)
  rownames(counts) <- rownames(obj)
  
  # subset the count matrix to the common gene panel
  counts <- counts[intersect(rownames(counts), common_genes), , drop = FALSE]
  
  # aggregate those gene counts
  common_gene_counts <- Matrix::colSums(counts)
  common_feature_counts <- Matrix::colSums(counts > 0)
  
  # add metadata
  obj <- Seurat::AddMetaData(
    obj, 
    metadata = data.frame(
      common_gene_counts = common_gene_counts,
      common_feature_counts = common_feature_counts
    )
  )
  
  return(obj)
} 

generate_sample_df <- function(obj, segmentation_clean, count_col) {
  # extract cell metadata
  dt <- data.table::as.data.table(obj@meta.data)
  has_common <- "common_gene_counts" %in% names(dt)
  
  # group cells by sample_id and compute summaries per sample
  # syntax dt[i, j, by]
  dt[, .(
    cell_count = .N, # number of rows in current group, faster than n()
    transcript_count = sum(get(count_col), na.rm = TRUE),
    common_gene_counts = if (has_common) sum(common_gene_counts, na.rm = TRUE) else sum(get(count_col), na.rm = TRUE),
    batch = batch[1],
    segmentation_clean = segmentation_clean
  ), by = sample_id]
}

generate_cell_df <- function(obj, segmentation_clean, count_col, feature_col) {
  # extract cell metadata
  dt <- data.table::as.data.table(obj@meta.data)
  dt[, .(
    sample = gsub("_", "", sample_id),
    nCounts = get(count_col),
    nFeatures = get(feature_col),
    common_gene_counts = if ("common_gene_counts" %in% names(dt)) common_gene_counts else get(count_col),
    common_feature_counts = if ("common_feature_counts" %in% names(dt)) common_feature_counts else get(feature_col),
    cell_type = ScType_label_res.0.6,
    cell_area = if("cell_area" %in% colnames(dt)) cell_area else NA_real_,
    aspect_ratio = if("aspect_ratio" %in% colnames(dt)) aspect_ratio else NA_real_,
    log10_signal_density = if("log10_signal_density" %in% colnames(dt)) log10_signal_density else NA_real_,
    log10_signal_density_outlier_sc = if("log10_signal_density_outlier_sc" %in% colnames(dt)) log10_signal_density_outlier_sc else NA_real_,
    segmentation_clean = segmentation_clean
  )]
}

extract_sample_and_cell_df <- function(
    file_path,
    segmentation_clean,
    count_col,
    feature_col
) {
  obj <- readRDS(file_path)
  
  if (grepl("Xenium", file_path)) {
    obj <- add_common_gene_counts(obj, common_genes)
  }
  
  sample_df <- generate_sample_df(
    obj,
    segmentation_clean = segmentation_clean,
    count_col = count_col
  )
  
  cell_df <- generate_cell_df(
    obj,
    segmentation_clean = segmentation_clean,
    count_col   = count_col,
    feature_col = feature_col
  )
  
  rm(obj)
  gc()
  
  list(
    sample_df = sample_df,
    cell_df = cell_df
  )
}

save_sample_and_cell_df <- function(sample_info, common_genes, out_path) {
  metadata <- sample_info %>% 
    mutate(
      out = purrr::pmap(
        list(file_path, segmentation_clean, count_col, feature_col),
        extract_sample_and_cell_df
      )
    ) %>% 
    tidyr::unnest_wider(out)
  
  sample_df <- metadata %>% 
    group_by(platform, model) %>% 
    group_modify(~data.table::rbindlist(.x$sample_df, fill = TRUE))
  
  cell_df <- metadata %>% 
    group_by(platform, model) %>% 
    group_modify(~{
      clean_list <- lapply(.x$cell_df, function(dt) {
        dt[, log10_signal_density_outlier_sc :=
             as.character(log10_signal_density_outlier_sc)]
        dt
      })
      data.table::rbindlist(clean_list, fill = TRUE)
    })
  
  
  write.csv(sample_df, file = file.path(out_path, "sample_df.csv"), row.names = FALSE)
  data.table::fwrite(cell_df, file = file.path(out_path, "cell_df.csv.gz"), row.names = FALSE)
}


# --- plots for sample-level quality metrics
colour_panel <- c(
  "Default" = "#F8766D",
  "Proseg" = "#00BA38",
  "Cellpose2" = "#619CFF"
)

sample_boxplot <- function(df, metric, label, jitter_width = NULL, jitter_height = NULL) {
  # set default jitter values if not provided
  if (is.null(jitter_width)) jitter_width <-  0.1
  if (is.null(jitter_height)) jitter_height <-  0.1
  
  # create jitter layer
  jitter_layer <- geom_jitter(width = jitter_width, height = jitter_height)
  
  ggplot(df, 
         aes(x = segmentation_clean, y = .data[[metric]], fill = segmentation_clean)) +
    geom_boxplot(outliers = FALSE) + 
    jitter_layer + 
    # geom_line(aes(group = sample_id), colour = "darkred", linewidth = 0.5, alpha = 0.7) + # uncomment to have lines connected by sample_id
    scale_fill_manual(values = colour_panel) +
    scale_y_continuous(labels = scales::comma) + 
    facet_wrap(~platform_version) + 
    labs(x = "Segmentation",
         y = label) + 
    # scale_y_continuous(n.breaks = 8, labels = comma) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "none"
    )
}

# --- plots for cell-level quality and morphological metrics
cell_violin_plot <- function(cell_df, metric, label) {
  ggplot(
    cell_df, 
    aes(x = sample, y = .data[[metric]], fill = segmentation_clean)
  ) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(aes(group = interaction(sample, segmentation_clean)), 
                 width = 0.15, outlier.shape = NA, fill = "white", color = "black",
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = colour_panel) +
    facet_wrap(~platform_version) + 
    labs(x = "Samples",
         y = label,
         fill = "segmentation_clean") + 
    scale_y_continuous(n.breaks = 8) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

cell_median_boxplot <- function(cell_df, metric, label, jitter_width = NULL, jitter_height = NULL) {
  median_df <- cell_df %>% 
    group_by(platform_version, segmentation_clean, sample) %>% 
    summarise(
      median_value = median(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # set default jitter values if not provided
  if (is.null(jitter_width)) jitter_width <-  0.1
  if (is.null(jitter_height)) jitter_height <-  0.1
  
  # create jitter layer
  jitter_layer <- geom_jitter(width = jitter_width, height = jitter_height)
  
  ggplot(
    median_df, 
    aes(x = segmentation_clean, y = .data[["median_value"]], fill = segmentation_clean)
  ) +
    geom_boxplot(outliers = FALSE) +
    jitter_layer +
    scale_fill_manual(values = colour_panel) +
    facet_wrap(~platform_version) + 
    labs(x = "Segmentation", 
         y = paste("Median", label)) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

# plot dimensional plots across platforms
combine_DimPlot <- function(seurat_list, group.by = "ScType_label_res.0.6", cols) {
  dim_plots <- lapply(seq_along(seurat_list), function(i) {
    p <- DimPlot(
      seurat_list[[i]], reduction = "umap", group.by = group.by, 
      cols = cols, pt.size = 0.5, raster = T
    ) + labs(title = "")
    
    # hide legend for all expect the last plot
    if (i != length(seurat_list)) {
      p <- p + NoLegend()
    }
    
    return(p)
  })
  
  # combine plots in one row
  combined_plot <- patchwork::wrap_plots(dim_plots, nrow = 1)
  
  return(combined_plot)
}

# customised function for cropping as latest Seurat (v5.4.0) does not support segmentation_clean overlaying 
# fov = Seurat FOV object, 
# x = range of x_coords after cropping, 
# y = range of y_coords after cropping
crop_fov <- function(fov, x, y) {
  stopifnot(length(x) == 2, length(y) == 2)
  
  # extract centroid coords
  coords <- fov$centroids@coords
  
  # find cells within the cropping region
  x_idx <- coords[, "x"] >= x[1] & coords[, "x"] <= x[2]
  y_idx <- coords[, "y"] >= y[1] & coords[, "y"] <= y[2]
  keep_cells <- fov$centroids@cells[x_idx & y_idx]
  
  # return the new FOV object manually
  return(subset(fov, cells = keep_cells))
}

# combine multiple image dimensional plots across platforms
combine_ImageDimPlot <- function(
    seurat_list, 
    fov, 
    group.by = "ScType_label_res.0.6", 
    cols, 
    rect_bounds = NULL,
    boundaries = NULL,
    crop = FALSE
) {
  # precompute cropped FOVs once for all Seurat objects
  roi_list <- NULL
  if (crop && !is.null(rect_bounds)) {
    roi_list <- lapply(seurat_list, function(obj) {
      crop_fov(
        obj[[fov]],
        x = c(rect_bounds$xmin, rect_bounds$xmax),
        y = c(rect_bounds$ymin, rect_bounds$ymax)
      )
    })
  }
  
  # create rectangle layer if provided
  rect_layer <- NULL
  if (is.null(rect_bounds)) {
    rect_layer <- geom_rect(
      aes(xmin = rect_bounds$xmin, xmax = rect_bounds$xmax,
          ymin = rect_bounds$ymin, ymax = rect_bounds$ymax),
      fill = NA,
      color = "black"
    )
  } 
  
  image_plots <- lapply(seq_along(seurat_list), function(i) {
    obj <- seurat_list[[i]]
    
    # update fov if roi is provided
    if (crop && !is.null(roi_list)) {
      obj[["roi"]] <- roi_list[[i]]
      fov_use <- "roi"
    } else {
      fov_use <- fov
    }
    
    # create the image dimensional plot
    p <- ImageDimPlot(
      obj, fov = fov_use, boundaries = boundaries, 
      group.by = group.by, cols = cols, axes = FALSE
    ) + labs(title = "") + rect_layer
    
    # keep legend only on the last plot
    if (i != length(seurat_list)) {
      p <- p + NoLegend()
    }
    
    return(p)
  })
  
  # combine plots in one row
  combined_images <- patchwork::wrap_plots(image_plots, nrow = 1)
  
  return(combined_images)
}

