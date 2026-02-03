# Required packages are Seurat, data.table, dplyr, ggplot2

# helper function to validate if the length of given seurat list matches the length of segmentation options
validate_segmentation <- function(segmentation, seg_levels, n_objs) {
  
  if (!all(segmentation %in% seg_levels)) {
    bad <- unique(segmentation[!segmentation %in% seg_levels])
    stop(
      "Invalid segmentation value(s): ",
      paste(bad, collapse = ", "),
      "\nAllowed values: ",
      paste(seg_levels, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!(length(segmentation) == n_objs)) {
    stop(
      "segmentation must have length 1 or match merged_seurat_list length (",
      n_objs, ").",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

generate_sample_df <- function(merged_seurat_list, segmentation, seg_levels, count_col) {
  # check if segmentation options match
  n <- length(merged_seurat_list)
  validate_segmentation(segmentation, seg_levels, n)
  
  segmentation <- factor(segmentation, levels = seg_levels)
  
  # use data.table for faster grouping and aggregation 
  data.table::rbindlist(
    # iterate over each seurat object and segmentation in parallel 
    # syntax Map(function(a, b) {...}, A, B)
    Map(function(obj, seg) {
      dt <- data.table::as.data.table(obj@meta.data)
      # group cells by sample_id and compute summaries per sample
      # syntax dt[i, j, by]
      dt[, .(
        cell_count = .N, # number of rows in current group, faster than n()
        transcript_count = sum(get(count_col)),
        batch = batch[1],
        segmentation = seg
      ), by = sample_id]
    }, merged_seurat_list, segmentation),
    use.names = TRUE
  )
}

# This function is only suitable for Xenium data
generate_transcripts_assigned <- function(dir_list, segmentation) {
  # check if segmentation options match
  seg_levels <- c("Nuclear expansion", "Proseg", "Cellpose2")
  n <- length(dir_list)
  validate_segmentation(segmentation, seg_levels, n)
  
  segmentation <- factor(segmentation, levels = seg_levels)
  
  data.table::rbindlist(
    Map(function(dir, seg) {
      summary_info <- read.csv(dir)
      if (grepl("Proseg", dir)) {
        summary_info <- summary_info %>% 
          filter(
            obj_name %in% c(
              "Batch24__0011456__Region_1.rds",
              "Batch24__0011456__Region_2.rds",
              "Batch24__0011456__Region_3.rds",
              "Batch24__0011456__Region_4.rds",
              "Batch27__0017329__Region_1.rds",
              "Batch27__0017329__Region_2.rds",
              "Batch27__0017329__Region_3.rds",
              "Batch27__0017329__Region_4.rds",
              "Batch27__0017329__Region_5.rds"
            )
          )
      }
      dt <- data.table::as.data.table(summary_info)
      dt[, .(
        pt_transcripts_assigned, 
        sample,
        segmentation = seg
      )]
    }, dir_list, segmentation),
    use.names = TRUE
  )
}

generate_cell_df <- function(merged_seurat_list, segmentation, seg_levels, count_col, feature_col) {
  # check if segmentation options match
  n <- length(merged_seurat_list)
  validate_segmentation(segmentation, seg_levels, n)
  
  segmentation <- factor(segmentation, levels = seg_levels)
  
  # extract cell information sample-by-sample
  data.table::rbindlist(
    Map(function(obj, seg) {
      # extract cell metadata
      dt <- data.table::as.data.table(obj@meta.data)
      dt[, .(
        sample = gsub("_", "", sample_id),
        nCounts = get(count_col),
        nFeatures = get(feature_col),
        cell_type = ScType_label_res.0.6,
        cell_area = if("cell_area" %in% colnames(dt)) cell_area else NA_real_,
        aspect_ratio = if("aspect_ratio" %in% colnames(dt)) aspect_ratio else NA_real_,
        log10_signal_density = if("log10_signal_density" %in% colnames(dt)) log10_signal_density else NA_real_,
        log10_signal_density_outlier_sc = if("log10_signal_density_outlier_sc" %in% colnames(dt)) log10_signal_density_outlier_sc else NA_real_,
        segmentation = seg
      )]
    }, merged_seurat_list, segmentation),
    use.names = TRUE
  )
}

# --- plots for sample-level quality metrics
sample_boxplot <- function(df, metric, label) {
  ggplot(df, 
         aes(x = segmentation, y = .data[[metric]], fill = segmentation)) +
    geom_boxplot(outliers = FALSE) + 
    geom_jitter(width = 0.1, height = 0.1) + 
    # geom_line(aes(group = sample_id), colour = "darkred", linewidth = 0.5, alpha = 0.7) +
    scale_fill_manual(values = scales::hue_pal()(3)) +
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
    aes(x = sample, y = .data[[metric]], fill = segmentation)
  ) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(aes(group = interaction(sample, segmentation)), 
                 width = 0.15, outlier.shape = NA, fill = "white", color = "black",
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = scales::hue_pal()(3)) +
    labs(x = "Samples",
         y = label,
         fill = "Segmentation") + 
    scale_y_continuous(n.breaks = 8) + 
    theme_bw() + 
    theme(legend.position = "none")
}

cell_median_boxplot <- function(cell_df, metric, label, jitter_width = NULL, jitter_height = NULL) {
  median_df <- cell_df %>% 
    group_by(segmentation, sample) %>% 
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
    aes(x = segmentation, y = .data[["median_value"]], fill = segmentation)
  ) +
    geom_boxplot(outliers = FALSE) +
    jitter_layer +
    scale_fill_manual(values = scales::hue_pal()(3)) +
    labs(x = "Segmentation", 
         y = paste("Median", label)) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}