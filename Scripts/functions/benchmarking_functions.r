# Required packages are Seurat, Matrix, data.table, dplyr, ggplot2

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

generate_sample_df <- function(obj, segmentation, count_col) {
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
    segmentation = segmentation
  ), by = sample_id]
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

generate_cell_df <- function(obj, segmentation, count_col, feature_col) {
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
    segmentation = segmentation
  )]
}

extract_sample_and_cell_df <- function(
    file_path,
    segmentation,
    count_col,
    feature_col
) {
  obj <- readRDS(file_path)
  
  if (grepl("Xenium", file_path)) {
    obj <- add_common_gene_counts(obj, common_genes)
  }
  
  sample_df <- generate_sample_df(
    obj,
    segmentation = segmentation,
    count_col = count_col
  )
  
  cell_df <- generate_cell_df(
    obj,
    segmentation = segmentation,
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
        list(file_path, segmentation, count_col, feature_col),
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
  "Cellpose1" = "#F8766D",
  "Proseg" = "#00BA38",
  "Cellpose2" = "#619CFF",
  "Nuclear expansion" = "#00BFC4"
)

sample_boxplot <- function(df, metric, label) {
  ggplot(df, 
         aes(x = segmentation, y = .data[[metric]], fill = segmentation)) +
    geom_boxplot(outliers = FALSE) + 
    geom_jitter(width = 0.1, height = 0.1) + 
    # geom_line(aes(group = sample_id), colour = "darkred", linewidth = 0.5, alpha = 0.7) +
    scale_fill_manual(values = colour_panel) +
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
    aes(x = sample, y = .data[[metric]], fill = segmentation)
  ) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(aes(group = interaction(sample, segmentation)), 
                 width = 0.15, outlier.shape = NA, fill = "white", color = "black",
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = colour_panel) +
    facet_wrap(~platform_version) + 
    labs(x = "Samples",
         y = label,
         fill = "Segmentation") + 
    scale_y_continuous(n.breaks = 8) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

cell_median_boxplot <- function(cell_df, metric, label, jitter_width = NULL, jitter_height = NULL) {
  median_df <- cell_df %>% 
    group_by(platform_version, segmentation, sample) %>% 
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
