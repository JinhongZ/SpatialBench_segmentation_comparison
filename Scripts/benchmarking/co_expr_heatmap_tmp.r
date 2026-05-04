library(Seurat)
library(Matrix)
library(scales)
library(dplyr)
library(ggplot2)
library(patchwork)

make_coexp_plot <- function(df, xlab, ylab, title, bins = 20) {
  ggplot(df, aes(x = x, y = y)) +
    stat_bin2d(bins = bins, aes(fill = after_stat(count))) +
    scale_fill_gradient2(
      low  = "blue", mid = "white", high = "red",
      midpoint = 1,   # adjust to taste
      limits = c(1, 10^7),
      oob = scales::squish,
      trans    = "log10",
      labels   = label_number(),             # shows actual counts on colorbar
      name = NULL
    ) +
    scale_x_continuous(limits = c(0, 185)) + 
    scale_y_continuous(limits = c(0, 270)) + 
    labs(title = title, x = xlab, y = ylab) +
    # theme_classic(base_size = 10) +
    theme(
      plot.title   = element_text(size = 9, hjust = 0.5),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width  = unit(0.35, "cm")
    )
}

# ── Pre-bin the data and compute percentages ──────────────────────────────────
make_coexp_pt_plot <- function(df, xlab, ylab, title, bins = 20,
                            xlim = NULL, ylim = NULL,
                            xbreaks = waiver(), ybreaks = waiver()) {
  
  # Manually bin and compute percentage
  parse_mid <- function(bin_factor) {
    sapply(as.character(bin_factor), \(b) {
      nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", b), ",")[[1]])
      mean(nums)
    })
  }
  
  df_binned <- df |>
    mutate(
      x_bin = cut(x, breaks = seq(0, max(x), length.out = 25 + 1), include.lowest = TRUE),
      y_bin = cut(y, breaks = seq(0, max(y), length.out = 25 + 1), include.lowest = TRUE)
    ) |>
    count(x_bin, y_bin) |>
    mutate(
      pct   = n / sum(n) * 100,
      x_mid = parse_mid(x_bin),
      y_mid = parse_mid(y_bin)
    )
  
  ggplot(df_binned, aes(x = x_mid, y = y_mid, fill = pct)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0.004,          # adjust: white appears at this % value
      limits   = c(0, 0.005),
      labels   = function(x) paste0(x, "%"),
      name     = NULL
    ) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = ylim, breaks = ybreaks) +
    labs(title = title, x = xlab, y = ylab)
}

setwd("/vast/scratch/users/zhang.ji/data/cellpose_segmentation")
merscope_cellpose <- readRDS("Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds")
xenium_multi_proseg <- readRDS("../proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds")
xenium_multi_cellpose <- readRDS("Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds")

mer_counts <- GetAssayData(merscope_cellpose, layer = "counts")
xen_proseg_counts <- GetAssayData(xenium_multi_proseg, layer = "counts")
xen_counts <- GetAssayData(xenium_multi_cellpose, layer = "counts")

gene_x <- "Cd19"
gene_y <- "S100a9"

mer_df <- data.frame(
  x = as.numeric(mer_counts[gene_x, ]),
  y = as.numeric(mer_counts[gene_y, ])
)
xen_proseg_df <- data.frame(
  x = as.numeric(xen_proseg_counts[gene_x, ]),
  y = as.numeric(xen_proseg_counts[gene_y, ])
)
xen_df <- data.frame(
  x = as.numeric(xen_counts[gene_x, ]),
  y = as.numeric(xen_counts[gene_y, ])
)

p1 <- make_coexp_plot(mer_df, gene_x, gene_y, "MERSCOPE", bins = 25)
p2 <- make_coexp_plot(xen_df, gene_x, gene_y, "Xenium multimodal Cellpose", bins = 25)
p3 <- make_coexp_plot(xen_proseg_df, gene_x, gene_y, "Xenium multimodal Proseg", bins = 25)
p4 <- make_coexp_pt_plot(xen_df, gene_x, gene_y, "Xenium multimodal Cellpose", bins = 25)
p5 <- make_coexp_pt_plot(xen_proseg_df, gene_x, gene_y, "Xenium multimodal Proseg", bins = 25)

p1 | p2
p3 | p2 # ensure two plots use same axis scale
p5 | p4

