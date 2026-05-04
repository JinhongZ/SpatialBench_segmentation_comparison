library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)


MECR_sc_df_filtered <- MECR_sc_df_filtered %>% 
  group_by(cell_pair) %>% 
  mutate(n_pairs = n(), method_combo = "scRNA-seq") %>% 
  ungroup()
MECR_mat <- MECR_df_filtered %>%
  unite(method_combo, platform_version, segmentation, sep = " ") %>%
  bind_rows(MECR_sc_df_filtered) %>% 
  filter(n_pairs > 10) %>% 
  group_by(method_combo, cell_pair) %>% 
  summarise(MECR_median = median(MECR, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(
    names_from = method_combo,
    values_from = MECR_median
  ) %>%
  column_to_rownames("cell_pair") %>% 
  t()

head(MECR_mat)

row_order <- c("scRNA-seq", "MERSCOPE Default", "MERSCOPE Proseg", "MERSCOPE Cellpose2",
               "Xenium multimodal Default", "Xenium multimodal Proseg", "Xenium multimodal Cellpose2")

# cluster manually
mat4clust <- MECR_mat[row_order, ]
d <- dist(t(mat4clust), method = "euclidean")
hc <- hclust(d, method = "ward.D2")

col_order <- colnames(mat4clust)[hc$order]

n_pairs_mat <- MECR_df_filtered %>% 
  group_by(cell_pair) %>% 
  summarise(n_pairs = first(n_pairs)) %>% 
  # Build a 1-row matrix matching column order of mecr_matrix
  { setNames(.$n_pairs, .$cell_pair) } %>%
  .[col_order] %>%
  matrix(nrow = 1, dimnames = list("N pairs", names(.))) 

n_max <- max(n_pairs_mat, na.rm = TRUE)
col_fun_n <- colorRamp2(
  c(0, n_max / 2, n_max),
  c("white", "#6BAED6", "#08306B")
)

# --- N pairs annotation heatmap (bottom) ---
npairs_heatmap <- Heatmap(
  n_pairs_mat,
  name = "N pairs",
  col  = col_fun_n,
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- round(n_pairs_mat[i, j])
    text_colour <- ifelse(value > 20, "white", "black")
    grid.text(
      label = value,
      x, y,
      gp = gpar(fontsize = 8, col = text_colour)
    )
  },
  show_row_names    = TRUE,
  show_column_names = FALSE,   # shared with main heatmap above
  row_names_side    = "right",
  column_order      = col_order,  # lock column order
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  height            = unit(1, "cm"),
  heatmap_legend_param = list(
    title          = "N pairs",
    at             = c(0, 10, 20, 30, 40),
    direction      = "horizontal",
    title_position = "lefttop"
  )
)

# --- Main MECR heatmap ---
mecr_heatmap <- Heatmap(
  MECR_mat[, col_order],
  name = "Median MECR",
  col  = colorRamp2(
    seq(0, 1, length.out = 5),
    c("white", "#FFE5CC", "#FFB347", "#FF8C00", "#FF6600")
  ),
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  row_names_side    = "right",
  column_names_side = "top",
  column_names_rot  = 45,
  column_names_gp   = gpar(fontsize = 10),
  row_names_gp      = gpar(fontsize = 10),
  rect_gp           = gpar(col = "white", lwd = 1.5),
  heatmap_legend_param = list(
    title          = "Median MECR",
    at             = c(0, 0.25, 0.5, 0.75, 1),
    direction      = "horizontal",
    title_position = "lefttop"
  )
)

# --- Draw stacked ---
draw(
  mecr_heatmap %v% npairs_heatmap,   # %v% stacks vertically
  heatmap_legend_side   = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE
)


