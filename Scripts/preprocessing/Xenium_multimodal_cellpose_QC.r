suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyverse)
  library(Seurat)
  library(ggplot2)
})

source("~/Scripts/functions/processing.r")

# Create a table of sample information 
main_dir <- "/vast/projects/SpatialBench/analysis/cellpose_segmentation/xenium_batch34_cellpose/seu_objects"
obj_dir <- list.files(main_dir, pattern = "rds$")
sample_info <- tibble(file = obj_dir) %>% 
  mutate(
    sample_id = case_when(
      grepl("Region_1", file) ~ "ctrl_173",
      grepl("Region_2", file) ~ "wt_710",
      grepl("Region_3", file) ~ "wt_709",
      grepl("Region_4", file) ~ "ctrl_174",
      grepl("Region_5", file) ~ "ctrl_172",
      grepl("Region_6", file) ~ "embryo"
    ),
    batch = "batch34",
    slide = "slide0032118"
  ) %>% 
  filter(sample_id != "embryo") %>% 
  tidyr::separate(col = sample_id, into = c("stim", "sample"), sep = "_", remove = FALSE)

# Create a seurat list of all batch 34 samples
multi_cellpose_list <- create_seurat_list(sample_info, main_dir)
names(multi_cellpose_list) <- sample_info$sample_id

nCounts <- lapply(seq_along(multi_cellpose_list), function(i) {
  obj <- multi_cellpose_list[[i]]
  obj@meta.data[, c("nCount_Xenium", "nCount_BlankCodeword", "nCount_ControlCodeword", "nCount_ControlProbe")] %>% 
    mutate(Sample = names(multi_cellpose_list)[i])
})
nCounts <- bind_rows(nCounts) 

summary_table <- nCounts %>% 
  group_by(Sample) %>% 
  summarise(nCount_BlankCodeword_90pt = quantile(nCount_BlankCodeword, 0.9),
            nCount_ControlCodeword_90pt = quantile(nCount_ControlCodeword, 0.9),
            nCount_ControlProbe_90pt = quantile(nCount_ControlProbe, 0.9),
            nCount_Xenium_10pt = quantile(nCount_Xenium, 0.1),
            nCount_Xenium_25pt = quantile(nCount_Xenium, 0.25),
            nCount_Xenium_50pt = quantile(nCount_Xenium, 0.5))
write.csv(summary_table, file = "/vast/scratch/users/zhang.ji/data/cellpose_segmentation/Xenium_seurat/xenium_qc_metrics.csv", row.names = F)

summary_table %>% 
  summarise(median_nBlankCodeword_90pt = median(nCount_BlankCodeword_90pt),
            median_nControlCodeword_90pt = median(nCount_ControlCodeword_90pt),
            median_nControlProbe_90pt = median(nCount_ControlProbe_90pt),
            median_nCount_10pt = median(nCount_Xenium_10pt))

# ggplot(nCounts, aes(x = Sample, y = nCount_BlankCodeword)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   geom_hline(yintercept = 141, linetype = "dashed", colour = "red") + 
#   annotate("text", label = "median nCount_Xenium = 141", colour = "red", x = 3.5, y = 125, size = 3) + 
#   theme_bw()
# 
# ggplot(nCounts, aes(x = Sample, y = nCount_ControlCodeword)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   geom_hline(yintercept = 141, linetype = "dashed", colour = "red") + 
#   annotate("text", label = "median nCount_Xenium = 141", colour = "red", x = 3.5, y = 135, size = 3) + 
#   theme_bw()
# 
# ggplot(nCounts, aes(x = Sample, y = nCount_ControlProbe)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   geom_hline(yintercept = 141, linetype = "dashed", colour = "red") + 
#   annotate("text", label = "median nCount_Xenium = 141", colour = "red", x = 3.5, y = 135, size = 3) + 
#   theme_bw()
# 
# ggplot(nCounts, aes(x = Sample, y = nCount_Xenium)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   # geom_hline(yintercept = 9, linetype = "dashed", colour = "red") + 
#   theme_bw()

# plot them on the same figure
nCounts_long <- nCounts %>% 
  pivot_longer(
    cols = starts_with("nCount_"),
    names_to = "Count type",
    values_to = "Counts"
  )

xenium_qc_plot <- ggplot(nCounts_long, aes(x = Sample, y = Counts)) + 
  geom_violin(trim = FALSE, scale = "width", adjust = 3) +
  geom_hline(yintercept = 9, linetype = "dashed", colour = "red") + 
  facet_wrap(~`Count type`, scales = "free_y") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/vast/scratch/users/zhang.ji/data/cellpose_segmentation/Xenium_seurat/xenium_qc_plot.pdf", plot = xenium_qc_plot, width = 8, height = 4)

sapply(multi_cellpose_list, function(obj) mean(obj$nCount_Xenium > 10))
