suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyverse)
  library(Seurat)
  library(ggplot2)
})

source("~/Scripts/functions/processing.r")

# create a table of sample directories and metadata
main_dir <- "/vast/projects/SpatialBench/segmentation/Proseg/updated_seurat_objects/Xenium"
batch34_dir <- list.files(main_dir, pattern = "^Batch34__0032118")
sample_info <- tibble(file = batch34_dir) %>% 
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

# create a list of seurat objects
xenium_multi_proseg_list <- create_seurat_list(sample_info, main_dir)
names(xenium_multi_proseg_list) <- sample_info$sample_id

# extract the gene counts
nCounts <- lapply(seq_along(xenium_multi_proseg_list), function(i) {
  obj <- xenium_multi_proseg_list[[i]]
  data.frame(
    nCount_Xenium = obj$nCount_Xenium,
    Sample = names(xenium_multi_proseg_list)[i]
  )
})
nCounts <- bind_rows(nCounts) 

summary_table <- nCounts %>% 
  group_by(Sample) %>% 
  summarise(
    nCount_Xenium_10pt = quantile(nCount_Xenium, 0.1),
    nCount_Xenium_25pt = quantile(nCount_Xenium, 0.25),
    nCount_Xenium_50pt = quantile(nCount_Xenium, 0.5)
  )
write.csv(
  summary_table, 
  file = "/vast/scratch/users/zhang.ji/data/proseg_segmentation/Xenium_seurat/multimodal/xenium_proseg_qc_metrics.csv", 
  row.names = F
)

summary_table %>% 
  summarise(
    median_nCount_10pt = median(nCount_Xenium_10pt)
  )

xenium_proseg_qc_plot <- ggplot(nCounts, aes(x = Sample, y = nCount_Xenium)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 3) +
  geom_hline(yintercept = 31, linetype = "dashed", colour = "red") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/vast/scratch/users/zhang.ji/data/proseg_segmentation/Xenium_seurat/multimodal/xenium_proseg_qc_plot.pdf", 
       plot = xenium_proseg_qc_plot, width = 4, height = 4)

sapply(xenium_multi_proseg_list, function(obj) mean(obj$nCount_Xenium > 10))
