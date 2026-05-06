suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyverse)
  library(Seurat)
  library(ggplot2)
})

source("~/Scripts/functions/processing.r")

# list file directories and their metadata
sample_info <- tribble(
  ~file,                            ~sample, ~stim,  ~batch,    ~slide,
  "Batch24__0011456__Region_1.rds", "166",   "ko",   "batch24", "slide0011456",
  "Batch24__0011456__Region_2.rds", "167",   "ko",   "batch24", "slide0011456",
  "Batch27__0017329__Region_3.rds", "168",   "ko",   "batch27", "slide0017329",
  "Batch27__0017329__Region_2.rds", "172",   "ctrl", "batch27", "slide0017329",
  "Batch24__0011456__Region_4.rds", "173",   "ctrl", "batch24", "slide0011456",
  "Batch27__0017329__Region_4.rds", "174",   "ctrl", "batch27", "slide0017329",
  "Batch27__0017329__Region_5.rds", "709",   "wt",   "batch27", "slide0017329",
  "Batch27__0017329__Region_1.rds", "710",   "wt",   "batch27", "slide0017329",
  "Batch24__0011456__Region_3.rds", "713",   "wt",   "batch24", "slide0011456"
) %>%
  mutate(
    sample_id = paste(stim, sample, sep = "_")
  )

main_dir <- "/vast/projects/SpatialBench/segmentation/Proseg/updated_seurat_objects/Xenium"

xenium_uni_proseg_list <- create_seurat_list(sample_info, main_dir)
names(xenium_uni_proseg_list) <- sample_info$sample_id

nCounts <- lapply(seq_along(xenium_uni_proseg_list), function(i) {
  obj <- xenium_uni_proseg_list[[i]]
  data.frame(
    nCount_Xenium = obj$nCount_Xenium,
    Sample = names(xenium_uni_proseg_list)[i]
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

summary_table %>% 
  summarise(
    median_nCount_10pt = median(nCount_Xenium_10pt)
  )

xenium_proseg_qc_plot <- ggplot(nCounts, aes(x = Sample, y = nCount_Xenium)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 3) +
  geom_hline(yintercept = 31, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 9, linetype = "dashed", colour = "blue") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/vast/scratch/users/zhang.ji/data/proseg_segmentation/Xenium_seurat/unimodal/xenium_proseg_qc_plot.pdf", plot = xenium_proseg_qc_plot, width = 4, height = 4)

sapply(xenium_uni_proseg_list, function(obj) mean(obj$nCount_Xenium > 32))
sapply(xenium_uni_proseg_list, function(obj) mean(obj$nCount_Xenium > 6))
sapply(xenium_uni_proseg_list, function(obj) mean(obj$nCount_Xenium > 11))
