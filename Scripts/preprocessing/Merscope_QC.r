suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyverse)
  library(Seurat)
  library(ggplot2)
})

# define files and their metadata
sample_info <- tribble(
  ~file,        ~sample, ~stim,  ~batch,    ~slide,
  "ko166.rds",  "166",   "ko",   "batch22", "slide1",
  "ko167.rds",  "167",   "ko",   "batch27", "slide1_2",
  "ko168.rds", "168",   "ko",   "batch26", "slide3_4",
  "ctrl172.rds","172",  "ctrl", "batch27", "slide1_2",
  "ctrl173.rds","173",   "ctrl", "batch26", "slide3_4",
  "ctrl174.rds","174",  "ctrl", "batch22", "slide1",
  "wt709.rds", "709",   "wt",   "batch33", "slide1",
  "wt710.rds", "710",   "wt",   "batch33", "slide1",
  "wt713.rds",  "713",   "wt",   "batch33", "slide1"
) %>%
  mutate(
    sample_id = paste(stim, sample, sep = "_")
  )

# main_dir <- "/vast/projects/SpatialBench/segmentation/Cellpose/resegmentation_cellpose_merscope/MERSCOPE_cellpose_reseg/seurat_objects/Vizgen"
main_dir <- "/vast/projects/SpatialBench/segmentation/Cellpose/resegmentation_cellpose_merscope/MERSCOPE_cellpose_reseg/seurat_objects/Cellpose"

# read each file and add metadata via AddMetaData()
seurat_list <- sample_info %>%
  # each row stores a Seurat object
  mutate(
    so = map(file, ~ readRDS(file.path(main_dir, .x))) # map each row of file into .x
  ) %>%
  mutate(
    # perform parallel map to iterate row-wise over columns in list 
    so = pmap(
      list(so, sample, sample_id, stim, batch, slide),
      function(so, sample, sample_id, stim, batch, slide) {
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
names(seurat_list) <- sample_info$sample_id

nCounts <- lapply(seq_along(seurat_list), function(i) {
  obj <- seurat_list[[i]]
  obj@meta.data[, c("nCount_Vizgen", "nCount_Blanks")] %>% 
    mutate(Sample = names(seurat_list)[i])
})
nCounts <- bind_rows(nCounts) 

summary_table <- nCounts %>% 
  group_by(Sample) %>% 
  summarise(nCount_Blanks_90pt = quantile(nCount_Blanks, 0.9),
            nCount_Vizgen_10pt = quantile(nCount_Vizgen, 0.1),
            nCount_Vizgen_25pt = quantile(nCount_Vizgen, 0.25),
            nCount_Vizgen_50pt = quantile(nCount_Vizgen, 0.5))

summary_table %>% 
  filter(nCount_Blanks_90pt != 0) %>% 
  summarise(median_nBlank_90pt = median(nCount_Blanks_90pt))

summary_table %>% 
  summarise(median_nCount_10pt = median(nCount_Vizgen_10pt))

ggplot(nCounts, aes(x = Sample, y = nCount_Blanks)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 3) +
  geom_hline(yintercept = 9, linetype = "dashed", colour = "red") + 
  annotate("text", label = "median nCount_Vizgen = 9", colour = "red", x = 3.5, y = 15, size = 3) + 
  theme_bw()

ggplot(nCounts, aes(x = Sample, y = nCount_Vizgen)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 3) +
  # geom_hline(yintercept = 9, linetype = "dashed", colour = "red") + 
  theme_bw()
