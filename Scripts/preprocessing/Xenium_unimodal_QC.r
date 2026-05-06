suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyverse)
  library(Seurat)
  library(ggplot2)
})

source("~/Scripts/functions/processing.r")

# define files and their metadata
sample_info <- tribble(
  ~file,        ~sample, ~stim,  ~batch,    ~slide,
  "ko166.rds",  "166",   "ko",   "batch24", "slide0011456",
  "ko167.rds",  "167",   "ko",   "batch24", "slide0011456",
  "ko168.rds", "168",   "ko",   "batch27", "slide0017329",
  "ctrl172.rds","172",  "ctrl", "batch27", "slide0017329",
  "ctrl173.rds","173",   "ctrl", "batch24", "slide0011456",
  "ctrl174.rds","174",  "ctrl", "batch27", "slide0017329",
  "wt709.rds", "709",   "wt",   "batch27", "slide0017329",
  "wt710.rds", "710",   "wt",   "batch27", "slide0017329",
  "wt713.rds",  "713",   "wt",   "batch24", "slide0011456"
) %>%
  mutate(
    sample_id = paste(stim, sample, sep = "_")
  )

# create a list of those samples
main_dir <- "~/data/default_segmentation/Xenium_seurat/unimodal"
xenium_uni_default_list <- create_seurat_list(sample_info, main_dir)
names(xenium_uni_default_list) <- sample_info$sample_id

nCounts <- lapply(seq_along(xenium_uni_default_list), function(i) {
  obj <- xenium_uni_default_list[[i]]
  obj@meta.data[, c("nCount_Xenium", "nCount_BlankCodeword", "nCount_ControlCodeword", "nCount_ControlProbe")] %>% 
    mutate(Sample = names(xenium_uni_default_list)[i])
})
nCounts <- bind_rows(nCounts) 

summary_table <- nCounts %>% 
  group_by(Sample) %>% 
  summarise(nCount_BlankCodeword_90pt = quantile(nCount_BlankCodeword, 0.9),
            nCount_BlankCodeword_max = max(nCount_BlankCodeword),
            nCount_ControlCodeword_90pt = quantile(nCount_ControlCodeword, 0.9),
            nCount_ControlCodeword_max = max(nCount_ControlCodeword),
            nCount_ControlProbe_90pt = quantile(nCount_ControlProbe, 0.9),
            nCount_ControlProbe_max = max(nCount_ControlProbe),
            nCount_Xenium_10pt = quantile(nCount_Xenium, 0.1),
            nCount_Xenium_25pt = quantile(nCount_Xenium, 0.25),
            nCount_Xenium_50pt = quantile(nCount_Xenium, 0.5))
write.csv(summary_table, file = file.path(main_dir, "xenium_qc_metrics.csv"), row.names = F)

summary_table %>% 
  summarise(median_nBlankCodeword_90pt = median(nCount_BlankCodeword_90pt),
            median_nControlCodeword_90pt = median(nCount_ControlCodeword_90pt),
            median_nControlProbe_90pt = median(nCount_ControlProbe_90pt),
            median_nCount_10pt = median(nCount_Xenium_10pt))

# ggplot(nCounts, aes(x = Sample, y = nCount_BlankCodeword)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   geom_hline(yintercept = 31, linetype = "dashed", colour = "red") + 
#   annotate("text", label = "median nCount_Xenium = 31", colour = "red", x = 3.5, y = 30, size = 3) + 
#   theme_bw()
# 
# ggplot(nCounts, aes(x = Sample, y = nCount_ControlCodeword)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   geom_hline(yintercept = 31, linetype = "dashed", colour = "red") + 
#   annotate("text", label = "median nCount_Xenium = 9", colour = "red", x = 3.5, y = 30, size = 3) + 
#   theme_bw()
# 
# ggplot(nCounts, aes(x = Sample, y = nCount_ControlProbe)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   geom_hline(yintercept = 31, linetype = "dashed", colour = "red") + 
#   annotate("text", label = "median nCount_Xenium = 9", colour = "red", x = 3.5, y = 30, size = 3) + 
#   theme_bw()
# 
# ggplot(nCounts, aes(x = Sample, y = nCount_Xenium)) +
#   geom_violin(trim = FALSE, scale = "width", adjust = 3) +
#   # geom_hline(yintercept = 31, linetype = "dashed", colour = "red") + 
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
  geom_hline(yintercept = 31, linetype = "dashed", colour = "red") + 
  geom_hline(yintercept = 9, linetype = "dashed", colour = "blue") + 
  facet_wrap(~`Count type`, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(main_dir, "xenium_qc_plot.pdf"), plot = xenium_qc_plot, width = 8, height = 4)

sapply(xenium_uni_default_list, function(obj) mean(obj$nCount_Xenium > 32))
sapply(xenium_uni_default_list, function(obj) mean(obj$nCount_Xenium > 5))
sapply(xenium_uni_default_list, function(obj) mean(obj$nCount_Xenium > 10))
