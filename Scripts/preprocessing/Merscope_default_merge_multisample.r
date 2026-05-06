suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(future)
})
options(future.globals.maxSize = 8000  * 1024^3)

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

main_dir <- "/vast/projects/SpatialBench/segmentation/Cellpose/resegmentation_cellpose_merscope/MERSCOPE_cellpose_reseg/seurat_objects/Vizgen"

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


# per‐sample QC
processed_list <- map(seurat_list, function(so){
  DefaultAssay(so) <- "Vizgen"
  
  # QC filter based on median 10th percentile of gene counts across 9 samples
  if (names(so@images) == "wt709") {
    so <- subset(so, subset = nCount_Vizgen > 3)
  } else {
    so <- subset(so, subset = nCount_Vizgen > 10)
  }
  
  # SCTransform normalisation
  so <- SCTransform(so, assay = "Vizgen", clip.range = c(-10, 10))
  
  so
})

# name your list by sample_id so merge() can use those prefixes
names(processed_list) <- sample_info$sample_id

# now merge using those names
merged_seurat <- merge(
  x           = processed_list[[1]],
  y           = processed_list[-1],
  add.cell.ids= names(processed_list),
  project     = "vizgen_ctrl_wt_ko"
)

# sanity‐check: cell barcodes should now start with your sample_id prefixes
head(colnames(merged_seurat))

# and your metadata columns are all preserved
head(
  merged_seurat@meta.data[, 
                          c("sample","sample_id","stim","batch","slide")
  ]
)


# choose your output directory
out_dir <- "/vast/scratch/users/zhang.ji/data/default_segmentation/Merscope_seurat"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# save as an RDS
saveRDS(merged_seurat,
        file = file.path(out_dir, "vizgen_ctrl_wt_ko.rds"))

