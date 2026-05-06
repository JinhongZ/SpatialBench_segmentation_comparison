suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dplyr)
})

setwd("~/data")

main.dir <- "/vast/projects/SpatialBench/data/xenium/G000218_Benchmarking"

obj.dir <- data.frame(
  obj = c("ctrl172", "ctrl173", "ctrl174", "ko166", "ko167", "ko168", "wt709", "wt710", "wt713"),
  batch = c("batch27", "batch34", "batch34", "batch24", "batch24", "batch27", "batch34", "batch27", "batch24"),
  sub.dir = c(
    "resegmented_xeniumranger_v2_0_1/resegmented_v2_0_1__0017329__Region_2/outs",
    "20241031__014144__20241031_G000218_Batch34_Xen/output-XETG00068__0032118__Region_1__20241031__014253",
    "20241031__014144__20241031_G000218_Batch34_Xen/output-XETG00068__0032118__Region_4__20241031__014253",
    "resegmented_xeniumranger_v2_0_1/resegmented_v2_0_1__0011456__Region_1/outs",
    "resegmented_xeniumranger_v2_0_1/resegmented_v2_0_1__0011456__Region_2/outs",
    "resegmented_xeniumranger_v2_0_1/resegmented_v2_0_1__0017329__Region_3/outs",
    "20241031__014144__20241031_G000218_Batch34_Xen/output-XETG00068__0032118__Region_3__20241031__014253",
    "resegmented_xeniumranger_v2_0_1/resegmented_v2_0_1__0017329__Region_1/outs",
    "resegmented_xeniumranger_v2_0_1/resegmented_v2_0_1__0011453__Region_3/outs"
  )
)
obj.dir$dir <- file.path(main.dir, obj.dir$sub.dir)

assigned_sum <- data.frame()
for (i in 1:nrow(obj.dir)) {
  name <- obj.dir$obj[i]
  dir <- obj.dir$dir[i]
  batch <- obj.dir$batch[i]
  
  # read transcript information
  if (grepl("resegmented", dir)) {
    transcript_df <- fread(file.path(dir, "transcripts.csv.gz"))
  } else {
    transcript_df <- read_parquet(file.path(dir, "transcripts.parquet"))
  }
  
  # create codeword_category if not exist
  if (!any(colnames(transcript_df) == "codeword_category")) {
    transcript_df <- transcript_df %>% 
      mutate(codeword_category = case_when(
        grepl("NegControlCodeword", feature_name) ~ "negative_control_codeword",
        grepl("NegControlProbe", feature_name) ~ "negative_control_probe",
        grepl("UnassignedCodeword", feature_name) ~ "unassigned_codeword",
        .default = "custom_gene"
      ))
  }
  
  pt_assigned <- transcript_df %>% 
    group_by(codeword_category) %>% 
    summarise(total_count = n(),
              assigned_count = sum(cell_id != "UNASSIGNED"),
              pt_assigned = assigned_count / total_count * 100,
              sample = name,
              batch = batch,
              segmentation = "Nuclear expansion")
  
  assigned_sum <- bind_rows(assigned_sum, pt_assigned)
}

write.csv(assigned_sum, file = "Xenium_assigned_summary.csv")
