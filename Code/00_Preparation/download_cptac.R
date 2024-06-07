library(here)
library(TCGAbiolinks)
library(httr)
library(tools)
library(dplyr)
library(arrow)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

output_data_dir <- output_data_base_dir
downloads_dir <- downloads_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(downloads_dir, recursive = TRUE)

api_url <- "https://api.gdc.cancer.gov/data/"

# For parameters see https://portal.gdc.cancer.gov/analysis_page?app=Downloads
## Cohort Builder -> Repository
query_CPTAC3 <- GDCquery(project = "CPTAC-3",
                         workflow.type = "AscatNGS",
                         data.category = "Copy Number Variation",
                         data.type = "Gene Level Copy Number")

# TODO: Remove submitter IDs not present in proteomics data before downloading
query_results <- query_CPTAC3$results[[1]]

cn_datasets <- list()

pb <- txtProgressBar(min = 0, max = nrow(query_results), style = 3)
for (i in seq_len(nrow(query_results))) {
  dataset <- query_results[i, ]
  download_result <- httr::GET(paste0(api_url, dataset$id))
  download_content <- httr::content(download_result, "raw")
  writeBin(download_content, here(downloads_dir, dataset$file_name))

  # Check MD5 sum
  download_md5 <- tools::md5sum(here(downloads_dir, dataset$file_name))[[1]]
  if (dataset$md5sum != download_md5)
    stop(paste("MD5 sum for dataset", dataset$id, "does not match!"))

  # Read downloaded file
  max_fields <- max(count.fields(here(downloads_dir, dataset$file_name)))
  cn_datasets[[dataset$id]] <- read.table(here(downloads_dir, dataset$file_name),
                                          header = FALSE, fill = TRUE,
                                          col.names = paste0("V", seq_len(max_fields))) %>%
    janitor::row_to_names(1) %>%
    mutate(Model.ID = dataset$cases.submitter_id)

  setTxtProgressBar(pb, pb$getVal() + 1)
}
close(pb)

df_cn_ascat <- bind_rows(cn_datasets) %>%
  write_parquet(here(output_data_dir, 'gdc_cptac-3_cnv_ascat.parquet'), version = "2.6")
