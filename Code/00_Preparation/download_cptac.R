
library(here)
library(TCGAbiolinks)

here::i_am("DosageCompensationFactors.Rproj")

download_dir <- here("Downloads")
dir.create(download_dir, recursive = TRUE)

# For parameters see https://portal.gdc.cancer.gov/analysis_page?app=Downloads
## Cohort Builder -> Repository
query_CPTAC3 <- GDCquery(project = "CPTAC-3",
                                       workflow.type = "AscatNGS",
                                       data.category = "Copy Number Variation",
                                       data.type = "Gene Level Copy Number")


setwd(download_dir)

GDCdownload(query, method = "client", directory = download_dir)

cn_data <- GDCprepare(query)

setwd(here())