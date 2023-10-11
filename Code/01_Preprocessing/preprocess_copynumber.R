library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

copynumber_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")

get_chromosome_arms <- function(df, mart = ensembl_mart, symbol_col = "Gene.Symbol") {
  chr_location <- biomaRt::getBM(attributes = c("hgnc_symbol", 'chromosome_name', 'band'),
                                 filters = 'hgnc_symbol',
                                 values = unique(df[[symbol_col]]),
                                 mart = ensembl_mart) %>%
    rename(Gene.Symbol = "hgnc_symbol",
           Gene.Chromosome = "chromosome_name",
           Gene.ChromosomeBand = "band") %>%
    mutate(Gene.ChromosomeArm = ifelse(str_detect(Gene.ChromosomeBand, "q"), "q",
                                       ifelse(str_detect(Gene.ChromosomeBand, "p"), "p",
                                              NA))) %>%
    unite("Gene.ChromosomeArm", Gene.Chromosome, Gene.ChromosomeArm,
          sep = "", remove = FALSE) %>%
    mutate_all(na_if, "") %>%
    drop_na() %>%
    distinct() %>%
    # Remove genes that may be located on multiple chromosomes
    add_count(Gene.Symbol) %>%
    filter(n == 1) %>%
    select(-n)

  df <- df %>%
    mutate(Gene.Symbol = .data[[symbol_col]]) %>% # Ensure correct column name
    left_join(y = chr_location, by = "Gene.Symbol",
              na_matches = "never", relationship = "many-to-one")

  return(df)
}

df_celllines <- read_parquet(here(output_data_dir, "celllines.parquet"))

# ToDo: NOTE! Chromosomes 13p, 14p, 15p, 21p, 22p, X, Y missing!
df_arm_level_cna <- read_csv_arrow(here(copynumber_data_dir, "Arm-level_CNAs.csv")) %>%
  rename(CellLine.DepMapModelId = 1) %>%
  pivot_longer(everything() & !CellLine.DepMapModelId,
               names_to = "Gene.ChromosomeArm", values_to = "ChromosomeArm.CNA") %>%
  mutate(ChromosomeArm.CNA = as.integer(ChromosomeArm.CNA))

df_aneuploidy <- read_csv_arrow(here(copynumber_data_dir, "Aneuploidy.csv")) %>%
  rename(CellLine.DepMapModelId = 1,
         CellLine.AneuploidyScore = "Aneuploidy score",
         CellLine.Ploidy = "Ploidy",
         CellLine.WGD = "Genome doublings")

copy_number <- read.csv(
  here(copynumber_data_dir, "Copy_Number_Public_23Q2.csv"),
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix() %>%
  as.table() %>%
  as.data.frame() %>%
  setNames(c("CellLine.DepMapModelId", "Gene.Symbol", "Gene.CopyNumber")) %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId",
             na_matches = "never", relationship = "many-to-one") %>%
  get_chromosome_arms() %>%
  # ToDo: Consider Left Join
  inner_join(y = df_arm_level_cna, by = c("CellLine.DepMapModelId", "Gene.ChromosomeArm"),
             na_matches = "never", relationship = "many-to-one") %>%
  inner_join(y = df_aneuploidy, by = "CellLine.DepMapModelId",
             na_matches = "never", relationship = "many-to-one") %>%
  # filter(CellLine.Ploidy < 3.5) %>%     # Removing near-tetraploid cell lines is detrimental for factor prediction
  write_parquet(here(output_data_dir, 'copy_number.parquet'),
                version = "2.6")

# ToDo: Use different copy number datasets for ProCan and DepMap