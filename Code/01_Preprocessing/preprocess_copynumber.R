library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "visualization.R"))

copynumber_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "CopyNumber")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

df_celllines <- read_parquet(here(output_data_dir, "celllines.parquet"))

# ToDo: Add sources

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

copy_number_absolute <- read.csv(
  here(copynumber_data_dir, "Copy_Number_(Absolute).csv"),
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix() %>%
  as.table() %>%
  as.data.frame() %>%
  setNames(c("CellLine.DepMapModelId", "Gene.Symbol", "Gene.CopyNumber"))

copy_number <- read.csv(
  here(copynumber_data_dir, "Copy_Number_Public_23Q2.csv"),
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix() %>%
  as.table() %>%
  as.data.frame() %>%
  setNames(c("CellLine.DepMapModelId", "Gene.Symbol", "Gene.CopyNumber.Relative")) %>%
  inner_join(y = copy_number_absolute, by = c("CellLine.DepMapModelId", "Gene.Symbol"),
             na_matches = "never", relationship = "one-to-one") %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId",
             na_matches = "never", relationship = "many-to-one") %>%
  get_chromosome_arms() %>%
  # ToDo: Consider Left Join
  inner_join(y = df_arm_level_cna, by = c("CellLine.DepMapModelId", "Gene.ChromosomeArm"),
             na_matches = "never", relationship = "many-to-one") %>%
  inner_join(y = df_aneuploidy, by = "CellLine.DepMapModelId",
             na_matches = "never", relationship = "many-to-one") %>%
  # filter(CellLine.Ploidy < 3.5) %>%     # Removing near-tetraploid cell lines is detrimental for factor prediction
  updateGeneSymbols() %>%
  write_parquet(here(output_data_dir, 'copy_number.parquet'),
                version = "2.6")

# Split dataset based on whole genome doubling status
copy_number_wgd <- copy_number %>%
  filter(CellLine.WGD > 0 & CellLine.Ploidy >= 2.9) %>%
  write_parquet(here(output_data_dir, 'copy_number_wgd.parquet'),
                version = "2.6")

copy_number_no_wgd <- copy_number %>%
  filter(CellLine.WGD == 0 | CellLine.Ploidy < 2.9) %>%
  write_parquet(here(output_data_dir, 'copy_number_no-wgd.parquet'),
                version = "2.6")

# ToDo: Use different copy number datasets for ProCan and DepMap

# === Evaluation ===
# TODO: Plot Copy Number, WGD, AS, CNA, ...

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet"))

## Copy Number distribution per WGD and Ploidy
copy_number %>%
  violin_plot(as.factor(CellLine.WGD), Gene.CopyNumber.Relative)
copy_number %>%
  violin_plot(as.factor(CellLine.WGD), Gene.CopyNumber)

## Aneuploidy Score Histogramm
df_aneuploidy %>%
  ggplot() +
  aes(x = CellLine.AneuploidyScore, fill = as.factor(CellLine.WGD)) +
  geom_histogram(position = "dodge")

# See description of copy number data: CN_rel = log2(CN_Ratio + 1)
## https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q4&filename=OmicsCNGene.csv
## https://github.com/broadinstitute/depmap_omics
# Assume that CN ratio is relative to rounded ploidy of cell line
copy_number_derived <- copy_number %>%
  mutate(CellLine.WGD = as.factor(CellLine.WGD),
         CellLine.Ploidy.Rounded = as.integer(round(CellLine.Ploidy)),
         Gene.CopyNumber.Ratio = (2^Gene.CopyNumber.Relative) - 1,
         Gene.CopyNumber.Derived = Gene.CopyNumber.Ratio * CellLine.Ploidy.Rounded)

copy_number_derived %>%
  group_by(CellLine.Ploidy.Rounded) %>%
  summarize(Gene.CopyNumber.Ratio.Median = median(Gene.CopyNumber.Ratio, na.rm = TRUE),
            Gene.CopyNumber.Derived.Median = median(Gene.CopyNumber.Derived, na.rm = TRUE))

# TODO: Plot Distances
copy_number_distance <- copy_number_derived %>%
  group_by(Gene.Symbol) %>%
  summarize(Gene.CopyNumber.Distance.Wasserstein = transport::wasserstein1d(Gene.CopyNumber.Derived, Gene.CopyNumber),
            Gene.CopyNumber.Distance.Mean = mean(Gene.CopyNumber.Derived - Gene.CopyNumber))
