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
# Plot Copy Number, WGD, AS, CNA, ...
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet"))

## Chromosome Arm CNA
### TODO: Do single CNA events occur more often than doubling/halving of chr number per chromosome?
### Events per chromosome
copy_number %>%
  # filter(CellLine.AneuploidyScore > 0) %>%
  distinct(Gene.Chromosome, Gene.ChromosomeArm, Gene.ChromosomeBand, ChromosomeArm.CNA, CellLine.CustomId) %>%
  mutate(Chromosome = fct_reorder(Gene.Chromosome, as.integer(Gene.Chromosome)),
         ChromosomeArm = factor(substr(Gene.ChromosomeBand, 1, 1), levels = c("p", "q"))) %>%
  count(Chromosome, ChromosomeArm, ChromosomeArm.CNA) %>%
  group_by(Chromosome, ChromosomeArm) %>%
  mutate(CNA = factor(ChromosomeArm.CNA),
         Share = n / sum(n),
         EventShare = as.integer(ChromosomeArm.CNA) * Share) %>%
  filter(CNA != 0) %>%
  ggplot() +
  aes(x = ChromosomeArm, y = EventShare, fill = CNA) +
  geom_bar(stat = "identity", position = "identity") +
  scale_y_continuous(limits = c(-0.5, 0.5),
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels = c("50%", "25%", "", "25%", "50%")) +
  facet_grid(~Chromosome)

### Total Number of CNA Events
copy_number %>%
  distinct(Gene.Chromosome, Gene.ChromosomeArm, Gene.ChromosomeBand, ChromosomeArm.CNA, CellLine.CustomId) %>%
  count(ChromosomeArm.CNA)

### Ploidy + CNA
copy_number %>%
  filter(CellLine.AneuploidyScore > 0) %>%
  distinct(Gene.ChromosomeArm, ChromosomeArm.CNA, CellLine.CustomId, CellLine.Ploidy) %>%
  mutate(ChromosomeArm.Count = as.factor(ChromosomeArm.CNA + round(CellLine.Ploidy))) %>%
  ggplot() +
  aes(x = ChromosomeArm.Count) +
  geom_bar()

### Gain/Loss Correlation
corr_matrix <- copy_number %>%
  distinct(Gene.Chromosome, Gene.ChromosomeArm, ChromosomeArm.CNA, CellLine.CustomId) %>%
  mutate(Gene.ChromosomeArm = fct_reorder(Gene.ChromosomeArm, as.integer(Gene.Chromosome))) %>%
  arrange(Gene.ChromosomeArm) %>%
  select(Gene.ChromosomeArm, ChromosomeArm.CNA, CellLine.CustomId) %>%
  pivot_wider(names_from = "Gene.ChromosomeArm", values_from = "ChromosomeArm.CNA", id_cols = "CellLine.CustomId") %>%
  select(-CellLine.CustomId) %>%
  psych::corr.test(method = "kendall", adjust = "BH")

corrplot(corr_matrix$r, p.mat = corr_matrix$p,
         type = "upper", tl.col = "black", tl.srt = 45,
         pch.cex = 1, pch.col = "darkgrey")

# TODO: Gain/Loss co-occurrence
## Step 1: Replace CNA with "Gain", "Neutral", "Loss" event
## Step 2: Drop neutral events
## Step 3: Combine event with chr arm
## Step 4: Pivot wider (cols: chr arm event, rows: cell lines, values: TRUE/FALSE)
## Step 5: Build Matrix - For each pair of ChrArmEvent increment by one if pair exists in a cell line

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
