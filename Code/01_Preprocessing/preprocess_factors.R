library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(readxl)
library(assertr)
library(purrr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "visualization.R"))

factor_data_dir <- here(external_data_dir, "Factors")
phosphositeplus_data_dir <- here(factor_data_dir, "PhosphoSitePlus")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)

# === Load Datasets ===
## Load PhosphoSitePlus datasets
phospho_sites <- read_excel(here(phosphositeplus_data_dir, "Phosphorylation_site_dataset.xlsx"), skip = 3)
ubi_sites <- read_excel(here(phosphositeplus_data_dir, "Ubiquitination_site_dataset.xlsx"), skip = 3)
sumo_sites <- read_excel(here(phosphositeplus_data_dir, "Sumoylation_site_dataset.xlsx"), skip = 3)
methyl_sites <- read_excel(here(phosphositeplus_data_dir, "Methylation_site_dataset.xlsx"), skip = 3)
ace_sites <- read_excel(here(phosphositeplus_data_dir, "Acetylation_site_dataset.xlsx"), skip = 3)
reg_sites <- read_excel(here(phosphositeplus_data_dir, "Regulatory_sites.xlsx"), skip = 3)

## Load other factor datasets
# Hausser, Jean (2019), “Central dogma rates and the trade-off between precision and economy”, Mendeley Data, V1, doi: 10.17632/2vbrg3w4p3.1
human_rates <- read_excel(here(factor_data_dir, "humanRates.xlsx"), sheet = "H. sapiens rates")
# URL: http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php, version: v2.3
hippie <- read_excel(here(factor_data_dir, "HIPPIE-v2.3.xlsx"))
# URL: https://genome.ucsc.edu/cgi-bin/hgTables, genome: Human, assembly: Dec. 2013 (GRCh38/hg38), track: NCBI RefSeq
ref_seq <- read_excel(here(factor_data_dir, "NCBI.Human.RefSeq.chr1.xlsx"))
# DOI: https://doi.org/10.1038/s41467-018-03106-1, Supplementary Table 2
half_life <- read_excel(here(factor_data_dir, "41467_2018_3106_MOESM5_ESM.xlsx"))
# CORUM, version 210512, DOI: 10.1093/nar/gkm936
corum <- read_tsv_arrow(here(factor_data_dir, "CORUMallComplexes.tsv"))
# URL: https://mobidb.bio.unipd.it/browse?proteome=UP000005640&limit=10&skip=0&projection=acc,name,organism,reviewed,prediction-disorder-mobidblite.contentfraction,gene,length
# Version: 5.0, Release: 2022_07, Accessed: 2023-06-23
mobidb <- read_tsv_arrow(here(factor_data_dir, "mobidb_result.tsv"))
# DOI: 10.1016/j.cell.2016.09.015, Supplementary Table 4: All AHA Pulse-Chase and Related Data for the Human RPE-1 Cells
ned_human <- read_excel(here(factor_data_dir, "NED-Human_RPE-1.xlsx"))
# DOI: 10.1101/gr.1272403, Supplementary Table 9: Decay Rates (hour^-1) for Accessions in HepG2 Experiments
mrna_decay <- read_excel(here(factor_data_dir, "Yang.2003.mRNADecay.Rates.xlsx"))
# DOI: 10.1016/j.celrep.2013.09.043, Supplementary Table 2: Supersaturation Database.
agg_score <- read_excel(here(factor_data_dir, "AggregationScore.xlsx"), skip = 18254)

# === Prepare Factor Datasets ===
## Prepare PhosphoSitePlus data by counting occurrance of PTM and regulatory sites for each gene
count_sites <- function(df, colname = "n") {
  df %>%
    filter(ORGANISM == "human") %>%
    count(ACC_ID, GENE, name = colname) %>%
    rename(Protein.Uniprot.Accession = ACC_ID,
           Gene.Symbol = GENE) %>%
    mutate(across(where(is.character), toupper))
}

## Prepare human_rates data
df_rates <- human_rates %>%
  rename(Gene.ENSEMBL.Id = "...1",
         `mRNA Abundance` = "m",
         `Protein Abundance` = "p",
         `Transcription Rate` = "bm",
         `Translation Rate` = "bp",
         `Protein Length` = "lp",
         `mRNA Length` = "lm") %>%
  mapIds("ENSEMBL", "UNIPROT",
         "Gene.ENSEMBL.Id", "Protein.Uniprot.Accession") %>%
  mapIds("ENSEMBL", "SYMBOL",
         "Gene.ENSEMBL.Id", "Gene.Symbol")

## Prepare Proten-Protein-Interaction data
hippie_filtered <- hippie %>%
  select("Gene Name Interactor A", "Gene Name Interactor B", "Confidence Value") %>%
  rename(Gene.Symbol = "Gene Name Interactor A") %>%
  filter(`Confidence Value` > 0.6) %>%
  count(Gene.Symbol, name = "Protein-Protein Interactions") %>%
  drop_na() %>%
  mapIds("SYMBOL", "ENSEMBL",
         "Gene.Symbol", "Gene.ENSEMBL.Id") %>%
  mapIds("SYMBOL", "UNIPROT",
         "Gene.Symbol", "Protein.Uniprot.Accession")

## Prepare Protein half-life data
half_life_sample_cols <- c("Bcells replicate 1 half_life", "Bcells replicate 2 half_life",
                           "NK cells replicate 1 half_life", "NK cells replicate 2 half_life",
                           "Monocytes replicate 1 half_life", "Monocytes replicate 2 half_life",
                           "Mouse Neurons, replicate 3 half_life", "Mouse Neurons, replicate 4 half_life")

half_life_avg <- half_life %>%
  select("gene_name", all_of(half_life_sample_cols)) %>%
  rename(Gene.Symbol = "gene_name") %>%
  pivot_longer(all_of(half_life_sample_cols), names_to = "sample", values_to = "half_life_values") %>%
  group_by(Gene.Symbol) %>%
  summarize(`Protein Half-Life` = mean(half_life_values, na.rm = TRUE)) %>%
  ungroup() %>%
  mapIds("SYMBOL", "ENSEMBL",
         "Gene.Symbol", "Gene.ENSEMBL.Id") %>%
  mapIds("SYMBOL", "UNIPROT",
         "Gene.Symbol", "Protein.Uniprot.Accession")

## Prepare 3'-/5'-UTR data
df_utr <- ref_seq %>%
  separate(name, c("curation", "name"), sep = "_") %>%
  filter(curation == "NM") %>%
  mutate(`5'-UTR Length` = ifelse(strand == "+", abs(cdsStart - txStart),
                                  ifelse(strand == "-", abs(txEnd - cdsEnd), NA))) %>%
  mutate(`3'-UTR Length` = ifelse(strand == "+", abs(txEnd - cdsEnd),
                                  ifelse(strand == "-", abs(cdsStart - txStart), NA))) %>%
  select("name2", "3'-UTR Length", "5'-UTR Length") %>%
  rename(Gene.Symbol = "name2") %>%
  group_by(Gene.Symbol) %>%
  mutate(`Median 3'-UTR Length` = median(`3'-UTR Length`, na.rm = TRUE),
         `Median 5'-UTR Length` = median(`5'-UTR Length`, na.rm = TRUE),
         `Mean 3'-UTR Length` = mean(`3'-UTR Length`, na.rm = TRUE),
         `Mean 5'-UTR Length` = mean(`5'-UTR Length`, na.rm = TRUE),
         `3'-UTR Length` = list(`3'-UTR Length`),
         `5'-UTR Length` = list(`5'-UTR Length`)) %>%
  ungroup() %>%
  distinct() %>%
  mapIds("SYMBOL", "ENSEMBL",
         "Gene.Symbol", "Gene.ENSEMBL.Id") %>%
  mapIds("SYMBOL", "UNIPROT",
         "Gene.Symbol", "Protein.Uniprot.Accession")

## Prepare CORUM data
df_complexes <- corum %>%
  filter(Organism == "Human") %>%
  select(ComplexID, ComplexName, "subunits(UniProt IDs)") %>%
  separate_longer_delim("subunits(UniProt IDs)", delim = ";") %>%
  count(`subunits(UniProt IDs)`, name = "Protein Complexes (CORUM)") %>%
  rename(Protein.Uniprot.Accession = "subunits(UniProt IDs)") %>%
  mapIds("UNIPROT", "ENSEMBL",
         "Protein.Uniprot.Accession", "Gene.ENSEMBL.Id") %>%
  mapIds("UNIPROT", "SYMBOL",
         "Protein.Uniprot.Accession", "Gene.Symbol")

## Prepare MobiDB data
mobidb_features <- c("prediction-disorder-mobidb_lite", "prediction-low_complexity-merge",
                     "homology-domain-merge", "prediction-lip-anchor",
                     "prediction-polyampholyte-mobidb_lite_sub", "prediction-polar-mobidb_lite_sub")

df_mobidb <- mobidb %>%
  select(acc, feature, content_count) %>%
  filter(feature %in% mobidb_features) %>%
  pivot_wider(id_cols = "acc", names_from = "feature", values_from = "content_count") %>%
  rename(Protein.Uniprot.Accession = "acc",
         `Intrinsic Protein Disorder` = "prediction-disorder-mobidb_lite",
         `Low Complexity Score` = "prediction-low_complexity-merge",
         `Homology Score` = "homology-domain-merge",
         `Loops In Protein Score` = "prediction-lip-anchor",
         `Protein Polyampholyte Score` = "prediction-polyampholyte-mobidb_lite_sub",
         `Protein Polarity` = "prediction-polar-mobidb_lite_sub") %>%
  mapIds("UNIPROT", "ENSEMBL",
         "Protein.Uniprot.Accession", "Gene.ENSEMBL.Id") %>%
  mapIds("UNIPROT", "SYMBOL",
         "Protein.Uniprot.Accession", "Gene.Symbol")

## Prepare NED data
delta_score_col <- paste0(enc2utf8("\u0394"), "-score")
df_ned <- ned_human %>%
  select("Protein IDs (Uniprot)", "Gene names", delta_score_col) %>%
  rename(Protein.Uniprot.Accession = "Protein IDs (Uniprot)",
         Gene.Symbol = "Gene names",
         `Non-Exponential Decay Delta` = delta_score_col) %>%
  separate_longer_delim(Gene.Symbol, delim = ";") %>%
  drop_na() %>%
  mapIds("SYMBOL", "ENSEMBL",
         "Gene.Symbol", "Gene.ENSEMBL.Id")

## Prepare mRNA Decay data
df_decay <- mrna_decay %>%
  mutate_at(c("Rate_1", "Rate_2", "StdDev"), as.numeric) %>%
  mutate(`mRNA Decay Rate` = ifelse(is.na(Rate_1), Rate_2, Rate_1)) %>%
  drop_na(`mRNA Decay Rate`) %>%
  select("Accession", "mRNA Decay Rate") %>%
  rename(Gene.GenBank.Accession = "Accession") %>%
  mapIds("ACCNUM", "ENSEMBL",
         "Gene.GenBank.Accession", "Gene.ENSEMBL.Id") %>%
  mapIds("ACCNUM", "SYMBOL",
         "Gene.GenBank.Accession", "Gene.Symbol") %>%
  mapIds("ACCNUM", "UNIPROT",
         "Gene.GenBank.Accession", "Protein.Uniprot.Accession") %>%
  drop_na() %>%
  group_by(Gene.Symbol, Gene.ENSEMBL.Id, Protein.Uniprot.Accession) %>%
  summarize(`Mean mRNA Decay Rate` = mean(`mRNA Decay Rate`, na.rm = TRUE))

## Prepare Aggregation Score data
df_agg <- agg_score %>%
  select("Uniprot ID", "ZaggSC") %>%
  mutate(`Aggregation Score` = as.numeric(if_else(ZaggSC == "-", NA, ZaggSC))) %>%
  drop_na() %>%
  separate_wider_delim(`Uniprot ID`,
           delim = "_",
           names = c("Gene.Symbol", "Species")) %>%
  filter(Species == "human") %>%
  select(-ZaggSC, -Species) %>%
  mutate(Gene.Symbol = toupper(Gene.Symbol)) %>%
  mapIds("SYMBOL", "UNIPROT",
         "Gene.Symbol", "Protein.Uniprot.Accession") %>%
  mapIds("SYMBOL", "ENSEMBL",
         "Gene.Symbol", "Gene.ENSEMBL.Id")

# === Combine Factor Datasets ===

ptm_factor_datasets <- list(
  count_sites(phospho_sites, colname = "Phosphorylation Sites"),
  count_sites(ubi_sites, colname = "Ubiquitination Sites"),
  count_sites(sumo_sites, colname = "Sumoylation Sites"),
  count_sites(methyl_sites, colname = "Methylation Sites"),
  count_sites(ace_sites, colname = "Acetylation Sites"),
  count_sites(reg_sites, colname = "Regulatory Sites")
)

other_factor_datasets <- list(df_rates, hippie_filtered, half_life_avg, df_utr,
                           df_complexes, df_mobidb, df_ned, df_decay, df_agg)

df_dc_factors_ptm <- ptm_factor_datasets %>%
  reduce(full_join, by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
         relationship = "many-to-one") %>%
  filter(!(is.na(Protein.Uniprot.Accession) & is.na(Gene.Symbol))) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mapIds("UNIPROT", "ENSEMBL",
         "Protein.Uniprot.Accession", "Gene.ENSEMBL.Id")

df_dc_factors_other <- other_factor_datasets %>%
  reduce(full_join, by = c("Protein.Uniprot.Accession", "Gene.ENSEMBL.Id", "Gene.Symbol"),
         relationship = "many-to-one") %>%
  filter(!(is.na(Protein.Uniprot.Accession) & is.na(Gene.Symbol))) %>%
  select(-Gene.ENSEMBL.Id) %>%
  group_by(Protein.Uniprot.Accession, Gene.Symbol) %>%
  summarize_if(is.numeric, ~first(na.omit(.))) %>%
  ungroup() %>%
  mutate_at(c("Protein-Protein Interactions", "Protein Complexes (CORUM)"), ~replace_na(., 0))

df_dc_factors <- df_dc_factors_ptm %>%
  full_join(y = df_dc_factors_other, by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            relationship = "many-to-one")

# === Quality Control ===
# Check for unmatched and ambiguous rows
ambig <- df_dc_factors_other %>%
  add_count(Gene.Symbol, Protein.Uniprot.Accession) %>%
  filter(n > 1)

## Determine amount of missing data
df_dc_factors_values <- df_dc_factors %>% select(where(is.numeric))
na_counts <- colSums(is.na(df_dc_factors_values))
total_na_count <- sum(na_counts)
na_share <- total_na_count / (nrow(df_dc_factors_values) * ncol(df_dc_factors_values))

## Check correlation between factors
png(here(plots_dir, "factors_correlation.png"), width = 300, height = 300, units = "mm", res = 200)
df_dc_factors %>%
  select(all_of(dc_factor_cols)) %>%
  plot_correlation()
dev.off()

# === Write combined factors to disk ===
write_parquet(df_dc_factors, here(output_data_dir, 'dosage_compensation_factors.parquet'),
              version = "2.6")