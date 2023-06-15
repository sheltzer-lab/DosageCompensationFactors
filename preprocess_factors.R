library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(openxlsx)
library(readxl)
library(assertr)

here::i_am("DosageCompensationFactors.Rproj")

external_data_dir <- here("Data", "External")
factor_data_dir <- here(external_data_dir, "Factors")
phosphositeplus_data_dir <- here(factor_data_dir, "PhosphoSitePlus")

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


# === Summarize Factor Datasets ===
dc_factors <- c(
  "Protein Complex (CORUM)", "Protein-Protein Interactions", "3'-UTR Length", "5'-UTR Length", "Protein Half-Life",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites", "Methylation Sites", "Acetylation Sites", "Regulatory Sites",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate", "Translation Rate", "Protein Length", "mRNA Length"
)

## Get UniProt Accession from ENSEMBL ID
ensembl_to_uniprot <- function(df, ensembl_column = "Gene.ENSEMBL.Id") {
  uniprot <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                   keys = df[[ensembl_column]],
                                   column = 'UNIPROT',
                                   keytype = 'ENSEMBL',
                                   multiVals = 'first')

  df <- df %>%
    mutate(MERGED_ACC_ID = names(uniprot),
           Protein.Uniprot.Accession = uniprot) %>%
    assertr::verify(MERGED_ACC_ID == .data[[ensembl_column]]) %>%      # Sanity check
    select(-MERGED_ACC_ID)

  return(df)
}

ensembl_to_genesymbol <- function(df, ensembl_column = "Gene.ENSEMBL.Id") {
  symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                   keys = df[[ensembl_column]],
                                   column = 'SYMBOL',
                                   keytype = 'ENSEMBL',
                                   multiVals = 'first')

  df <- df %>%
    mutate(MERGED_ACC_ID = names(symbols),
           Gene.Symbol = symbols) %>%
    assertr::verify(MERGED_ACC_ID == .data[[ensembl_column]]) %>%      # Sanity check
    select(-MERGED_ACC_ID)

  return(df)
}

genesymbol_to_ensembl <- function(df, symbol_column = "Gene.Symbol") {
  gene_IDs <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                   keys = df[[symbol_column]],
                                   column = 'ENSEMBL',
                                   keytype = 'SYMBOL',
                                   multiVals = 'first')

  df <- df %>%
    mutate(MERGED_ACC_ID = names(gene_IDs),
           Gene.ENSEMBL.Id = gene_IDs) %>%
    assertr::verify(MERGED_ACC_ID == .data[[symbol_column]]) %>%      # Sanity check
    select(-MERGED_ACC_ID)

  return(df)
}

uniprot_to_ensembl <- function(df, acc_column = "Protein.Uniprot.Accession") {
  gene_IDs <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                    keys = df[[acc_column]],
                                    column = 'ENSEMBL',
                                    keytype = 'UNIPROT',
                                    multiVals = 'first')

  df <- df %>%
    mutate(MERGED_ACC_ID = names(gene_IDs),
           Gene.ENSEMBL.Id = gene_IDs) %>%
    assertr::verify(MERGED_ACC_ID == .data[[acc_column]]) %>%      # Sanity check
    select(-MERGED_ACC_ID)

  return(df)
}


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
  ensembl_to_uniprot()

## Prepare Proten-Protein-Interaction data
hippie_filtered <- hippie %>%
  select("Gene Name Interactor A", "Gene Name Interactor B", "Confidence Value") %>%
  rename(Gene.Symbol = "Gene Name Interactor A") %>%
  filter(`Confidence Value` > 0.6) %>%
  group_by(Gene.Symbol) %>%
  count(name = "Protein-Protein Interactions") %>%
  # ToDo: Code breaks here
  genesymbol_to_ensembl() %>%
  ensembl_to_uniprot()

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
  genesymbol_to_ensembl() %>%
  ensembl_to_uniprot()

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
  summarize(`3'-UTR Length` = list(`3'-UTR Length`),
            `5'-UTR Length` = list(`5'-UTR Length`)) %>%
  ungroup()

# === Combine Factor Datasets ===

# Add phosphorylation site counts
df_dc_factors <- count_sites(phospho_sites, colname = "Phosphorylation Sites") %>%
  mutate(`Phosphorylation Sites` = replace_na(`Phosphorylation Sites`, 0)) %>%
  # Add ubiquitination site counts
  full_join(y = count_sites(ubi_sites, colname = "Ubiquitination Sites"),
            by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Ubiquitination Sites` = replace_na(`Ubiquitination Sites`, 0)) %>%
    # Add sumorylation site counts
  full_join(y = count_sites(sumo_sites, colname = "Sumoylation Sites"),
            by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Sumoylation Sites` = replace_na(`Sumoylation Sites`, 0)) %>%
  # Add methylation site counts
  full_join(y = count_sites(methyl_sites, colname = "Methylation Sites"),
            by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Methylation Sites` = replace_na(`Methylation Sites`, 0)) %>%
  # Add acetylation site counts
  full_join(y = count_sites(ace_sites, colname = "Acetylation Sites"),
            by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Acetylation Sites` = replace_na(`Acetylation Sites`, 0)) %>%
  # Add regulatory site counts
  full_join(y = count_sites(reg_sites, colname = "Regulatory Sites"),
            by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Regulatory Sites` = replace_na(`Regulatory Sites`, 0)) %>%
  # Transform Uniprot Accessions to ENSEMBL IDs
  uniprot_to_ensembl() %>%
  # Add central dogma rates data
  full_join(y = df_rates, by = c("Protein.Uniprot.Accession", "Gene.ENSEMBL.Id"),
            na_matches = "never", relationship = "many-to-one") %>%
  # Add protein-protein-interactions data
  # ToDo: Merge by ENSEMBL Id
  full_join(y = hippie_filtered, by = "Gene.Symbol",
            na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Protein-Protein Interactions` = replace_na(`Protein-Protein Interactions`, 0)) %>%
  # Add Protein Half-Life data
  full_join(y = half_life_avg, by = "Gene.Symbol",
            na_matches = "never", relationship = "many-to-one")