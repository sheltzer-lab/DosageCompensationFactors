library(tidyr)
library(dplyr)

ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")

get_chromosome_arms <- function(df, mart = ensembl_mart, symbol_col = "Gene.Symbol") {
  chr_location <- biomaRt::getBM(attributes = c("hgnc_symbol",
                                                'chromosome_name', 'band',
                                                "start_position", "end_position"),
                                 filters = 'hgnc_symbol',
                                 values = unique(df[[symbol_col]]),
                                 mart = ensembl_mart) %>%
    rename(Gene.Symbol = "hgnc_symbol",
           Gene.Chromosome = "chromosome_name",
           Gene.ChromosomeBand = "band",
           Gene.StartPosition = "start_position",
           Gene.EndPosition = "end_position") %>%
    mutate(Gene.ChromosomeArm = ifelse(str_detect(Gene.ChromosomeBand, "q"), "q",
                                       ifelse(str_detect(Gene.ChromosomeBand, "p"), "p",
                                              NA))) %>%
    unite("Gene.ChromosomeArm", Gene.Chromosome, Gene.ChromosomeArm,
          sep = "", remove = FALSE) %>%
    mutate(across(all_of(c("Gene.Symbol", "Gene.Chromosome", "Gene.ChromosomeBand", "Gene.ChromosomeArm")),
                  ~ na_if(., ""))) %>%
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

## Keys/Values: UNIPROT, ENSEMBL, SYMBOL, ...
## See documentation of org.Hs.eg.db package on Bioconductor
mapIds <- function(df, key_type, value_type, key_col, value_col) {
    key_values <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                   keys = df[[key_col]],
                                   column = value_type,
                                   keytype = key_type,
                                   multiVals = 'first')
  df <- df %>%
    mutate(Key = names(key_values),
           !!value_col := key_values) %>%
    assertr::verify(Key == .data[[key_col]]) %>%      # Sanity check
    select(-Key)
  return(df)
}