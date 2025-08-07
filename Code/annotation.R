library(tidyr)
library(dplyr)
library(stringr)

ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl") # GRCh38.p14

get_chromosome_arms <- function(df, mart = ensembl_mart, id_col = "Gene.Symbol", id_type = "hgnc_symbol") {
  chr_location <- biomaRt::getBM(attributes = c(id_type,
                                                'chromosome_name', 'band',
                                                'start_position', 'end_position'),
                                 filters = id_type,
                                 values = unique(df[[id_col]]),
                                 mart = ensembl_mart) %>%
    rename(ID.Temp = id_type,
           Gene.Chromosome = "chromosome_name",
           Gene.ChromosomeBand = "band",
           Gene.StartPosition = "start_position",
           Gene.EndPosition = "end_position") %>%
    mutate(Gene.ChromosomeArm = ifelse(str_detect(Gene.ChromosomeBand, "q"), "q",
                                       ifelse(str_detect(Gene.ChromosomeBand, "p"), "p",
                                              NA))) %>%
    unite("Gene.ChromosomeArm", Gene.Chromosome, Gene.ChromosomeArm,
          sep = "", remove = FALSE) %>%
    mutate(across(all_of(c("ID.Temp", "Gene.Chromosome", "Gene.ChromosomeBand", "Gene.ChromosomeArm")),
                  ~ na_if(., ""))) %>%
    drop_na() %>%
    distinct() %>%
    # Remove genes that may be located on multiple chromosomes
    add_count(ID.Temp) %>%
    filter(n == 1) %>%
    select(-n)

  df <- df %>%
    mutate(ID.Temp = .data[[id_col]]) %>% # Ensure correct column name
    left_join(y = chr_location, by = "ID.Temp",
              na_matches = "never", relationship = "many-to-one") %>%
    select(-ID.Temp)

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

# For possible ID types look into ensembl_mart$filters
# E.g. uniprot_gn_id, hgnc_symbol
id2uniprot_acc <- function(df, id_col, id_type, mart = ensembl_mart) {
  df_acc <- biomaRt::getBM(attributes = c(id_type, 'uniprot_gn_id',
                                          'uniprotswissprot'),
                           filters = id_type,
                           values = unique(df[[id_col]]),
                           mart = mart) %>%
    # Use swissprot ID if available
    mutate(uniprotswissprot = if_else(uniprotswissprot == "", NA, uniprotswissprot),
           uniprot_gn_id = if_else(uniprot_gn_id == "", NA, uniprot_gn_id)) %>%
    group_by_at(id_type) %>%
    summarize(Protein.Uniprot.Accession = if_else(any(!is.na(uniprotswissprot)),
                                                  first(uniprotswissprot, na_rm = TRUE),
                                                  first(uniprot_gn_id, na_rm = TRUE))) %>%
    rename(!!id_col := id_type) %>%
    drop_na()

  df <- df %>%
    left_join(y = df_acc, by = id_col,
              na_matches = "never", relationship = "many-to-one")
}

# For possible ID types look into ensembl_mart$filters
# E.g. uniprot_gn_id, hgnc_symbol
uniprot_acc2id <- function(df, id_col, id_type, acc_col = "Protein.Uniprot.Accession", mart = ensembl_mart) {
  df_acc <- biomaRt::getBM(attributes = c(id_type, 'uniprot_gn_id'),
                           filters = 'uniprot_gn_id',
                           values = unique(df[[acc_col]]),
                           mart = mart) %>%
    group_by_at("uniprot_gn_id") %>%
    summarize_at(id_type, ~first(.x, na_rm = TRUE)) %>%
    rename(!!id_col := id_type,
           !!acc_col := "uniprot_gn_id") %>%
    drop_na()

  df <- df %>%
    left_join(y = df_acc, by = acc_col,
              na_matches = "never", relationship = "many-to-one")
}

updateGeneSymbols <- function(df, gene_col = "Gene.Symbol") {
  # new.hgnc.table <- getCurrentHumanMap()
  updated_genes <- HGNChelper::checkGeneSymbols(unique(df[[gene_col]]), unmapped.as.na = FALSE, species = "human") %>%
    rename(!!gene_col := x)

  message("Unapproved Gene Symbols: ", updated_genes %>%
    filter(Approved == FALSE) %>%
    count())

  df %>%
    left_join(y = updated_genes, by = gene_col,
              relationship = "many-to-one", unmatched = "error") %>%
    mutate(!!gene_col := toupper(Suggested.Symbol)) %>%
    select(-Approved, -Suggested.Symbol)
}

get_oncotree_parent <- function(df_tumor_types = NULL, start_code = NULL, target_level = 3) {
  require(mskcc.oncotree)

  if (is.null(start_code) || is.na(start_code)) return(NA)

  if (is.null(df_tumor_types))
    df_tumor_types <- mskcc.oncotree::get_tumor_types()

  current <- df_tumor_types %>%
    filter(oncotree_code == start_code)

  if (nrow(current) == 0) return(NA)

  if (current$level <= target_level) return(current$oncotree_code)
  else get_oncotree_parent(df_tumor_types, current$parent, target_level)
}
