library(tidyr)
library(dplyr)

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