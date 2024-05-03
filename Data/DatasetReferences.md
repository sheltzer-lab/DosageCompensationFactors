# Dataset References

## Broad DepMap

URL: https://depmap.org/portal/download/custom/

---

From CCLE arm-level CNAs
arm_call_scores.csv

Arm-level CNA calls scores from Cohen-Sharir_Table_S1.xlsx, formatted as a matrix of scores versus cell lines.

Each chromosome arm is called as having copy gain (+1), loss (-1), or no change/no call (0) relative to background ploidy. Calls are derived from ABSOLUTE copy number profiles published in Ghandi et al., (Nature, 2019), by computing segment length-weighted median modal copy number for each chromosome arm, and comparing these with background ploidy estimates.

#### Data Usage

##### Please cite the following when using these data

Cohen-Sharir, Yael, et al. "Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition." Nature (2021): 1-6.

Source: Broad Institute

---

From DepMap Public 23Q2
aneuploidy_scores.csv

Aneuploidy scores from Cohen-Sharir_Table_S1.xlsx and Ploidy and Genome doubling scores from Ghandi et al., 2019 CCLE_ABSOLUTE_combined_20181227.xlsx. Formatted as a matrix of scores versus cell line.

Source: Broad Institute

---

From DepMap Public 23Q2
OmicsCNGene.csv

Gene-level copy number data that is log2 transformed with a pseudo-count of 1; log2(CN ratio + 1) . Inferred from WGS, WES or SNP array depending on the availability of the data type. Values are calculated by mapping genes onto the segment level calls and computing a weighted average along the genomic coordinate. Additional copy number datasets are available for download as part of the full DepMap Data Release. More information on the DepMap Omics processing pipeline is available at https://github.com/broadinstitute/depmap_omics.

Genes: 25368,
Cell Lines: 1814,
Primary Diseases: 82,
Lineages: 31,
Source: Broad Institute

---

From Proteomics
protein_quant_current_normalized.csv

This is the core table of normalized protein expression data for all cell lines

Proteins: 12399
Cell Lines: 375
Primary Diseases: 24
Lineages: 27
Source: Broad Institute

Proteomics

Quantitative profiling of thousands of proteins by mass spectrometry across 375 cell lines from the Gygi lab. The normalized protein quantitation data from https://gygi.med.harvard.edu/publications/ccle was imported into the portal. For more information about this dataset released in Nusinow et al., 2020, please see this supplemental guide.

Data Usage
Please cite the following when using these data

David P. Nusinow, John Szpyt, Mahmoud Ghandi, Christopher M. Rose, E. Robert McDonald III, Marian Kalocsay, Judit Jané-Valbuena, Ellen Gelfand, Devin K. Schweppe, Mark Jedrychowski, Javad Golji, Dale A. Porter, Tomas Rejtar, Y. Karen Wang, Gregory V. Kryukov, Frank Stegmeier, Brian K. Erickson, Levi A. Garraway, William R. Sellers, Steven P. Gygi (2020). Quantitative Proteomics of the Cancer Cell Line Encyclopedia. Cell 180, 2. https://doi.org/10.1016/j.cell.2019.12.023

---

From PRISM Repurposing Public 23Q2
Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv

The final LFC values in Repurposing_Public_23Q2_LFC_collapsed.csv table are cast into a matrix with rows corresponding to individual treatments and columns are depmap_id's. This final matrix also contains the compounds from the original PRISM Repurposing Primary screen for convenience. Please use Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv to get the meta-data for each row.

Source: Broad Institute

---

From PRISM Repurposing Public 23Q2

Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv

Metadata for the perturbations represented by each row of  Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv. Columns of  this table are given below:

- IDs: Unique identifier to join with the rows of
- Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv
- screen: Identifier for in which Repurposing screen this particular  profile is created.
- dose: Dose of the treatment (in uM).  Drug.Name: Human read-able name for the drug.
- repurposing_target: Annotated target of the drug.
- MOA: Annotated mechanism of action for the drug.
- Synonyms: Known synonyms of the drug name.

Source: Broad Institute

---

From DepMap Public 23Q2
CRISPRGeneEffect.csv

Post-Chronos Gene effect estimates for all models, integrated using Harmonia. - Columns: Gene - Rows: ModelID

Source: Broad Institute

---

From DepMap Public 23Q2

CRISPRGeneDependency.csv

Post-Chronos Gene dependency probability estimates for all models in the integrated gene effect. - Columns: Gene - Rows: ModelID

Source: Broad Institute

## Sanger ProCan

Goncalves, Emanuel; Poulos, Rebecca C.; Cai, Zhaoxiang; Barthorpe, Syd; Manda, Srikanth; Lucas, Natasha; et al. (2022). Pan-cancer proteomic map of 949 human cell lines. figshare. Dataset. https://doi.org/10.6084/m9.figshare.19345397.v1

## PhosphoSitePlus

Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res. 2015 43:D512-20.

URL: https://www.phosphosite.org/staticDownloads.action

## Sanger DepMap

URL: https://cellmodelpassports.sanger.ac.uk/downloads

## CORUM

Version: 210512

## UniProt ID Mapping

https://www.uniprot.org/help/id_mapping
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

Data Source URL: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

## Other

### HIPPIE

URL: http://cbdm-01.zdv.uni-mainz.de/\~mschaefer/hippie/download.php, version: v2.3

### Central Dogma Rates

Hausser, Jean (2019), “Central dogma rates and the trade-off between precision and economy”, Mendeley Data, V1, doi: 10.17632/2vbrg3w4p3.1

### NCBI RefSeq (UTR data)

URL: https://genome.ucsc.edu/cgi-bin/hgTables, genome: Human, assembly: Dec. 2013 (GRCh38/hg38), track: NCBI RefSeq

### Protein Half-Life data

DOI: 10.1038/s41467-018-03106-1, Supplementary Table 2

### MobiDB (Protein Disorder, etc.)

URL: https://mobidb.bio.unipd.it/browse?proteome=UP000005640&limit=10&skip=0&projection=acc,name,organism,reviewed,prediction-disorder-mobidb_lite.content_fraction,gene,length

Version: 5.0, Release: 2022_07, Accessed: 2023-06-23

### mRNA Decay Rates

DOI: 10.1101/gr.1272403, Supplementary Table 9: Decay Rates (hour^-1) for Accessions in HepG2 Experiments

### Non-Exponential Decay Deltas

DOI: 10.1016/j.cell.2016.09.015

* Supplementary Table 1: All AHA Pulse-Chase and Related Data for the Mouse NIH 3T3 Cells
* Supplementary Table 4: All AHA Pulse-Chase and Related Data for the Human RPE-1 Cells

### Aggregation Scores

DOI: 10.1016/j.celrep.2013.09.043

* Supplementary Table 2: Supersaturation Database.

### Haploinsufficiency Scores

High ranks (e.g. 0-10%) indicate a gene which is more likely to exhibit the features of haploinsufficient genes, low ranks (e.g. 90-100%) indicate a gene is more likely to NOT exhibit haploinsufficiency. The manuscript describing the generation and validation of these haploinsufficiency predictions (Huang et al., 2010) is published in PLoS Genetics. Updated predictions of haploinsufficiency can be downloaded from our data download page.

Download URL: https://www.deciphergenomics.org/about/downloads/data

FAQ URL: https://www.deciphergenomics.org/about/faqs

DOI: 10.1371/journal.pgen.1001154

### Proteomics Reproducibility - Protein-Protein Correlation - Aggregated Ranks

* DOI: 10.1016/j.crmeth.2022.100288
* File: 1-s2.0-S2667237522001709-mmc3.xlsx 
  * Table S2. B.: Computed correlation between experimental proteomic replicates of ovarian and colon tumor samples and cancer cell lines using the standard pipeline and the aggregated normalized ranks, related to Figures 1 and S2 and STAR Methods. Computed protein-protein reproducibility ranks