# Dataset References

## Protein Expression / Proteomics

### Broad DepMap

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=Proteomics&filename=protein_quant_current_normalized.csv

**Location:** `./Data/External/Expression/Broad-DepMap-Proteomics.csv`

#### Source Notes

From Proteomics
protein_quant_current_normalized.csv

This is the core table of normalized protein expression data for all cell lines

Proteins: 12399
Cell Lines: 375
Primary Diseases: 24
Lineages: 27
Source: Broad Institute

Quantitative profiling of thousands of proteins by mass spectrometry across 375 cell lines from the Gygi lab. The normalized protein quantitation data from https://gygi.med.harvard.edu/publications/ccle was imported into the portal. For more information about this dataset released in [Nusinow et al., 2020](https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0020), please see [this supplemental guide](https://www.biorxiv.org/content/10.1101/2020.02.03.932384v1).

#### Citation
David P. Nusinow, John Szpyt, Mahmoud Ghandi, Christopher M. Rose, E. Robert McDonald III, Marian Kalocsay, Judit Jané-Valbuena, Ellen Gelfand, Devin K. Schweppe, Mark Jedrychowski, Javad Golji, Dale A. Porter, Tomas Rejtar, Y. Karen Wang, Gregory V. Kryukov, Frank Stegmeier, Brian K. Erickson, Levi A. Garraway, William R. Sellers, Steven P. Gygi (2020). Quantitative Proteomics of the Cancer Cell Line Encyclopedia. Cell 180, 2. https://doi.org/10.1016/j.cell.2019.12.023

### Sanger ProCan

**URL:** https://cellmodelpassports.sanger.ac.uk/downloads

**Location:** `./Data/External/Expression/ProCan-DepMapSanger_protein_matrix_6692_averaged.txt`

#### Citation
Goncalves, Emanuel; Poulos, Rebecca C.; Cai, Zhaoxiang; Barthorpe, Syd; Manda, Srikanth; Lucas, Natasha; et al. (2022). Pan-cancer proteomic map of 949 human cell lines. figshare. Dataset. https://doi.org/10.6084/m9.figshare.19345397.v1

### CPTAC

**URL:** https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer

**Location:** `./Data/External/Expression/CPTAC/Proteome_BCM_GENCODE_v34_harmonized_v1`

#### Source Notes

File:  Proteome_BCM_GENCODE_v34_harmonized_v1.zip

Description:  CPTAC Pan-Cancer Proteome data processed by the University of Michigan team's pipeline, and then post-processed by the Baylor College of Medicine's pipeline. Details can be found in the STAR Methods of 'Proteogenomic Data and Resources for Pan-Cancer Analysis' (i.e., 'BCM pipeline for pan-cancer multi-omics data harmonization').

Cohorts:  BRCA, ccRCC, COAD, GBM, HGSC, HNSCC, LSCC, LUAD, PDAC, UCEC

#### Citation
Li, Y., Dou, Y., da Veiga Leprevost, F., Geffen, Y., et al. (2023). Proteogenomic data and resources for pan-cancer analysis, Cancer Cell, 41, 1397-1406. https://doi.org/10.1016/j.ccell.2023.06.009

### Lab-Engineered Cell Lines

#### P0211

Published as supplemental data to this project.

* **Cell Lines:** RPE1 p53-/-, RM13, Rtr13
* **Location:** `./Data/P0211/`

#### Chunduri

* **Cell Lines:** RPE1 p53 KO, RM X, RM 10;18, RM 13, RPE1 p53 KD, RM 19p
* **URL:** https://www.ebi.ac.uk/pride/archive/projects/PXD018440
* **DOI:** https://doi.org/10.1038/s41467-021-25288-x
* **Citation:** Chunduri, N.K., Menges, P., Zhang, X. et al. Systems approaches identify the consequences of monosomy in somatic human cells. Nat Commun 12, 5576 (2021). https://doi.org/10.1038/s41467-021-25288-x
* **Location:**  `./Data/External/Expression/Chunduri/proteinGroups.txt`

## Copy Number

### Arm-Level CNAs

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=CCLE%20arm-level%20CNAs&filename=arm_call_scores.csv

**Location:** `./Data/External/CopyNumber/DepMap/Arm-level_CNAs.csv`

#### Source Notes

From CCLE arm-level CNAs
arm_call_scores.csv

Arm-level CNA calls scores from Cohen-Sharir_Table_S1.xlsx, formatted as a matrix of scores versus cell lines.

Each chromosome arm is called as having copy gain (+1), loss (-1), or no change/no call (0) relative to background ploidy. Calls are derived from ABSOLUTE copy number profiles published in Ghandi et al., (Nature, 2019), by computing segment length-weighted median modal copy number for each chromosome arm, and comparing these with background ploidy estimates.

Source: Broad Institute

#### Citation
Cohen-Sharir, Yael, et al. "Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition." Nature (2021): 1-6.

### Aneuploidy Scores

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=CCLE%20arm-level%20CNAs&filename=aneuploidy_scores.csv

**Location:** `./Data/External/CopyNumber/DepMap/Aneuploidy.csv`

#### Source Notes

From DepMap Public 23Q2
aneuploidy_scores.csv

Aneuploidy scores from Cohen-Sharir_Table_S1.xlsx and Ploidy and Genome doubling scores from Ghandi et al., 2019 CCLE_ABSOLUTE_combined_20181227.xlsx. Formatted as a matrix of scores versus cell line.

Source: Broad Institute

#### Citation
Ghandi, Mahmoud, et al. “Next-Generation Characterization of the Cancer Cell Line Encyclopedia.” Nature 569, no. 7757 (May 2019): 503–8. https://doi.org/10.1038/s41586-019-1186-3.

### Relative Gene Copy Number

**URL:** https://depmap.org/portal/download/custom/

**Location:** `./Data/External/CopyNumber/DepMap/Copy_Number_Public_23Q2.csv`

#### Source Notes

From DepMap Public 23Q2
OmicsCNGene.csv

Gene-level copy number data that is log2 transformed with a pseudo-count of 1; log2(CN ratio + 1) . Inferred from WGS, WES or SNP array depending on the availability of the data type. Values are calculated by mapping genes onto the segment level calls and computing a weighted average along the genomic coordinate. Additional copy number datasets are available for download as part of the full DepMap Data Release. More information on the DepMap Omics processing pipeline is available at https://github.com/broadinstitute/depmap_omics.

Genes: 25368,
Cell Lines: 1814,
Primary Diseases: 82,
Lineages: 31,
Source: Broad Institute

#### Citation
DepMap, Broad Institute. “DepMap 23Q2 Public.” figshare, 2023. https://doi.org/10.6084/M9.FIGSHARE.22765112.V4.

### Absolute Gene Copy Number

**URL:** https://depmap.org/portal/download/custom/

**Location:** `./Data/External/CopyNumber/DepMap/Copy_Number_(Absolute).csv`

#### Citation
DepMap, Broad Institute. “DepMap 23Q2 Public.” figshare, 2023. https://doi.org/10.6084/M9.FIGSHARE.22765112.V4.

### CPTAC

Downloaded using the Genomic Data Commons Downloads API. 
For more information see https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/. 

#### Citation
Li, Y., Dou, Y., da Veiga Leprevost, F., Geffen, Y., et al. (2023). Proteogenomic data and resources for pan-cancer analysis, Cancer Cell, 41, 1397-1406. https://doi.org/10.1016/j.ccell.2023.06.009


## Metadata

### Broad DepMap

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap+Public+24Q2&filename=Model.csv

**Location:** `./Data/External/CopyNumber/DepMap/Model.csv`

#### Source Notes

From DepMap Public 23Q2
Model.csv

Metadata describing all cancer models/cell lines which are referenced by a dataset contained within the DepMap portal. Columns: - ModelID: Static primary key assigned by DepMap to each cell model - PatientID: Indicates relationships of patients from which a model was derived - CellLineName: Original cell line name, including punctuation - StrippedCellLineName: Cell line name with alphanumeric characters only - Age: Age of tissue donor at time of sample collection, if known - SourceType: Source of the cell model used by DepMap - SangerModelID: The corresponding Sanger Model ID for this model - RRID: Cellosaurus research resource identifier - DepmapModelType: Internal model type used by DepMap project (Identical to OncotreeCodes but with additional codes added for different types of 'normal' tissue) - AgeCategory: Classification of age of patient into 'Adult', 'Unknown', 'Pediatric' or 'Fetus' - GrowthPattern: Typical growth pattern of the cell model when it was onboarded - LegacyMolecularSubtype: Historical molecular subtype annotations (Annotated prior to adoption of Oncotree) - PrimaryOrMetastasis: Indicates whether tissue sample is from primary or metastatic site - SampleCollectionSite: Tissue collection site - Sex: Sex of tissue donor, if known - SourceDetail: Details of the source of the cell model used by DepMap - LegacySubSubtype: Historical sub-subtype annotations (Annotated prior to adoption of Oncotree) - CatalogNumber: The catalog number for lines obtained from a commercial vendor - CCLEName: Previous naming system that used the stripped cell line name followed by the lineage; no longer assigned to new cell lines - COSMICID: Cell line ID used in Cosmic cancer database - PublicComments: Free text notes about the Model - WTSIMasterCellID: ID of corresponding record in Sanger Drug dataset - EngineeredModel: Indicates a genetically modified cell model - TreatmentStatus: Status of patient treatment at the time the model was derived from tissue. Either 'Unknown', 'Post-treatment' or 'Pre-treatment' - OnboardedMedia: The media that was initially used to expand cell lines. - PlateCoating: The type of coating on a plate when a cell model is grown, if any - OncotreeCode: Oncotree classification of disease (Nomenclature taken from http://oncotree.mskcc.org) - OncotreeSubtype: Intermediate Oncotree classification (More specific than OncotreeLineage, but often less specific than OncotreeCode) - OncotreePrimaryDisease: The disease name corresponding to OncotreeCode - OncotreeLineage: Name of top level tissue from Oncotree nomenclature

Source: Broad Institute

#### Citation
DepMap, Broad Institute. “DepMap 23Q2 Public.” figshare, 2023. https://doi.org/10.6084/M9.FIGSHARE.22765112.V4.

### Sanger ProCan

**URL:** https://cellmodelpassports.sanger.ac.uk/downloads

**Location:** `./Data/External/CopyNumber/ProCan/model_list_20240110.csv`

### CPTAC

**URL:** https://proteomic.datacommons.cancer.gov/pdc/browse

**Location:** `./Data/External/PDC_biospecimen_manifest.tsv`

How to Download:
1. Visit URL
2. Select "Proteome" in filter "General" -> "Analytical Fraction"
3. Open the "Biospecimens" page
4. Select all entries on all pages
5. Click on "Export Biospecimen Manifest (TSV)"

## Dosage Compensation Factors

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

### Haploinsufficiency & Triplosensitivity

**URL:** https://www.sciencedirect.com/science/article/pii/S0092867422007887?via%3Dihub#mmc7

**Location:** `./Data/External/Factors/1-s2.0-S0092867422007887-mmc7.xlsx`

**DOI:** [10.1016/j.cell.2022.06.036](https://doi.org/10.1016/j.cell.2022.06.036)

Table S7. Haploinsufficiency and triplosensitivity predictions for all autosomal protein-coding genes, related to Figure 6.

### Random Allelic Expression

**URL:** https://ars.els-cdn.com/content/image/1-s2.0-S2211124722018460-mmc3.xlsx

**Location:** `./Data/External/Factors/random_allelic_expression.xlsx`

#### Citation

Kravitz, Stephanie N., Elliott Ferris, Michael I. Love, Alun Thomas, Aaron R. Quinlan, and Christopher Gregg. “Random Allelic Expression in the Adult Human Body.” Cell Reports 42, no. 1 (January 2023): 111945. https://doi.org/10.1016/j.celrep.2022.111945.

### Transcription Factors

Obtained DoRothEA regulons using the `decoupleR` package in R.

#### Citation
Garcia-Alonso L, Holland C, Ibrahim M, Turei D, Saez-Rodriguez J (2019). “Benchmark and integration of resources for the estimation of human transcription factor activities.” Genome Research. doi:10.1101/gr.240663.118.

Badia-i-Mompel P, Santiago JV, Braunger J, Geiss C, Dimitrov D, Müller-Dott S, Taus P, Dugourd A, Holland CH, Flores ROR, Saez-Rodriguez J (2022). “decoupleR: ensemble of computational methods to infer biological activities from omics data.” Bioinformatics Advances. doi:10.1093/bioadv/vbac016.

Müller-Dott S, Tsirvouli E, Vázquez M, Flores ROR, Badia-i-Mompel P, Fallegger R, Lægreid A, Saez-Rodriguez J (2023). “Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities.” bioRxiv. doi:10.1101/2023.03.30.534849.


### CORUM

Version: 210512

### PhosphoSitePlus

**URL:** https://www.phosphosite.org/staticDownloads.action

Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res. 2015 43:D512-20.


## Screens

### Drug Screens

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=PRISM%20Repurposing%20Public%2023Q2&filename=Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv

**Location:** `./Data/External/Screens/PRISM_Repurposing_Public_23Q2.csv`

From PRISM Repurposing Public 23Q2
Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv

The final LFC values in Repurposing_Public_23Q2_LFC_collapsed.csv table are cast into a matrix with rows corresponding to individual treatments and columns are depmap_id's. This final matrix also contains the compounds from the original PRISM Repurposing Primary screen for convenience. Please use Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv to get the meta-data for each row.

Source: Broad Institute

---

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=PRISM%20Repurposing%20Public%2023Q2&filename=Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv

**Location:** `./Data/External/Screens/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv`

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

### CRISPR-KO Screens

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2023Q2&filename=CRISPRGeneEffect.csv

**Location:** `./Data/External/Screens/CRISPRGeneEffect.csv`

From DepMap Public 23Q2
CRISPRGeneEffect.csv

Post-Chronos Gene effect estimates for all models, integrated using Harmonia. - Columns: Gene - Rows: ModelID

Source: Broad Institute

---

**URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2023Q2&filename=CRISPRGeneDependency.csv

**Location:** `./Data/External/Screens/CRISPRGeneDependency.csv`

From DepMap Public 23Q2
CRISPRGeneDependency.csv

Post-Chronos Gene dependency probability estimates for all models in the integrated gene effect. - Columns: Gene - Rows: ModelID

Source: Broad Institute

#### Citations

Robin M. Meyers, Jordan G. Bryan, James M. McFarland, Barbara A. Weir, ... David E. Root, William C. Hahn, Aviad Tsherniak. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nature Genetics 2017 October 49:1779-1784. doi:10.1038/ng.3984

Dempster, J. M., Rossen, J., Kazachkova, M., Pan, J., Kugener, G., Root, D. E., & Tsherniak, A. (2019). Extracting Biological Insights from the Project Achilles Genome-Scale CRISPR Screens in Cancer Cell Lines. BioRxiv, 720243.

Dempster, J. M., Boyle, I., Vazquez, F., Root, D., Boehm, J. S., Hahn, W. C., Tsherniak, A., & McFarland, J. M. (2021). Chronos: a CRISPR cell population dynamics model. BioRxiv, 432728.

Pacini, C., Dempster, J. M., Boyle, I., Goncalves, E., Najgebauer, H., Karakoc, E., van der Meer, D., Barthorpe, A., Lightfoot, H., Jaaks, P., McFarland, J. M., Garnett, M. J., Tsherniak, A., & Iorio, F. (2021) Integrated cross-study datasets of genetic dependencies in cancer. Nature Commmunications, 12-1661.

### GDSC Cell Line Growth Rates 

**URL:** https://cog.sanger.ac.uk/cmp/download/growth_rate_20220907.csv

**Location:** `./Data/External/Screens/growth_rate_20220907.csv`

**Source Notes:** https://depmap.sanger.ac.uk/documentation/datasets/growth-rate/

### CRISPR-Inferred Model Growth Rate

**URL:** https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2&filename=CRISPRInferredModelGrowthRate.csv

**Location:** `./Data/External/Screens/CRISPRInferredModelGrowthRate.csv`

#### Source Notes

From DepMap Public 23Q2
CRISPRInferredModelGrowthRate.csv

Post-Chronos The estimates for the growth rate of all models in the different libraries-screen types, computed from the Chronos runs. Columns: - ScreenID - Achilles-Avana-2D - Project-Score-KY

Source: Broad Institute

#### Citation

Robin M. Meyers, Jordan G. Bryan, James M. McFarland, Barbara A. Weir, ... David E. Root, William C. Hahn, Aviad Tsherniak. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nature Genetics 2017 October 49:1779-1784. doi:10.1038/ng.3984

Dempster, J. M., Rossen, J., Kazachkova, M., Pan, J., Kugener, G., Root, D. E., & Tsherniak, A. (2019). Extracting Biological Insights from the Project Achilles Genome-Scale CRISPR Screens in Cancer Cell Lines. BioRxiv, 720243.

Dempster, J. M., Boyle, I., Vazquez, F., Root, D., Boehm, J. S., Hahn, W. C., Tsherniak, A., & McFarland, J. M. (2021). Chronos: a CRISPR cell population dynamics model. BioRxiv, 432728.

Pacini, C., Dempster, J. M., Boyle, I., Goncalves, E., Najgebauer, H., Karakoc, E., van der Meer, D., Barthorpe, A., Lightfoot, H., Jaaks, P., McFarland, J. M., Garnett, M. J., Tsherniak, A., & Iorio, F. (2021) Integrated cross-study datasets of genetic dependencies in cancer. Nature Commmunications, 12-1661.

## Other

### UniProt ID Mapping

https://www.uniprot.org/help/id_mapping
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

Data Source URL: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz


### Proteomics Reproducibility - Protein-Protein Correlation - Aggregated Ranks

* DOI: 10.1016/j.crmeth.2022.100288
* File: 1-s2.0-S2667237522001709-mmc3.xlsx 
  * Table S2. B.: Computed correlation between experimental proteomic replicates of ovarian and colon tumor samples and cancer cell lines using the standard pipeline and the aggregated normalized ranks, related to Figures 1 and S2 and STAR Methods. Computed protein-protein reproducibility ranks

### Cancer Driver Genes

* **Source:** OncoKB
* **URL:** https://www.oncokb.org/cancer-genes
* **Location:** `./Data/External/cancerGeneList.tsv`

#### Citation

* https://doi.org/10.1158/2159-8290.CD-23-0467
* https://doi.org/10.1200/PO.17.00011

### Damaging Mutations (DepMap)

* **Source:** DepMap
* **URL:** https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap+Public+24Q2&filename=OmicsSomaticMutationsMatrixDamaging.csv
* **Location:** `./Data/External/OmicsSomaticMutationsMatrixDamaging.csv`
