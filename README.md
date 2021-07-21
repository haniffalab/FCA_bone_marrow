# FCA bone marrow study
## Jardine and Webb et al., 2021, 'Blood and immune development in human fetal bone marrow and in Down syndrome' (in revision)

This repository contains all analysis scripts used to explore single cell datasets for: **Jardine and Webb et al., 2021, 'Blood and immune development in human fetal bone marrow and in Down syndrome' (in revision)**.

## Single cell datasets 

Single cell datasets used in this study include: 10X data (YS, FL, FBM, DS FBM, thymus, ABM, CB), SS2 data (FBM), and BCR-/TCR-enriched VDJ data (FBM) and CITE-seq data (CD34+ FBM/FL/CB and FBM total). 

#### Novel data accessibility
There are no restrictions on data availability for novel data presented in this study. The raw and processed DS/non-DS FBM GEX/BCR-/TCR-enriched droplet-based scRNA-seq data for this study are deposited at EMBL-EBI ArrayExpress, EMBL-EBI ENA, and NCBI GEO, with accession codes as follows: [E-MTAB-9389](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9389/), [E-MTAB-10042](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10042/) and [ERP125305](https://www.ebi.ac.uk/ena/browser/view/PRJEB41514). Related accession codes for this novel data in this study (including plate-based scRNA-seq and CITE-seq) are linked to the main accession at [E-MTAB-9389](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9389/) and include: [E-MTAB-9801](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9801/) for FBM Smart-seq2 scRNA-seq; [GSE166895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166895) for CD34+ (FBM, FL and CB) CITE-seq; [GSE166895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166895) for FBM (total) CITE-seq. 

Our webportal (https://fbm.cellatlas.io/) makes all novel single cell datasets generated for this study (except FBM SS2) publicly available for exploration (and download as Scanpy h5ad objects with transformed counts) with a searchable database of genes implicated in inherited blood and immune cell disorders. 

Metadata for datasets described above are provided in manuscript Supplementary Table 1. All source data are available in the accompanying source data file, unless manuscript or figure legend refers to a Supplementary Table. 

#### External data accessibility
External datasets incorporated in this study include: i)  ABM and CB scRNA-seq data (Human Cell Atlas Data Coordination Portal [‘Census of Immune Cells’ project](https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79)); ii) FL and YS scRNA-seq data were sourced from a [previous publication](https://doi.org/10.1038/s41586-019-1652-y) from the Haniffa Lab ([E-MTAB-7407](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7407/)), whose initial analysis scripts can be found on the [Haniffa Lab Github](https://github.com/haniffalab/FCA_liver); iii) [Monocyte-DC blood SS2 data](https://science.sciencemag.org/content/356/6335/eaah4573) ([GSE94820](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820) ‘GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt’ file); iv) [Mouse BM scRNA-seq data](https://doi.org/10.1038/s41556-019-0439-6) data in all files in [GSE122467](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122467)); v) [Thymus scRNA-seq data](https://science.sciencemag.org/content/367/6480/eaay3224/tab-article-info). At the time of submission, there are no known accessibility restrictions on these external datasets. All single cell datasets (except the FBM Smart-seq2 scRNA-seq) are also available to download as Scanpy h5ad objects with transformed counts via our interactive webportal: https://fbm.cellatlas.io/. 

YS - yolk sac. FL - fetal liver. FBM - fetal bone marrow. ABM - adult bone marrow. CB - cord blood. DS - Down syndrome.

## Code availability

Cellranger count matrix files for FBM 10X data were loaded into one object as described in the [Haniffa Lab Github](https://github.com/haniffalab/FCA_liver) for downstream analysis. 

Generalisable analysis scripts are saved in the 'pipelines' directory. Please refer to pipelines/readme for further information on the methods used for each pipeline. The figure panels created through use of each pipeline are detailed in readme.

Custom scripts for each figure (that do not fall under 'generalisable' pipeline scripts) are saved in the 'figures' directory, where you will see all analysis performed on raw count matrix files in order to produce each figure panel for manuscript, with exception of those detailed in 'pipelines' directory. Please refer to methods section of fetal BM manuscript for further information on methods used in each custom script. 

The authors of custom scripts (whether novel approaches or adapted from published workflows) are noted by initials appended to script filename as follows: [Simone Webb](https://github.com/simonewebb) (SW), [Issac Goh](https://github.com/Issacgoh) (IG), [Mariana Quiroga Londono](https://github.com/marianaql) (MQL), [Gary Reynolds](https://github.com/greynolds81) (GR), [Emma Dann](https://github.com/emdann) (ED), [Iwo Kucinski](https://github.com/Iwo-K) (IK).
