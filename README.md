# FCA bone marrow study
## [Jardine and Webb et al., 2021, 'Blood and immune development in human fetal bone marrow and in Down syndrome' (Nature, 2021)](https://www.nature.com/articles/s41586-021-03929-x)

This repository contains all analysis scripts used to explore single cell datasets for: **Jardine and Webb et al., 2021, 'Blood and immune development in human fetal bone marrow and in Down syndrome' (in revision)**.

## Single cell datasets 

Single cell datasets used in this study include: 10X data (YS, FL, FBM, DS FBM, thymus, ABM, CB), SS2 data (FBM), and BCR-/TCR-enriched VDJ data (FBM) and CITE-seq data (CD34+ FBM/FL/CB and FBM total). 

#### Novel data accessibility
There are no restrictions on data availability for novel data presented in this study. FASTQ and raw count matrices for DS and non-DS FBM droplet-based scRNA-seq data are deposited at EMBL-EBI ArrayExpress and ENA, with accession codes as follows: [E-MTAB-9389](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9389/) (DS and non-DS FBM), [E-MTAB-10042](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10042/) (DS FBM) and [ERP125305](https://www.ebi.ac.uk/ena/browser/view/PRJEB41514) (non-DS FBM). FASTQ and raw count matrices for all other novel data in this study are deposited at EMBL-EBI ArrayExpress and GEO with accession codes [E-MTAB-9801](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9801/) (FBM Smart-seq2 scRNA-seq); [E-MTAB-9389](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9389/) (BCR-/TCR-enriched VDJ FBM scRNA-seq- FASTQs only); [GSE166895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166895) (CD34+ FBM, FL and CB CITE-seq) and [GSE166895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166895) (FBM total CITE-seq). The following data are also available to download as Scanpy h5ad objects with transformed counts via our interactive webportal: https://fbm.cellatlas.io/: i) DS FBM scRNA-seq, ii) non-DS FBM scRNA-seq, iii) CD34+ FBM, FL and CB CITE-seq, iv) FBM total CITE-seq. Our webportal also contains a searchable database of genes implicated in inherited blood and immune cell disorders. All source data are available in the accompanying source data file, unless manuscript or figure legend refers to a Supplementary Table. 

Metadata for single cell datasets described above are provided in manuscript Supplementary Tables (see manuscript Supplementary Information guide), with an overview provided in Supplementary Table 1. Source data for graphs in main and extended figures are provided as excel files. These include Fig. 1c; Fig. 2a,c,e; Fig. 3c,d; Extended Data Fig. 3c; Extended Data Fig. 4b,d,j; Extended Data Fig. 5g; Extended Data Fig. 6d,e,f; Extended Data Fig. 7c; Extended Data Fig. 8d; Extended Data Fig. 9g.   

#### External data accessibility
External datasets incorporated in this study include: i) Human FL and YS scRNA-seq data (EMBL-EBI ArrayExpress accession: [E-MTAB-7407](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7407/)); ii) Human blood monocyte-DC scRNA-seq data (NCBI GEO accession: [GSE94820](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820)); iii) Mouse BM scRNA-seq data (NCBI GEO accession: [GSE122467](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122467)); iv) Fetal and pediatric thymus scRNA-seq data (EMBL-EBI ArrayExpress accession: [E-MTAB-8581](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8581/)); v) Adult BM and CB scRNA-seq data from the Human Cell Atlas Data Coordination Portal [‘Census of Immune Cells’ project](https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79). At the time of submission, there are no known accessibility restrictions on these external datasets. 

YS - yolk sac. FL - fetal liver. FBM - fetal bone marrow. ABM - adult bone marrow. CB - cord blood. DS - Down syndrome.

## Code availability

Cellranger count matrix files for FBM 10X data were loaded into one object as described in the [Haniffa Lab Github](https://github.com/haniffalab/FCA_liver) for downstream analysis. 

Generalisable analysis scripts are saved in the 'pipelines' directory. Please refer to pipelines/readme for further information on the methods used for each pipeline. The figure panels created through use of each pipeline are detailed in readme.

Custom scripts for each figure (that do not fall under 'generalisable' pipeline scripts) are saved in the 'figures' directory, where you will see all analysis performed on raw count matrix files in order to produce each figure panel for manuscript, with exception of those detailed in 'pipelines' directory. Please refer to methods section of fetal BM manuscript for further information on methods used in each custom script. 

The authors of custom scripts (whether novel approaches or adapted from published workflows) are noted by initials appended to script filename as follows: [Simone Webb](https://github.com/simonewebb) (SW), [Issac Goh](https://github.com/Issacgoh) (IG), [Mariana Quiroga Londono](https://github.com/marianaql) (MQL), [Gary Reynolds](https://github.com/greynolds81) (GR), [Emma Dann](https://github.com/emdann) (ED), [Iwo Kucinski](https://github.com/Iwo-K) (IK).
