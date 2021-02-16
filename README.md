# FCA bone marrow project

This repository contains all analysis scripts used to explore single cell datasets for: Jardine and Webb et al., 2021, 'Intrinsic and extrinsic regulation of human fetal bone marrow haematopoiesis and perturbations in Down syndrome', (manuscript in revision).

## Single cell datasets 

Single cell datasets used in this study include: 10X data (YS, FL, fetal BM, DS fetal BM, thymus, ABM, CB), SS2 data (fetal BM), and BCR-/TCR-enriched VDJ data (fetal BM) and CITE-seq data (CD34+ fetal BM/FL/CB and fetal BM total).

10X data (YS and FL) were sourced from a previous publication from the Haniffa Lab titled ['Decoding fetal liver haematopoiesis'](https://doi.org/10.1038/s41586-019-1652-y), whose initial analysis scripts can be found on the [Haniffa Lab Github](https://github.com/haniffalab/FCA_liver). 10x data (thymus) were sourced from ['A cell atlas of human thymic development defines T cell repertoire formation'](https://science.sciencemag.org/content/367/6480/eaay3224). 10X data (ABM and CB)  were sourced from the publicly available [Immune Cell Atlas](https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79). See Human Cell Atlas for more details.

Multiple datasets were novel contributions from this study, including: 10x data (DS and non-DS fetal BM), SS2 data (fetal BM), and BCR-/TCR-enriched VDJ data (fetal BM) and CITE-seq data (CD34+ fetal BM/FL/CB and fetal BM total). The raw and processed DS/non-DS fetal bone marrow 10x GEX/BCR-/TCR-enriched scRNA-seq data for this study are deposited at EMBL-EBI ArrayExpress, EMBL-EBI ENA, and NCBI GEO, with accession codes as follows: E-MTAB-9389, E-MTAB-10042 and ERP125305. Related accession codes for this study (including SS2 and CITE-seq data) are linked to the main accession at E-MTAB-9389 and include: E-MTAB-9801 for fetal BM SS2; XXX for CD34+ fetal BM, fetal liver and cord blood CITE-seq; XXX  for fetal BM CITE-seq.

YS - yolk sac. FL - fetal liver. BM - fetal bone marrow. ABM - adult bone marrow. CB - cord blood. DS - Down syndrome.

## Analysis scripts

Cellranger count matrix files for FBM 10X data were loaded into one object as described in the [Haniffa Lab Github](https://github.com/haniffalab/FCA_liver) for downstream analysis. 

Generalisable analysis scripts are saved in the 'pipelines' directory. Please refer to pipelines/readme for further information on the methods used for each pipeline. The figure panels created through use of each pipeline are detailed in readme.

Custom scripts for each figure (that do not fall under 'generalisable' pipeline scripts) are saved in the 'figures' directory, where you will see all analysis performed on raw count matrix files in order to produce each figure panel for manuscript, with exception of those detailed in 'pipelines' directory. Please refer to methods section of fetal BM manuscript for further information on methods used in each custom script. 

The authors of custom scripts (whether novel approaches or adapted from published workflows) are noted by initials appended to script filename as follows: Simone Webb (SW), Gary Reynolds (GR), Issac Goh (IG), Mariana Quiroga Londono (MQL), Emma Dann (ED), Iwo Kucinski (IK).
