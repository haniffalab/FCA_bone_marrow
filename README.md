# FCA bone marrow project

This repository contains all analysis scripts used to explore single cell datasets for: Jardine and Webb et al., 2021, 'Intrinsic and extrinsic regulation of human fetal bone marrow haematopoiesis and perturbations in Down syndrome', (manuscript in revision).

## Datasets 

Datasets analysed include: 10X data (YS, FL, FBM, thymus, DS FBM, ABM, CB), FBM SS2 data, and FBM BCR-/TCR-enriched VDJ data.

YS and FL 10X data were sourced from a previous publication from the Haniffa Lab (Popescu, Botting, Stephenson et al 2019, "Decoding fetal liver haematopoiesis"), whose initial analysis scripts can be found here: https://github.com/haniffalab/FCA_liver

Thymus 10X data were sourced from Park et al 2020, "A cell atlas of human thymic development defines T cell repertoire formation".

ABM and CB 10X data were sourced from the publicly available Immune Cell Atlas: https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79

The raw and processed DS/non-DS fetal bone marrow 10x GEX/BCR-/TCR-enriched scRNA-seq data for this study are deposited at ArrayExpress or European Nucleotide Archive, with accession codes as follows: E-MTAB-9389, E-MTAB-10042 and ERP125305. Related accession codes for this study (including Smart-seq2 and CITE-seq data) are linked to the main accession at E-MTAB-9389 and include: E-MTAB-9801 for fetal BM Smart-seq2; XXX for BG fetal BM CD34+ CITE-seq / CB CD34+ CITE-seq / fetal BM MNC CITE-seq

See Human Cell Atlas for more details.

YS - yolk sac. FL - fetal liver. FBM - fetal bone marrow. ABM - adult bone marrow. CB - cord blood. DS - Downs syndrome.

## Analysis scripts

Cellranger count matrix files for FBM 10X data were loaded into one object as described in https://github.com/haniffalab/FCA_liver for downstream analysis. 

Generalisable scripts are saved in the 'pipelines' directory. Please refer to pipelines/readme file for further information on the methods used for each pipeline. The figure panels created through use of each pipeline are detailed in readme file.

Custom scripts for each figure (that do not fall under 'generalisable' pipeline scripts) are saved in the 'figures' directory, where you will see all analysis performed on raw count matrix files in order to produce each figure panel for manuscript, with exception of those detailed in 'pipelines' directory. Please refer to methods section of fetal BM manuscript for further information on methods used in each custom script. The authors of custom scripts (whether novel approaches or adapted from published workflows) are noted by initials appended to script filename as follows: Simone Webb (SW), Gary Reynolds (GR), Issac Goh (IG), Mariana Quiroga Londono (MQL), Emma Dann (ED), Iwo Kucinski (IK).
