# FCA bone marrow project

This repository contains all analysis scripts used to explore single cell datasets for: Jardine et al., 2020, "Intrinsic and extrinsic regulation of human fetal bone marrow haematopoiesis and perturbations in Down syndrome" (manuscript submitted).

## Datasets 

Datasets analysed include: 10X data (YS, FL, FBM, thymus, Downs FBM, ABM, CB), FBM SS2 data, and FBM BCR-/TCR-enriched VDJ data.

Yolk sac and fetal liver 10X data were sourced from a previous publication from the Haniffa Lab (Popescu, Botting, Stephenson et al 2019, "Decoding fetal liver haematopoiesis"), whose initial analysis scripts can be found here: https://github.com/haniffalab/FCA_liver

Thymus 10X data were sourced from Park et al 2020, "A cell atlas of human thymic development defines T cell repertoire formation".

Downs and non-Downs FBM 10X were novel datasets produced for this study. Downs and non-Downs FBM SS2 were novel datasets produced for this study. Healthy FBM BCR-/TCR-enriched VDJ data were novel datasets produced for this study. The non DS and DS 10X data (as FASTQ files) can be accessed at ArrayExpress (under accession code: E-MTAB-9389).

ABM and CB 10X data were sourced from the publicly available Immune Cell Atlas: https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79

See Human Cell Atlas for more details.

YS - yolk sac. FL - fetal liver. FBM - fetal bone marrow. ABM - adult bone marrow. CB - cord blood.

## Analysis scripts

Cellranger count matrix files for FBM 10X data were loaded into one object as described in https://github.com/haniffalab/FCA_liver, then cells were clustered into major groups and analysed. Raw counts and cell type annotations were uploaded to ArrayExpress (accession code: E-MTAB-9389).

Generalisable scripts are saved in the 'pipelines' directory. Please refer to pipelines/readme for further information on the methods used for each pipeline.

Custom scripts for each figure (that do not fall under 'generalisable' pipeline scripts) are saved in the 'figures' directory, where you will see all analysis performed on raw count matrix files in order to produce figures for manuscript (with methodology not detailed in readme file - please refer to methods section of fetal BM manuscript for further information). Custom scripts written by SW unless otherwise stated (with 'IG'/'GR' in title of script).

The web portal (https://developmentcellatlas.ncl.ac.uk/) files were made using the single cell analysis bundle (v1.0.0, https://github.com/haniffalab/scRNA-seq_analysis), set up on the School of Computing Science HPC Facility, Newcastle University.
