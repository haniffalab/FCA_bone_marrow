# FCA_bone_marrow

This repository contains all analysis scripts used to explore single cell datasets for: Jardine and Webb et al., 2020, "Emergent haematopoiesis in the human fetal bone marrow" (manuscript not submitted).

## Datasets 

Datasets analaysed include: 10X data (YS, FL, FBM, thymus, Down's FBM, ABM, CB), FBM SS2 data, and FBM BCR-/TCR-enriched VDJ data.

Yolk sac and fetal liver 10X data were sourced from a previous publication from the Haniffa Lab (Popescu, Botting, Stephenson et al 2019, "Decoding fetal liver haematopoiesis"), whose initial analysis scripts can be found here: https://github.com/haniffalab/FCA_liver

Thymus 10X data were sourced from Park et al 2020, "A cell atlas of human thymic development defines T cell repertoire formation".

Down's FBM and healthy FBM 10X were novel datasets produced for this study. Down's/healthy FBM SS2 data and healthy FBM BCR-/TCR-enriched VDJ data were novel datasets produced for this study. These data (as FASTQ files) can be accessed at ArrayExpress (E-MTAB-####).

ABM and CB 10X data were sourced from the publicly available Immune Cell Atlas: https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79

See Human Cell Atlas for more details.

YS - yolk sac
FL - fetal liver
FBM - fetal bone marrow
ABM - adult bone marrow
CB - cord blood

## Analysis scripts

Cellranger count matrix files were loaded into one object as described in the 'data' directory, then cells were subset into four major groups and analysed. Finally, cell type annotations and raw counts were exported and uploaded to ArrayExpress (E-MTAB-####).

The web portal (https://developmentcellatlas.ncl.ac.uk/####/) files were made using the single cell analysis bundle (v1.0.0, https://github.com/haniffalab/scRNA-seq_analysis), set up on the School of Computing Science HPC Facility, Newcastle University.





