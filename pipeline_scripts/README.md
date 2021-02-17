# This directory contains generalisable scripts re-used throughout the FBM manuscript for analysis of single cell data.

## pipeline_scrublet.ipynb
For all 10x and CITE-seq (mRNA lane) data shown in manuscript figures.
This script was used to perform doublet detection on all 10x data (including DS/non-DS fetal BM, ABM, CB) and CITE-seq (including total fetal BM and CD34+ fetal BM, FL, CB) using the scrublet package which estimates a scrublet score for each cell in the input data.

## pipeline_add_dr_harmony_clus_degs_annot.ipynb
For all 10x data shown in figures.
This script was used to perform batch correction on all data and computes harmony corrected PC coordinates for an input scRNA-seq dataset.

## pipeline_Negative_binomial_Quasi_binomial_test_barplot_stats.R
For Fig 1D, 2B, 3C, 5B, E5E, E8A.
This script was used to compute statistical significance of proportion flux between cellstates of interest in different tissues and developmental stages in all barplots. Written by IG.
The script takes metadata from sc objects as an input and requires it's user to define:
- The preferred model to model the data, options are either "nb" or "quasi" 
- the categorical variable to test between (e.g developmental_stage/age_groups)
- the categorical variable which is in flux (e.g cell.labels, lineage)
- the categorical variable for batch (e.g orig.ident/lane)

The output is formated as two .csv files which contain the statistical output for the liklihood ratio test for signficance of proportion flux and spearman's rho test for characterising the trend direction. 

## pipeline_logist_general.ipynb
The notebook describes a method for integrating data by label transfering based on ridge regulariased logistic regression. This approach enables us to fit a model on the annotated training/landscape/reference dataset to predict labels of a new dataset. Logistic regression can (1) be used to classify samples, (2) use different types of data (continuous and descrete measurements)and (3) also be used to assess what variables are useful for classifying samples. 
This script was used to transfer labels to all single cell datastes used in this study, with fetal BM 10x data as reference. 

## pipeline_logistic_regression_alignment.ipynb
For Fig E1D, E5D, E7B
This script was used to compare analogous annotations between datasets to ensure accuracy and consistency of annotation and produces a heatmap weighted by probability of alignment derived from the binary assignment of all cells to a category by logistic regression. Written by IG.
The script takes as input:
- the analogous categorical variable in each dataset to be compared (e.g cell.labels)

## pipeline_logit_regression_PCA_model_stability_validation
Results not currently shown in current figures/described in manuscript, but this script was used to optimise selection of clustering parameters used in datasets. 
This script was used to compare the stability of cluster partitions and thus annotations across different PCA inpits across 5-70 PCs. The script assigns label predictions to newly computed clusters based on differring PC values by logistic regression and computes a stability score by RAND index, mutual information score and silhouette score comparing the original annotations to automatically assigned annotations. Written by IG.
The script takes as input:
- The range of PCs to test
- The categorical variable to project (e.g cell.labels)

## Interactive heatmap dotplot
This script was used to compute the interactive heatmap portals located (https://developmentcellatlas.ncl.ac.uk/fbm_index). Written by IG.
The script takes an option file as input containing:
- a vector of categorical variables to be used (e.g cell.labels)

An example of the options file may be found at (https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data)

## pipeline_web_portal
This script was used to compute the gene expression portals located (https://developmentcellatlas.ncl.ac.uk/fbm_index). Written by IG.
The script takes an option file as input containing:
- a vector of categorical variables to be used (e.g cell.labels)
- a vector of the dimension reduction slots in the data (e.g X-umamp, X-PCA)

An example of the options file may be found at (https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data)
