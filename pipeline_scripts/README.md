# This directory contains generalisable scripts re-used throughout the FBM manuscript for analysis of single cell data.

These scripts include: (1) Scrublet doublet removal, (2) Pre-processing including: transforming count matrix, batch correction, adding dimensional reduction, clustering (3) Gene-set over-representation analysis and (4) Negative binomial barplot statistics and (5) logistic regression.


## pipeline_scrublet.ipynb
Script to --

## pipeline_add_dr_harmony_clus_degs_annot.ipynb
script to --

## pipeline_Negative_binomial_Quasi_binomial_test_barplot_stats.R
This script was used to compute statistical significance of proportion flux between cellstates of interest in different tissues and developmental stages in all barplots. 
The script takes metadata from sc objects as an input and requires it's user to define:
- The preferred model to model the data, options are either "nb" or "quasi" 
- the categorical variable to test between (e.g developmental_stage/age_groups)
- the categorical variable which is in flux (e.g cell.labels, lineage)
- the categorical variable for batch (e.g orig.ident/lane)

The output is formated as two .csv files which contain the statistical output for the liklihood ratio test for signficance of proportion flux and spearman's rho test for characterising the trend direction. 

## pipeline_Gene_set_over_representation_analysis.R
The script was used to produce the chord plots in Fig4E and Fig4F.


## pipeline_logistic_regression_alginment.ipynb

## pipeline_logit_regression_PCA_model_stability_validation

## Interactive heatmap dotplot

## pipeline_pseudotime_webportal

## pipeline_web_portal
