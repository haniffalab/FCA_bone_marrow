# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# container rpy_v3.1

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import tables as tb
import scipy as scipy
import dotscore

#Specifying random seed
import random

sc.set_figure_params(color_map = 'viridis', dpi_save = 350)
sc.settings.verbosity = 3

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = (4, 4)
mpl.rcParams['pdf.fonttype'] = 42 #Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42

# %%
# Setting up target directories
sc.settings.figdir = './figures/02script/'
base_figures = './figures/02script/'
base_procdata = './procdata/02script/'
import pathlib
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## DoTscore using abm data from Simone

# %%
abm = sc.read('./procdata/abm_annot.h5ad')

# %% [markdown]
# ### FL vs FBM pairwise

# %%
abmFL = abm.copy()
meansFL = abm[abm.obs.FL_proj > 0.05,:].X.mean(axis = 0).copy()
dotscore.custom_scale(abmFL, mean = meansFL)

de1 = pd.read_csv("./data/080221_mq224_DEGWilcoxon_BHcorrection_FL-FBM_groups_mRNA_HSCMPP-I_only_ccincluded_P4_sorted_annotations_FBMS2-3-H3-F3-E5_CB-G7_FL-A7-A6-C5.csv", index_col=0)
de1sig = de1.loc[de1['HSC/MPP I_FL_pvals_a'] < 0.05,:]
de1sig = de1sig.loc[abs(de1sig['HSC/MPP I_FL_logfold']) > 0.5,:]
print(de1sig.shape)

allfolds = de1sig['HSC/MPP I_FL_logfold'].values
allgenes = de1['HSC/MPP I_FL_names']

abmFL.obs["FL_v_FBM_dotscore"] = dotscore.get_DoTscore(abmFL, 
                                                      de = de1sig, 
                                                      allgenes = allgenes, 
                                                      allfolds = allfolds,
                                                      zscore = True, 
                                                      id_col = 'HSC/MPP I_FL_names',
                                                   weight_col='HSC/MPP I_FL_logfold')
abmFL.obs["FL_v_FBM_dotscore"] = dotscore.qfilt(abmFL.obs["FL_v_FBM_dotscore"], 0.999, 0.001)

c1 = dotscore.cmap_RdBu(abmFL.obs["FL_v_FBM_dotscore"])
sc.pl.umap(abmFL, color=["FL_v_FBM_dotscore"], alpha = 0.7, color_map = c1, frameon = False, save = 'ambFL_FL_v_FBM_DoTscore.pdf')

scores = dotscore.get_genescore_pergroup(abmFL, de1sig, group = 'anno', sortby = 'early erythroid', id_col = 'HSC/MPP I_FL_names', weight_col='HSC/MPP I_FL_logfold')
scores = scores.loc[scores.sum(axis = 1) != 0, :]
scores.to_csv(base_procdata + 'ambFL_FL_v_FBM_DoTscore_per_cluster.csv')

# %% [markdown]
# ### Each tissue vs the others

# %%
de2 = pd.read_csv('./data/080221_mq224_DEGWilcoxon_BHcorrection_FL-FBM-CBgroups_mRNA_HSCMPP-I_only_ccincluded_P4_sorted_annotations_FBMS2-3-H3-F3-E5_CB-G7_FL-A7-A6-C5.csv', index_col=0)
for tissue, i in zip(['FL', 'FBM', 'CB'], ['HSCMPPI_FL_', 'HSCMPPI_FBM_', 'HSCMPPI_CB_']):
    abm_temp = abm.copy()
    means_temp = abm[abm.obs[tissue + '_proj'] > 0.05,:].X.mean(axis = 0).copy()
    dotscore.custom_scale(abm_temp, mean = means_temp)
    
    de = de2.filter(regex = i)
    
    desig = de.loc[de[i + 'pvals_a'] < 0.05,:]
    desig = desig.loc[abs(desig[i + 'logfold']) > 1,:]
    print(desig.shape)

    allfolds = desig[i + 'logfold'].values
    allgenes = de[i + 'names']
    
    abm_temp.obs[i + "_dotscore"] = dotscore.get_DoTscore(abm_temp, 
                                                          de = desig, 
                                                          allgenes = allgenes, 
                                                          allfolds = allfolds,
                                                          zscore = True, 
                                                          id_col = i + 'names',
                                                       weight_col= i + 'logfold')
    abm_temp.obs[i + "_dotscore"] = dotscore.qfilt(abm_temp.obs[i + "_dotscore"], 0.999, 0.001)

    if min(abm_temp.obs[i + "_dotscore"])<0 and max(abm_temp.obs[i + "_dotscore"])>0:
        c1 = dotscore.cmap_RdBu(abm_temp.obs[i + "_dotscore"])
    elif min(abm_temp.obs[i + "_dotscore"])>0:
        from matplotlib.colors import LinearSegmentedColormap
        c1 = LinearSegmentedColormap.from_list('mycmap', [(0, 'white'),
                                                       (1, 'red')])
    else:
        c1 = 'viridis'
    sc.pl.umap(abm_temp, color=[i + "_dotscore"], alpha = 0.7, color_map = c1, frameon = False, save = 'abm' + tissue + '_' + i + '_DoTscore.pdf')

    
    scores = dotscore.get_genescore_pergroup(abm_temp, desig,
                                                    group = 'anno',
                                                    sortby = 'early erythroid',
                                                    id_col = i + 'names',
                                                    weight_col= i + 'logfold')
    scores = scores.loc[scores.sum(axis = 1) != 0, :]
    scores.to_csv(base_procdata + 'abm' + tissue + '_' + i + '_DoTscore_per_cluster.csv')

