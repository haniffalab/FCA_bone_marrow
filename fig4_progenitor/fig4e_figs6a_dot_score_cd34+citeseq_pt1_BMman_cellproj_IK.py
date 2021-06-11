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

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import tables as tb
import scipy as scipy
import dotscore

sc.set_figure_params(color_map = 'viridis', dpi_save = 350)
sc.settings.verbosity = 3

import cellproj

# %%
#To ensure 1:1 aspect ratio
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4, 4)

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# ## DoTscore using abm data from Simone

# %%
abm = sc.read('./data/abm_for_issac_20200717.h5ad')

# %%
sc.pp.normalize_total(abm, target_sum = 10000) # normalize with total UMI count per cell
sc.pp.log1p(abm)

# %%
coords = pd.read_csv('./data/figs2d_abm_umap_20210213.csv', index_col = 0)
coords = coords.loc[abm.obs.index,]
abm.obsm['X_umap'] = np.array(coords)

# %%
anno = pd.read_csv('./data/S Table 8_ abm_annot_20200717.csv', index_col = 0)
anno = anno.loc[abm.obs.index,]
abm.obs['anno'] = anno

anno_broad_key = pd.read_csv('./data/Stab5_anno_key.csv', index_col = 0)
abm.obs['anno_broad'] = anno_broad_key.loc[abm.obs.anno,'broad_cell.labels'].values

# %%
sc.pp.highly_variable_genes(abm, flavor='seurat', n_top_genes = 5000)
sc.pl.highly_variable_genes(abm)

# %%
fetdata = sc.read('./data/080221_mq224_part3_RAW_mRNA_Progenitors_only_ccincluded_FBMS2-3-H3-F3-E5_CB-G7_FL-A7-A6-C5.h5ad')

# %%
fethsc = fetdata[fetdata.obs.index[fetdata.obs['Cell.labels.P4.sorted'] == 'HSC/MPP I'],:].copy()

# %%
target, ref = cellproj.process_together(fethsc, abm, use_vargenes = True, lognorm = False)

# %%
pops = ['CB', 'FBM', 'FL']
for i in pops:
    print(i)
    abm.obs[i + '_proj'] = cellproj.project_to_ref(target[target.obs.Tissue == i,:], ref, k = 15)[1]

# %%
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#DDD3D3'),
                                                        (0.001, 'tan'),
                                                        (1, 'blue')])
sc.pl.umap(abm, color = ['CB_proj', 'FBM_proj', 'FL_proj'], cmap = cmap2, save = '_fethsc_proj_onto_abm.pdf', frameon = False)

# %%
for i in ['CB', 'FBM', 'FL']:
    print(sum(abm.obs[i + '_proj'] >0.05))

# %%
mpl.rcParams['figure.figsize'] = (36, 4)
sc.pl.violin(abm, keys = ['CB_proj'], groupby='anno', save = 'CB_proj.png')
sc.pl.violin(abm, keys = ['FL_proj'], groupby='anno', save = 'FL_proj.png')
sc.pl.violin(abm, keys = ['FBM_proj'], groupby='anno', save = 'FBM_proj.png')
mpl.rcParams['figure.figsize'] = (4, 4)

# %%
abm.write('./procdata/abm_annot.h5ad', compression = 'lzf')

# %%
sc.pl.umap(abm, color = ['anno'], save = '_abm_annotation.pdf')
sc.pl.umap(abm, color = ['anno_broad'], save = '_abm_annotation_broad.pdf')

