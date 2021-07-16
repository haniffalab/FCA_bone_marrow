import pandas as pd
import scanpy as sc
import numpy as np

#Takes in an array and returns a boolean array with True set to k smallest values
def kfilter(x, k = 5):
    tr = np.sort(x)[k-1]
    #print(tr)
    return(x <= tr)

#Function which runs kfilter on a 2 dimensional array
def kneighbour_mask(array, k, axis = -1):
    array01 = np.apply_along_axis(kfilter, axis= axis, arr=array, k = k)
    return(array01)

# Function for processing together two adata objects (normalise, log2 and scale)
def process_together(target_adata, ref_adata, pca_n_components = 50, use_vargenes = True, lognorm = True):

    from scipy.sparse import issparse

    #Subsetting for common genes
    combo_adata = ref_adata.concatenate(target_adata)
    target_adata = target_adata[:, combo_adata.var.index]
    ref_adata = ref_adata[:, combo_adata.var.index]

    target_adata.obs['source'] = 'target'
    ref_adata.obs['source'] = 'ref'

    combo_adata = ref_adata.concatenate(target_adata)
    if lognorm:
        sc.pp.normalize_per_cell(combo_adata)
        sc.pp.log1p(combo_adata)
    sc.pp.scale(combo_adata)

        # print(combo_adata)
    if use_vargenes:
        combo_adata = combo_adata[:,combo_adata.var['highly_variable-0']]

    target_adata = combo_adata[combo_adata.obs.source == "target",:].copy()
    ref_adata = combo_adata[combo_adata.obs.source == "ref",:].copy()

    if issparse(target_adata.X):
         target_adata.X = target_adata.X.todense()
    if issparse(ref_adata.X):
         ref_adata.X = ref_adata.X.todense()

    return(target_adata, ref_adata)

#Function to project cells onto a reference dataset
# Returns a tuple, friends is a matrix indicating nearest neighbors and friendusms the number of nearest neighbors per cell in the reference dataset
def project_to_ref(target_adata, ref_adata, pca_n_components = 50, k = 10):

    from sklearn.decomposition import PCA
    from sklearn.metrics import pairwise_distances

    pca1 = PCA(n_components=pca_n_components, svd_solver='auto', random_state=0)
    pca1.fit(ref_adata.X)

    ref_pca = pca1.transform(ref_adata.X)
    target_pca = pca1.transform(target_adata.X)

    combo_pca = np.concatenate((ref_pca, target_pca), 0)

    dists = pairwise_distances(target_pca, ref_pca)

    friends = kneighbour_mask(dists, k = k)
    friendsums = np.sum(friends, axis = 0) / target_adata.shape[0] #The sums are normalised to the number of target cell

    return(friends, friendsums)

