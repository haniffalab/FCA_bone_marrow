# Import libraries
library(scran)
library(RColorBrewer)
library(slingshot)
library(monocle)
library(gam)
library(clusterExperiment)
library(ggplot2)
library(plyr)
library(magrittr) 
library(dplyr)

# Download the b cell data as a CDS object from the monocle 2 obj (is raw and has HVGs imported from scanpy)
ie_regions_cds <- readRDS("/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/data/fig5e_5f_fbm_eryth_cds_20200625.RDS")
ie_regions_cds

# set variables
figs.addr = "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/figs/"
fig_prefix="fbm_eryth_20200625"
progenitor = 'HSC'
celltype.trajectory = c('HSC', 'MEP_MEMP', 'early erythroid')
palette.colours = c('#d200d2', '#00a500', '#cbcbff')

# Load data into a monocle3 object
fd <- fData(ie_regions_cds)
fd["gene_short_name"] = rownames(fd)
pd <- pData(ie_regions_cds)
exp <- exprs(ie_regions_cds)
#devtools::install_github('cole-trapnell-lab/leidenbase')
#usethis::browse_github_pat()
#usethis::edit_r_environ()
#devtools::install_url('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz')
#devtools::install_github('cole-trapnell-lab/monocle3')
# pre-process data 
library(monocle3)
cds <- new_cell_data_set(exp, gene_metadata = fd, cell_metadata = pd)
set.seed(26)
cds <- preprocess_cds(cds, num_dim = 50, method="PCA", verbose=TRUE)

# Import FDG coordinates from scanpy, or run FDG/batch correction in Monocle
fdg_coordinates <- read.csv("/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/fbm_eryth_fdg_20200625.csv", row.names = 1)
cds@int_colData@listData[["reducedDims"]][["UMAP"]] <- fdg_coordinates
plot_cells(cds, color_cells_by="orig.ident")
plot_cells(cds, color_cells_by="cell.labels")
plot_cells(cds, color_cells_by="sequencing.type")

# Fit a principle graph within each partition/cluster using the learn_graph() function
cds <- cluster_cells(cds, reduction_method = "UMAP", verbose=TRUE, resolution = 1e-07)
# Maximal modularity is 0 for resolution parameter 1e-07

# plot pst path by celltype (before ordering and assigning pst vals for each cell)
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "cell.labels", cell_size = 2, trajectory_graph_color="black", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE) 
p + scale_colour_manual(values=palette.colours) 
ggsave(file.path(figs.addr, "fbm_eryth_monocle3_pst_by_celltype_20200625.pdf"), width=12, height=10)
p + scale_colour_manual(values=palette.colours) 

# Order cells in pseudotime (to colour by pst val)- here, an RShiny app will pop up and I will manually choose a start point
#cds <- order_cells(cds)
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, cell.labels="HSC"){
  cell_ids <- which(colData(cds)[, "cell.labels"] == cell.labels)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

# Plot the final 2D trajectory
p = plot_cells(cds, color_cells_by = "pseudotime", cell_size = 3, trajectory_graph_color="black", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
p
ggsave(file.path(figs.addr, "fbm_eryth_monocle3_pst_by_pst_20200625_cellsize3.pdf"), width=12, height=10)
p

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Plot a 3d trajectory
#cds_3d <- reduce_dimension(cds, max_components = 3)
#cds_3d <- cluster_cells(cds_3d)
#cds_3d <- learn_graph(cds_3d)
#cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
#cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")

# ----------------------------------------------------------------------------------------------------------------------------------------------
# DEG analysis

# OVERALL DEGS
# run graph test with neighbour_graph="principal_graph" tells it to test whether cells nearby on traj have 
# correlated expression
# the graph test runs a Moran’s I test, a measure of multi-directional/dimensional spatial autocorrelation. 
# The statistic tells you whether cells at nearby positions on a trajectory will have similar (or dissimilar) expression levels for 
# the gene being tested. +1 means that nearby cells will have perfectly similar expression; 0 represents no correlation, and -1 
# means that neighboring cells will be anti-correlated.
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
cds_pr_test_res
pr_deg_ids <- subset(cds_pr_test_res, q_value < 0.05)
pr_deg_ids
write.csv(pr_deg_ids, "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/fbm_eryth_monocle3_overall_degs_20200625.csv")

# OVERALL DEG MODULE PLOTS
# collect trajectory-variable genes into modules
# gene_module_df is a df showing genes and the module they belong to
# Cluster genes into modules that are co-expressed across cells.
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100), verbose=TRUE, random_seed = 26)
write.csv(gene_module_df, "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/fbm_eryth_monocle3_overall_degs_gene_modules_20200625.csv")

# Maximal modularity is 0.89035665176761 for resolution parameter 0.01
#Run kNN based graph clustering DONE.
#  -Number of clusters: 38

# combine the DEG info df (pr_deg_ids) with the infor on the clusters (gene_module_df)
pr_deg_ids <- subset(cds_pr_test_res, q_value < 0.05)
pr_deg_ids$id = rownames(pr_deg_ids)
gene_module_df_final <- as.data.frame(gene_module_df)
merged <- merge(pr_deg_ids, gene_module_df_final, by="id")
head(merged)
write.csv(merged, "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/fbm_eryth_monocle3_overall_merged_deg_info_20200625.csv")

# plot the aggregate GEX of DEGs across cell types (which cells dominate which modules?)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$cell.labels)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
agg_mat
as.matrix(agg_mat)
write.csv(as.matrix(agg_mat), "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/fbm_eryth_monocle3_overall_agg_mat_20200625.csv")

# plot the module expression
#plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(2, 12, 14, 17, 4, 16, 21, 11, 13, 5, 18, 23, 19)), label_cell_groups=FALSE, show_trajectory_graph=FALSE)
p = plot_cells(cds, genes=gene_module_df, label_cell_groups=FALSE, show_trajectory_graph=FALSE)
p
ggsave(file.path(figs.addr, "fbm_eryth_monocle3_gene_modules_20200625.pdf"), width=10, height=10)
p
dev.off()

# plot a heatmap for the module scores by cell type
agg_mat
pheatmap::pheatmap(agg_mat, scale="row", clustering_method="ward.D2", cluster_rows = TRUE, cluster_cols = FALSE)
# now, save as "monocle3_heatmap_fbm_eryth_20200625.pdf" 15x10

final_order = c("Module 7", "Module 2", 
                "Module 6", "Module 12",  "Module 5", "Module 10",  "Module 3", 
                "Module 8", "Module 16",  
                 "Module 1", "Module 4", "Module 15",  "Module 11",  
                "Module 9", "Module 17", "Module 13", "Module 14")

celltype_order=c('HSC', 'MEP_MEMP', 'early erythroid')
                                                                                               
final_order
final_agg_mat <- agg_mat[final_order,]
final_agg_mat <- final_agg_mat[,celltype_order]
final_agg_mat
pheatmap::pheatmap(final_agg_mat, scale="row", clustering_method="ward.D2", cluster_rows = FALSE, cluster_cols = FALSE)
write.csv(as.matrix(final_agg_mat), "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/fbm_eryth_monocle3_overall_final_agg_mat_20200625.csv")
# now, save as "monocle3_heatmap_fbm_eryth_20200625.pdf" 15x10

# get pst values
p = plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1.5, trajectory_graph_color="black", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
p
point.data <- ggplot_build(p)[["plot"]][["data"]]
colnames(point.data)
point.data["cell_color"]
write.csv(point.data["cell_color"], "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/data/fbm_eryth_monocle3_pst_metadata_20200625.csv")
# save cds
saveRDS(cds, "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/data/fbm_eryth_cds_final_20200625.RDS")

# ----------------------------------------------------------------------------------------------------------------------------------------------

# RELOAD THE HEATMAP REAL QUICK AND REPLOT.....
ie_regions_cds <- readRDS("/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/data/fbm_eryth_cds_final_20200625.RDS")
ie_regions_cds
agg_mat=read.csv("/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/fbm_eryth_monocle3_overall_agg_mat_20200625.csv", row.names="X")
agg_mat
celltype_order=c('HSC', 'MEP_MEMP', 'early.erythroid')
final_order = c("Module 7", "Module 2", 
                "Module 6", "Module 12",  "Module 5", "Module 10",  "Module 3", 
                "Module 8", "Module 16",  
                 "Module 1", "Module 4", "Module 15",  "Module 11",  
                "Module 9", "Module 17", "Module 13", "Module 14")
final_order
final_agg_mat <- agg_mat[final_order,]
final_agg_mat <- final_agg_mat[,celltype_order]
final_agg_mat
pheatmap::pheatmap(final_agg_mat, scale="row", clustering_method="ward.D2", cluster_rows = FALSE, cluster_cols = FALSE)

# BRANCH ANALYSIS FOR NEGATIVELY SELECTED PRE B (select all of pre b quadrant - pst high=neg pre b here)
#cds_subset <- choose_cells(cds)
# run degs
#subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
#subset_pr_test_res
#pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
#pr_deg_ids
#gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100), verbose=TRUE, random_seed = 26)
#write.csv(gene_module_df, "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs/resources_for_pipelines/gene_module_df_neg_pre_b_20200327.csv")
# plot the module expression
#plot_cells(cds_subset, genes=gene_module_df, label_cell_groups=FALSE, show_trajectory_graph=FALSE)
#plot_cells(cds_subset, genes=gene_module_df %>% filter(module %in% c(4, 2, 3, 1, 5)), label_cell_groups=FALSE, show_trajectory_graph=FALSE)
# plot a heatmap for the modules
#cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$cell.labels)
#agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
#row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")

