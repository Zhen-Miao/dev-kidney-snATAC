Sys.setenv(RETICULATE_PYTHON = "~/.myvenv2/bin/python")
Sys.setenv(RETICULATE_PYTHON_ENV = "~/.myvenv2/bin/python")
library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(monocle3)
library(htmlwidgets)
library(OneR)
set.seed(123)
setwd("/yourdirectoryhere/P0")


#=========================================================================
#=========================================================================
#========= Random sampling of cells from interesting clusters ============
#=========================================================================
#=========================================================================

### STEP 1: Pick random sample of cells with equal nCells from a number of interesting clusters 
#https://github.com/satijalab/seurat/issues/243
#read in clustered information
P0_adult <- readRDS(file = "P0_adult_final.rds")

#subset P0 only.
Idents(P0_adult) <- P0_adult@meta.data$exp.cond
P0_adult <- subset(P0_adult, idents="P0")

#check nCells in clusters
Idents(P0_adult) <- P0_adult@meta.data$RNA_snn_res.0.4
#A) NP=cluster14
WhichCells(P0_adult, idents = "14")
length(WhichCells(P0_adult, idents = "14"))
#407 cells
#B) LOH=cluster9
WhichCells(P0_adult, idents = "9")
length(WhichCells(P0_adult, idents = "9"))
#1336 cells
#C) Podo=cluster11
WhichCells(P0_adult, idents = "11")
length(WhichCells(P0_adult, idents = "11"))
#928 cells
#D) PTS1=cluster3
WhichCells(P0_adult, idents = "3")
length(WhichCells(P0_adult, idents = "3"))
#684 cells
#E) PT_1
WhichCells(P0_adult, idents = "1")
length(WhichCells(P0_adult, idents = "1"))
#1513 cells
#G) Prolif1
WhichCells(P0_adult, idents = "7")
length(WhichCells(P0_adult, idents = "7"))
#1749 cells
#H) earlyPT
WhichCells(P0_adult, idents = "8")
length(WhichCells(P0_adult, idents = "8"))
#906 cells
#I) Stroma-like
WhichCells(P0_adult, idents = "2")
length(WhichCells(P0_adult, idents = "2"))
#3993 cells

# Pick 684 cells per cluster
cells.to.sample <- length(WhichCells(P0_adult, idents = "3"))

# Sample from other clusters as many cells as there are cells in cluster12
# For reproducibility, set a random seed
set.seed(cells.to.sample)
sampled.cells_LOH <- sample(x = WhichCells(P0_adult, idents = "9"), size = cells.to.sample, replace = F)
sampled.cells_Podo <- sample(x = WhichCells(P0_adult, idents = "11"), size = cells.to.sample, replace = F)
sampled.cells_PTS1 <- sample(x = WhichCells(P0_adult, idents = "3"), size = cells.to.sample, replace = F)
sampled.cells_PT_1 <- sample(x = WhichCells(P0_adult, idents = "1"), size = cells.to.sample, replace = F)
sampled.cells_Prolif1 <- sample(x = WhichCells(P0_adult, idents = "7"), size = cells.to.sample, replace = F)
sampled.cells_earlyPTDCT <- sample(x = WhichCells(P0_adult, idents = "8"), size = cells.to.sample, replace = F)
sampled.cells_Stromalike <- sample(x = WhichCells(P0_adult, idents = "2"), size = cells.to.sample, replace = F)

# Create a vector of cells on which variable genes will be computed
cells.to.subset <- c(WhichCells(P0_adult, idents = "14"), 
                     sampled.cells_LOH, 
                     sampled.cells_Podo, 
                     sampled.cells_PTS1, 
                     sampled.cells_PT_1,
                     sampled.cells_Prolif1,
                     sampled.cells_earlyPTDCT,
                     sampled.cells_Stromalike)

P0_adult.subset <- subset(P0_adult, cells = cells.to.subset)

### STEP 3: create CDS object
## Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(P0_adult.subset@assays[["RNA"]]@counts@Dimnames[[1]], 
                                 row.names = P0_adult.subset@assays[["RNA"]]@counts@Dimnames[[1]])
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information
cell_metadata <- as.data.frame(P0_adult.subset@assays[["RNA"]]@counts@Dimnames[[2]], 
                               row.names = P0_adult.subset@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix (from fresh sample)
New_matrix <- P0_adult.subset@assays[["RNA"]]@counts
expression_matrix <- New_matrix

## Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

#=========================================================================
#=========================================================================
#================= Clustering and classifying your cells =================
#=========================================================================
#=========================================================================
cds <- cds_from_seurat

#================= Step 1: Normalize & pre-process the data =================
# CPU intensive, time-consuming.
cds <- preprocess_cds(cds, num_dim = 100, cores=8)

#ElbowPlot
pdf(file = "/plots_1_Dimensions_ElbowPlot.pdf")
plot_pc_variance_explained(cds)
dev.off()

#================= Step 2: Reduce dimensionality and visualize the cells (UMAP by default) =================
# Note: reduce_dimension will produce slightly different output each time you run it unless you set 
# 'umap.fast_sgd = FALSE' and 'cores = 1
cds <- reduce_dimension(cds, umap.fast_sgd = FALSE, cores=1)

pdf(file = "/plots_2_Dimensions_UMAPall.pdf")
plot_cells(cds)
dev.off()

# Visualize how individual genes vary along the trajectory.
NP_genes <- c("Cited1", "Uncx", "Spock2", "Slc12a2")
pdf(file = "/plots_2_Dimensions_UMAP_NPgenes.pdf")
plot_cells(cds, genes=NP_genes,
           label_cell_groups=FALSE, norm_method="size_only", alpha=0.1, cell_size=0.5) + scale_colour_gradient(low = "grey85", high = "darkblue", na.value = "grey95")
dev.off()

#================= Step 3: Cluster the cells =================
cds = cluster_cells(cds, resolution=0.5e-3, k=29)
clusters(cds)

#by clusters
pdf(file = "/plots_3_Dimensions_UMAPclusters_res0.5e-3_k29.pdf")
plot_cells(cds)
dev.off()

#by partitions
pdf(file = "/plots_3_Dimensions_UMAPpartitions_res0.5e-3_k29.pdf")
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
dev.off()

#================= Step 4: Find marker genes expressed by each cluster =================
#by cluster
marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)

# The data frame marker_test_res contains a number of metrics for how specifically expressed 
# each gene is in each cluster We could group the cells according to cluster, partition, 
# or any categorical variable in colData(cds). You can rank the table according to one or more 
# of the specificity metrics and take the top gene for each cluster. For example, pseudo_R2 is
# one such measure. We can rank markers according to pseudo_R2 like this:
top50_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(50, pseudo_R2)

fwrite(x = top50_specific_markers, row.names = TRUE, file = "/top50_specific_markers.csv")

top1_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

fwrite(x = top1_specific_markers, row.names = TRUE, file = "/top1_specific_markers.csv")

top1_specific_marker_ids <- unique(top1_specific_markers %>% pull(gene_id))

# Now, we can plot the expression and fraction of cells that express each marker in each group 
# with the plot_genes_by_group function:
#by cluster
pdf(file = "/plots_4_DiffExpr_Dotplot_top1_bycluster.pdf")
plot_genes_by_group(cds,
                    top1_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)
dev.off()

#by partition
pdf(file = "/plots_4_DiffExpr_Dotplot_top1_bypartition.pdf")
plot_genes_by_group(cds,
                    top1_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)
dev.off()

# It's often informative to look at more than one marker, which you can do just by changing the 
# first argument to top_n():
top3_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top3_specific_marker_ids <- unique(top3_specific_markers %>% pull(gene_id))

#by cluster
pdf(file = "/plots_4_DiffExpr_Dotplot_top3_bycluster.pdf")
plot_genes_by_group(cds,
                    top3_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
dev.off()

#by partition
pdf(file = "/plots_4_DiffExpr_Dotplot_top3_bypartition.pdf")
plot_genes_by_group(cds,
                    top3_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3)
dev.off()

#================= Annotate your cells according to type =================
# To assign cell types based on clustering, we begin by creating a new column
# in colData(cds) and initialize it with the values of clusters(cds_from_seurat):
colData(cds)$assigned_cell_type <- as.character(clusters(cds))

# Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type:
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                "1"="PT",
                                                "2"="NP",
                                                "3"="Prolif",
                                                "4"="LOH",
                                                "5"="Podo", 
                                                "6"="Stroma-like")



#=========================================================================
#=========================================================================
#================= Constructing single-cell trajectories =================
#=========================================================================
#=========================================================================

#================= Learn the trajectory graph =================
#reduce partitions to 1
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions == "2"] <- "1"

# Next, we will fit a principal graph within each partition using the learn_graph() function:
# CPU, time consuming.
cds <- learn_graph(cds)
#cds <- learn_graph(cds, use_partition = FALSE)

pdf(file = "/plots_5.1_Trajectory_assignedcelltype.pdf")
plot_cells(cds,
           color_cells_by = "assigned_cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

#================= Order the cells in pseudotime =================
# Finding spots in the UMAP space that are occupied by cells from early time points:
pdf(file = "/plots_5.2_Trajectory_assignedcelltype.pdf",
    width=10, height=5)
plot_cells(cds,
           color_cells_by = "assigned_cell_type",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           graph_label_size=1.5)
dev.off()

# Define NP cluster as earliest timepoint & order cells according to pseudotime.
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, assigned_cell_type="NP"){
  cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

pdf(file = "/plots_5.3_Trajectory_pseudotime.pdf",
    width=10, height=5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           graph_label_size=1.5)
dev.off()

#store pseudotime in colData
colData(cds)$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime



#=========================================================================
#=========================================================================
#=================== Differential expression analysis ====================
#=========================================================================
#=========================================================================

#================= FIRST approach: fit_models function =================
### 1) Podo trajectory
#store pseudotime in colData of cds
colData(cds)$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime

# subset only clusters from Podo trajectory (clusters NP, Prolif, Stroma-like, Podo)
cds_subsetPodotraj_clusters <- c("NP", "Prolif", "Podo")
cds_subsetPodotraj <- cds[,colData(cds)$assigned_cell_type %in% cds_subsetPodotraj_clusters]

#Fit model to identify genes that vary with pseudotime.
#CPU intensive
all.genes <- cds_subsetPodotraj@rowRanges@partitioning@NAMES
all.genes_cds_Podotraj <- cds[rowData(cds_subsetPodotraj)$gene_short_name %in% all.genes,]
gene_fits_Podotraj <- fit_models(all.genes_cds_Podotraj, model_formula_str = "~pseudotime")

#Now let's see which of these genes have time-dependent expression. 
#First, we extract a table of coefficients from each model using the coefficient_table() function:
fit_coefs_Podotraj <- coefficient_table(gene_fits_Podotraj)
#fwrite(x = fit_coefs_Podotraj[,-c(2,3)] , file = "/fit_coefs_Podotraj.csv")
emb_time_terms_Podotraj <- fit_coefs_Podotraj %>% filter(term == "pseudotime")
fwrite(x = emb_time_terms_Podotraj[,-c(2,3)] , file = "/emb_time_terms_Podotraj.csv")

#Pull out the genes that have a significant pseudotime component
emb_time_terms_Podotraj_sign_q <- emb_time_terms_Podotraj %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

fwrite(x = emb_time_terms_Podotraj_sign_q, file = "/emb_time_terms_Podotraj_sign_q.csv")



### 2) PT trajectory
#store pseudotime in colData of cds
colData(cds)$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime

# subset only clusters from PT trajectory (clusters NP, Prolif, Stroma-like, PT)
cds_subsetPTtraj_clusters <- c("NP", "Prolif", "Stroma-like", "PT")
cds_subsetPTtraj <- cds[,colData(cds)$assigned_cell_type %in% cds_subsetPTtraj_clusters]
#unique(colData(cds_subsetPTtraj)$assigned_cell_type)

#Fit model to identify genes that vary with pseudotime.
#CPU intensive
all.genes <- cds_subsetPTtraj@rowRanges@partitioning@NAMES
all.genes_cds_PTtraj <- cds[rowData(cds_subsetPTtraj)$gene_short_name %in% all.genes,]
gene_fits_PTtraj <- fit_models(all.genes_cds_PTtraj, model_formula_str = "~pseudotime")

#Now let's see which of these genes have time-dependent expression. 
#First, we extract a table of coefficients from each model using the coefficient_table() function:
fit_coefs_PTtraj <- coefficient_table(gene_fits_PTtraj)
#fwrite(x = fit_coefs_PTtraj[,-c(2,3)] , file = "/fit_coefs_PTtraj.csv")
emb_time_terms_PTtraj <- fit_coefs_PTtraj %>% filter(term == "pseudotime")
fwrite(x = emb_time_terms_PTtraj[,-c(2,3)] , file = "/emb_time_terms_PTtraj.csv")

#Pull out the genes that have a significant pseudotime component
emb_time_terms_PTtraj_sign_q <- emb_time_terms_PTtraj %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

fwrite(x = emb_time_terms_PTtraj_sign_q, file = "/emb_time_terms_PTtraj_sign_q.csv")



### 3) LOH trajectory
#store pseudotime in colData of cds
colData(cds)$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime

# subset only clusters from LOH trajectory (clusters NP, Prolif, Stroma-like, LOH)
cds_subsetLOHtraj_clusters <- c("NP", "Prolif", "Stroma-like", "LOH")
cds_subsetLOHtraj <- cds[,colData(cds)$assigned_cell_type %in% cds_subsetLOHtraj_clusters]
#unique(colData(cds_subsetLOHtraj)$assigned_cell_type)

#Fit model to identify genes that vary with pseudotime.
#CPU intensive
all.genes <- cds_subsetLOHtraj@rowRanges@partitioning@NAMES
all.genes_cds_LOHtraj <- cds[rowData(cds_subsetLOHtraj)$gene_short_name %in% all.genes,]
gene_fits_LOHtraj <- fit_models(all.genes_cds_LOHtraj, model_formula_str = "~pseudotime")

#Now let's see which of these genes have time-dependent expression. 
#First, we extract a table of coefficients from each model using the coefficient_table() function:
fit_coefs_LOHtraj <- coefficient_table(gene_fits_LOHtraj)
#fwrite(x = fit_coefs_LOHtraj[,-c(2,3)] , file = "/fit_coefs_LOHtraj.csv")
emb_time_terms_LOHtraj <- fit_coefs_LOHtraj %>% filter(term == "pseudotime")
fwrite(x = emb_time_terms_LOHtraj[,-c(2,3)] , file = "/emb_time_terms_LOHtraj.csv")

#Pull out the genes that have a significant pseudotime component
emb_time_terms_LOHtraj_sign_q <- emb_time_terms_LOHtraj %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

fwrite(x = emb_time_terms_LOHtraj_sign_q, file = "/emb_time_terms_LOHtraj_sign_q.csv")




#================= SECOND approach: graph_test function =================
##### A) in all cells

# We turn to graph_test(), this time passing it neighbor_graph="principal_graph", which tells it to test 
# whether cells at similar positions on the trajectory have correlated expression:
#CPU intensive! ca. 2min@1888cells
#The data frame pr_graph_test_res has the Moran's I test results for each gene in the cell_data_set.
# effect size: morans_I (ranges from -1 to +1). A value of 0 indicates no effect, while +1 indicates
#perfect positive autocorrelation and suggests that a nearby cells have very similar values of a gene's expression. 
#Significant values much less than zero are generally rare.
pr_graph_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
fwrite(x = pr_graph_test_res, file = "/pr_graph_test_res.csv")
saveRDS(pr_graph_test_res, file = "/_OUT_pr_graph_test_res.rds")

# Show genes varying significantly along pseudotime:
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

#Here are a couple of interesting genes that score as highly significant according to graph_test():
#Aldob        PT
#Slc34a1      PT
#Spink1       PT
#Sult1d1      PT
#Podxl        Podo
#Nphs1        Podo
#Umod         LOH
#Slc12a1      LOH
#Sostdc1      LOH
#Cldn10       LOH
pdf(file = "/plots_6.2_DEG_sigvarmarkers.pdf",
    width=20, height=10)
plot_cells(cds, genes=c("Aldob", "Slc34a1", "Spink1", "Sult1d1", "Podxl", "Nphs1", "Umod", "Slc12a1", 
                        "Sostdc1", "Cldn10"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
dev.off()

# We can collect the trajectory-variable genes into modules:
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
fwrite(x = gene_module_df, row.names = TRUE, file = "/gene_module_df.csv")

# Here we plot the aggregate module scores within each group of cell types as annotated by assigned_cell_type:
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$assigned_cell_type)
fwrite(x = cell_group_df, row.names = TRUE, file = "/cell_group_df.csv")

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pdf(file = "/plots_6.3_DEG_modules.pdf", width=10, height=20)
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
dev.off()

##### B) in PT trajectory
cds_subsetPTtraj_clusters <- c("NP", "Prolif", "Stroma-like", "PT")
cds_subsetPTtraj <- cds[,colData(cds)$assigned_cell_type %in% cds_subsetPTtraj_clusters]

pr_graph_test_res_PTtraj <- graph_test(cds_subsetPTtraj, neighbor_graph="principal_graph", cores=8)
fwrite(x = pr_graph_test_res_PTtraj , file = "/pr_graph_test_res_PTtraj.csv")
saveRDS(pr_graph_test_res_PTtraj, file = "/_OUT_pr_graph_test_res_PTtraj.rds")

pr_deg_ids_PTtraj <- row.names(subset(pr_graph_test_res_PTtraj, q_value < 0.05))
write.table(x = pr_deg_ids_PTtraj , file = "/pr_deg_ids_PTtraj.csv")

gene_module_df_PTtraj <- find_gene_modules(cds_subsetPTtraj[pr_deg_ids_PTtraj,], resolution=1e-2)
fwrite(x = gene_module_df_PTtraj, row.names = TRUE, file = "/gene_module_df_PTtraj.csv")

# Here we plot the aggregate module scores within each group of cell types as annotated by assigned_cell_type:
cell_group_df_PTtraj <- tibble::tibble(cell=row.names(colData(cds_subsetPTtraj)), 
                                       cell_group=colData(cds_subsetPTtraj)$assigned_cell_type)
fwrite(x = cell_group_df_PTtraj, row.names = TRUE, file = "/cell_group_df_PTtraj.csv")

agg_mat_PTtraj <- aggregate_gene_expression(cds_subsetPTtraj, gene_module_df_PTtraj, cell_group_df_PTtraj)
row.names(agg_mat_PTtraj) <- stringr::str_c("Module ", row.names(agg_mat_PTtraj))

pdf(file = "/plots_6.3_DEG_modules_PTtraj.pdf", width=10, height=20)
pheatmap::pheatmap(agg_mat_PTtraj,
                   scale="column", clustering_method="ward.D2")
dev.off()


##### C) in LOH trajectory
cds_subsetLOHtraj_clusters <- c("NP", "Prolif", "Stroma-like", "LOH")
cds_subsetLOHtraj <- cds[,colData(cds)$assigned_cell_type %in% cds_subsetLOHtraj_clusters]

pr_graph_test_res_LOHtraj <- graph_test(cds_subsetLOHtraj, neighbor_graph="principal_graph", cores=8)
fwrite(x = pr_graph_test_res_LOHtraj , file = "/pr_graph_test_res_LOHtraj.csv")
saveRDS(pr_graph_test_res_LOHtraj, file = "/_OUT_pr_graph_test_res_LOHtraj.rds")

pr_deg_ids_LOHtraj <- row.names(subset(pr_graph_test_res_LOHtraj, q_value < 0.05))

gene_module_df_LOHtraj <- find_gene_modules(cds_subsetLOHtraj[pr_deg_ids_LOHtraj,], resolution=1e-2)
fwrite(x = gene_module_df_LOHtraj, row.names = TRUE, file = "/gene_module_df_LOHtraj.csv")

cell_group_df_LOHtraj <- tibble::tibble(cell=row.names(colData(cds_subsetLOHtraj)), 
                                         cell_group=colData(cds_subsetLOHtraj)$assigned_cell_type)
fwrite(x = cell_group_df_LOHtraj, row.names = TRUE, file = "/cell_group_df_LOHtraj.csv")

agg_mat_LOHtraj <- aggregate_gene_expression(cds_subsetLOHtraj, gene_module_df_LOHtraj, cell_group_df_LOHtraj)
row.names(agg_mat_LOHtraj) <- stringr::str_c("Module ", row.names(agg_mat_LOHtraj))

pdf(file = "/plots_6.3_DEG_modules_LOHtraj.pdf", width=10, height=20)
pheatmap::pheatmap(agg_mat_LOHtraj,
                   scale="column", clustering_method="ward.D2")
dev.off()


##### D) in Podo trajectory
cds_subsetPodotraj_clusters <- c("NP", "Prolif", "Podo")
cds_subsetPodotraj <- cds[,colData(cds)$assigned_cell_type %in% cds_subsetPodotraj_clusters]

pr_graph_test_res_Podotraj <- graph_test(cds_subsetPodotraj, neighbor_graph="principal_graph", cores=8)
fwrite(x = pr_graph_test_res_Podotraj , file = "/pr_graph_test_res_Podotraj.csv")
saveRDS(pr_graph_test_res_Podotraj, file = "/_OUT_pr_graph_test_res_Podotraj.rds")

pr_deg_ids_Podotraj <- row.names(subset(pr_graph_test_res_Podotraj, q_value < 0.05))

gene_module_df_Podotraj <- find_gene_modules(cds_subsetPodotraj[pr_deg_ids_Podotraj,], resolution=1e-2)
fwrite(x = gene_module_df_Podotraj, row.names = TRUE, file = "/gene_module_df_Podotraj.csv")

cell_group_df_Podotraj <- tibble::tibble(cell=row.names(colData(cds_subsetPodotraj)), 
                                         cell_group=colData(cds_subsetPodotraj)$assigned_cell_type)
fwrite(x = cell_group_df_Podotraj, row.names = TRUE, file = "/cell_group_df_Podotraj.csv")

agg_mat_Podotraj <- aggregate_gene_expression(cds_subsetPodotraj, gene_module_df_Podotraj, cell_group_df_Podotraj)
row.names(agg_mat_Podotraj) <- stringr::str_c("Module ", row.names(agg_mat_Podotraj))

pdf(file = "/plots_6.3_DEG_modules_Podotraj.pdf", width=10, height=20)
pheatmap::pheatmap(agg_mat_Podotraj,
                   scale="column", clustering_method="ward.D2")
dev.off()


#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)



