library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)

#================= Load datasets and create Seurat objects =================
Ksp.data <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
Ksp <- CreateSeuratObject(counts = Ksp.data, project = "Ksp")
Ksp

P0_batch2.data <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
P0_batch2 <- CreateSeuratObject(counts = P0_batch2.data, project = "P0_batch2")
P0_batch2

#================= Merge =================
#To merge more than two Seurat objects, simply pass a vector of multiple Seurat objects to the y parameter for merge.
P0_adult <- merge(Ksp, y = c(P0_batch2), 
                  add.cell.ids = c("Ksp", "P0_batch2"),
                  project = "P0_adult")
P0_adult

#================= Add experimental setup to metadata =================
#add exp.cond
P0_adult@meta.data$exp.cond <- plyr::mapvalues(
  x = P0_adult@meta.data$orig.ident,
  from = c("Ksp", "P0_batch2"),
  to = c("Adult", "P0")
)

unique(sapply(X = strsplit(colnames(P0_adult), split = "_"), FUN = "[", 1))
cells_prefilter <- table(P0_adult$orig.ident)
write.table(x = cells_prefilter, file = "/cells_prefilter.csv", 
            append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = NA)

#================= Pre-processing =================
P0_adult[["percent.mt"]] <- PercentageFeatureSet(P0_adult, pattern = "^mt-")

#subsetting high quality cells
P0_adult <- subset(P0_adult, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 50)
cells_postfilter <- table(P0_adult$orig.ident)
write.table(x = cells_postfilter, file = "/cells_postfilter.csv", 
            append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = NA)

proportion_postfilter <- cells_postfilter/cells_prefilter
write.table(x = proportion_postfilter, file = "/proportion_postfilter.csv", 
            append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = NA)

P0_adult <- NormalizeData(P0_adult, normalization.method = "LogNormalize", scale.factor = 10000)
P0_adult <- FindVariableFeatures(P0_adult, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(P0_adult)
P0_adult <- ScaleData(P0_adult, features = all.genes)
P0_adult <- RunPCA(P0_adult, features = VariableFeatures(object = P0_adult))

#============= Run Harmony =============
library(harmony)

#Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
pdf(file = "/plots_7_Harmony.pdf", height=2.5, width=6)
P0_adult <- P0_adult %>% RunHarmony("orig.ident", plot_convergence = TRUE)
dev.off()

harmony_embeddings <- Embeddings(P0_adult, 'harmony')
harmony_embeddings[1:5, 1:5]

#============= Downstream analysis =============
P0_adult <- P0_adult %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), save.SNN=TRUE) %>% 
  identity()

#proceed with analysis with cluster resolution 0.4
Idents(P0_adult) <- P0_adult@meta.data$RNA_snn_res.0.4
P0_adult$seurat_clusters <- P0_adult$RNA_snn_res.0.4

#================= Average feature expression across cell cluster identity, invert log and write out =================
cluster.averages <- AverageExpression(P0_adult, return.seurat = FALSE, 
                                      slot = "scale.data", 
                                      verbose = TRUE)
write.table(x = cluster.averages, file = "/P0_adult_cluster.average_outfile.txt", 
            append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = NA)

#================= Finding differentially expressed features (cluster biomarkers) =================
P0_adult.markers <- FindAllMarkers(P0_adult, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



