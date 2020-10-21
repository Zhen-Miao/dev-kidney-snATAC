library(SCENIC)
library(AUCell)
library(RcisTarget)
library(doMC)
library(R2HTML)
library(rbokeh)
set.seed(123)
options(width=200)
setwd("/yourdirectoryhere/P0")

### ================================================================================
### ========================== Load data / Setup ===================================
### ================================================================================

#====================== Load Seurat object & extract expression matrix ===================
library(Seurat)
library(data.table)
P0_adult <- readRDS(file = "P0_adult_final.rds")
#subset P0
Idents(P0_adult) <- P0_adult$exp.cond
P0_adult_subsetP0 <- subset(P0_adult, idents = "P0")
#extract data matrix
exprMat <- as.matrix(GetAssayData(P0_adult_subsetP0, slot = "data"))

#====================== extract cell annotation info ===================
Idents(P0_adult_subsetP0) <- P0_adult_subsetP0$seurat_clusters
cellInfo <- data.frame(seuratCluster=Idents(P0_adult_subsetP0))

### Initialize SCENIC settings
#In order to keep consistent settings across the multiple steps of SCENIC, most functions in SCENIC package use a common object where the options for the current run are stored. 
#This object replaces the "arguments" for most functions, and should be created at the begining of a SCENIC run with the function `initializeScenic()`. 
#The default settings should be valid for most analyses. The parameters that need to be specified in all runs is the organism (`mgi` for mouse, `hgnc` for human, or `dmel` for fly), 
#and the directory where the RcisTarget databases are stored (you may create a link in the current directory to avoid duplicating them, e.g. in linux:
#` system("ln -s ~/path/to/dbs databases")`).
#For details on the options that can be modified check the help of `?initializeScenic` or of the specific function that takes it as input.
library(SCENIC)
org <- "mgi" # or hgnc, or dmel
dbDir <- "YOURDIRECTORY" # RcisTarget databases location
myDatasetTitle <- "P0" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=16) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


### ====================================================================================
### ========================== Co-expression network ===================================
### ====================================================================================
# 1. (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
# 2. Before proceeding to the network inference, check whether any known relevant genes are filtered-out
#(if any relevant gene is missing, double-check whether the filters are appropiate):
interestingGenes <- c("Six2", "Cited1", "Nphs1")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

# 3. We can now **filter the expression matrix** to contain only these `r length(genesKept)` genes. 
#This matrix is now ready for the co-expression analysis.
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

# 4. Correlation
#GENIE3/GRNBoost can detect both positive and negative associations. In order to distinguish potential activation from repression, we will split the targets into 
#positive- and negative-correlated targets (i.e. Spearman correlation between the TF and the potential target).
#*(This step can be run either before/after or simultaneously to GENIE3/GRNBoost)*
runCorrelation(exprMat_filtered, scenicOptions)

# 5. Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# 6. Build and score the GRN (runSCENIC_...)
#Once the results from GENIE3/GRNBoost (and the correlation) are ready, the remaining steps of SCENIC can be run. 
#The easiest/fastest way is to use the following *wrapper* functions, each of them corresponding to one of the main steps in SCENIC workflow:

#Build the *gene regulatory network*: 
#1. Get co-expression modules
#2. Get regulons (with [RcisTarget](http://bioconductor.org/packages/RcisTarget)): TF motif analysis)

#Identify *cell states*:
#3. Score GRN (regulons) in the cells (with [AUCell](http://bioconductor.org/packages/AUCell))
#4. Cluster cells according to the GRN activity

#>An overview of the steps workflow is explained in [*@aibar2017*](http://dx.doi.org/10.1038/nmeth.4463).
#Detailed tutorials/notebooks explaining these functions in detail are also available (see `vignette(package="SCENIC")`). These might be useful for users who want to know 
#the details of the implementation, understand the results more in depth, or to modify or run only some of the steps of the workflow. 
#We recommend to check `vignette("detailedStep_2_createRegulons", package="SCENIC")` to understand how the Gene Regulatory Networks are build. 
#For more info on the scoring of the networks on the cells see `vignette("detailedStep_3_scoreCells", package="SCENIC")`.

#Run the remaining steps using the *wrapper* functions: 
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 16
scenicOptions@settings$seed <- 123

### Build and score the GRN
exprMat_filtered <- as.matrix(exprMat_filtered)
exprMat_log <- log2(exprMat_filtered+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Binarize activity?
runSCENIC_4_aucell_binarize(scenicOptions)



### ====================================================================================
### ========================== Exploring output  =======================================
### ====================================================================================
# Check files in folder 'output'

##Plot AUC
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
#savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Show TF expression:
Cairo::CairoPDF("output/Step4_TFexpression(selected).pdf", width=20, height=15)
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_log, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Uncx", "Pax2")],], plots="Expression")
dev.off()

# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

#Density plot to detect most likely stable states (higher-density areas in the t-SNE):
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
Cairo::CairoPDF("output/Step4_Density plot.pdf", width=20, height=15)
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)
dev.off()

#Show several regulons simultaneously:
Cairo::CairoPDF("output/Step4_Regulons.pdf", width=10, height=5)
par(mfrow=c(1,2))
regulonNames <- c( "Ppargc1a","Hnf4a")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6 )
text(0, 10, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-20,-10, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
regulonNames <- list(red=c("Ppargc1a"),
                     blue=c("Hnf4a"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
text(5, 15, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(5, 15-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)
dev.off()

#Average Regulon Activity by cluster (Clusters could also be used instead of “cell types”, e.g. with Seurat: cellInfo <- data.frame(seuratCluster=Idents(seuratObject)))
cellInfo <- readRDS("int/cellInfo.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

Cairo::CairoPDF("output/Step4_regulonActivity_byCellType_Scaled.pdf", height=30, width=10)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
write.table(topRegulators, 'output/Step4_regulonActivity_byCellType_Scaled.txt', sep='\t', quote=F)

#Binarized version (~ percentage of cells of that cell type/cluster with the regulon active)
minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
Cairo::CairoPDF("output/Step4_regulonActivity_byCellType_Binarized.pdf", height=30, width=10)
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5, 
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
write.table(topRegulators, 'output/Step4_regulonActivity_byCellType_Binarized.txt', sep='\t', quote=F)

#Visualizing the regulon activities on embeddings/trajectories calculated with other methods…
dr_coords <- readRDS(file="int/dr_coords.Rds")

tfs <- c("Uncx","Pax2","Esrra")
par(mfrow=c(2,2))
Cairo::CairoPDF("output/Step4_regulonActivity_inSeuratUMAP.pdf", height=10, width=10)
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
dev.off()

#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)



