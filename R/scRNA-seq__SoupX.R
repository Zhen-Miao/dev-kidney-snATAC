library(SoupX)

#========= Loading the data ==========
#sc = load10X("/outs")
toc = Seurat::Read10X("/outs/filtered_feature_bc_matrix")
tod = Seurat::Read10X("/outs/raw_feature_bc_matrix")

#rename tod barcodes
colnames(tod)
newcol_tod <- paste0("Ksp_",colnames(tod)) #add Ksp_ to barcode
newcol2_tod <- sub("\\-.*", "", newcol_tod) #keep only the barcode without the "-1"
tod2 <- tod
colnames(tod2) <- newcol2_tod

#rename toc barcodes to match final seurat barcodes
colnames(toc)
newcol <- paste0("Ksp_",colnames(toc)) #add Ksp_ to barcode
newcol2 <- sub("\\-.*", "", newcol) #keep only the barcode without the "-1"
toc2 <- toc
colnames(toc2) <- newcol2

#FILTER TOC for remaining final barcodes
P0_adult <- readRDS("/P0_adult_final.rds")
bc <- colnames(P0_adult)
bc2 <- bc[grepl("Ksp_", bc)] #keep only Ksp barcodes
toc3 <- toc2[,bc2] #this is the new toc matrix that only keeps the barcodes from the final Seurat object

sc = SoupChannel(tod2, toc3)

#========= Adding extra meta data to the SoupChannel object ==========
P0_adult@meta.data$seurat_clusters2 <- plyr::mapvalues(
  x = P0_adult@meta.data$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
           "16", "17"),
  to = c("DCT",
         "PST",
         "Stroma 2",
         "PCT",
         "PST",
         "PC",
         "Endo",
         "Proliferating",
         "Early PT",
         "LOH",
         "Macro",
         "Podo",
         "IC",
         "Stroma 1", #Stroma
         "NP",
         "Neutro",
         "DCT",
         "Unnamed")
)
Idents(P0_adult) <- P0_adult@meta.data$seurat_clusters2
#sort
P0_adult@active.ident <- factor(P0_adult@active.ident, 
                                levels=c("NP",
                                         "Proliferating",
                                         "Podo",
                                         "Early PT",
                                         "PCT",
                                         "PST",
                                         "PT",
                                         "LOH",
                                         "DCT",
                                         "PC",
                                         "IC",
                                         "Endo",
                                         "Macro",
                                         "Neutro",
                                         "Stroma 1",
                                         "Stroma 2",
                                         "Unnamed"))
table(Idents(P0_adult))

metadata <- data.frame(Idents(P0_adult), P0_adult@reductions$umap@cell.embeddings)
metadata$barcode <- rownames(metadata)
metadata2 <- data.frame(metadata[grepl('Ksp_',metadata$barcode),]) #keep only rows of Ksp_ cells
colnames(metadata2) <- c("cluster", "UMAP_1", "UMAP_2", "barcode") #rename colnames of dataframe metadata2
head(metadata2)
tail(metadata2)

#add metadata to sc object
sc = setClusters(sc, setNames(metadata2$cluster, rownames(metadata2)))
sc = setDR(sc, metadata2[colnames(sc$toc), c("UMAP_1", "UMAP_2")])
str(sc)

#========= Visual sanity checks ==========
library(ggplot2)
dd = metadata2 #needed to change this
rownames(dd) <- colnames(sc$toc) #needed to add this
mids = aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = dd, FUN = mean)
gg = ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = cluster), size = 0.2) + 
  geom_label(data = mids, aes(label = cluster)) + ggtitle("Ksp cluster annotation") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))

pdf('/Ksp_plot_UMAP.pdf')
plot(gg)
dev.off()

dd$Aqp1 = sc$toc["Aqp1", ]
dd$Aqp2 = sc$toc["Aqp2", ]
dd$Aqp3 = sc$toc["Aqp3", ]
gg = ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = Aqp1 > 0))
pdf('/Ksp_plot_UMAP_Aqp1.pdf')
plot(gg)
dev.off()

gg = ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = Aqp2 > 0))
pdf('/Ksp_plot_UMAP_Aqp2.pdf')
plot(gg)
dev.off()

gg = ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = Aqp3 > 0))
pdf('/Ksp_plot_UMAP_Aqp3.pdf')
plot(gg)
dev.off()

#========= Estimating the contamination fraction ==========
pdf('/Ksp_plot_autoEstCont.pdf')
sc = autoEstCont(sc)
dev.off()

#========= Correcting expression profile ==========
out = adjustCounts(sc)

#=====Investigating changes in expression
cntSoggy = rowSums(as.array((sc$toc > 0)))
cntStrained = rowSums(as.array(out > 0))
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed
tail(sort(rowSums(as.array(sc$toc > out))/rowSums(as.array(sc$toc > 0))), n = 20)

#=====Visualising expression distribution
pdf('/Ksp_plot_ChangeMap_Aqp1.pdf')
plotChangeMap(sc, out, "Aqp1")
dev.off()

pdf('/Ksp_plot_ChangeMap_Aqp2.pdf')
plotChangeMap(sc, out, "Aqp2")
dev.off()

pdf('/Ksp_plot_ChangeMap_Aqp3.pdf')
plotChangeMap(sc, out, "Aqp3")
dev.off()

#========= Integrating with downstream tools ==========
DropletUtils:::write10xCounts("/Ksp_strainedCounts", out)





