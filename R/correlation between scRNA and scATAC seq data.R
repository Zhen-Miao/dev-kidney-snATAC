## pairwise comparison of scATAC and scRNA data

library("Seurat")

pbmc.merge = readRDS("...")

x.sp2 = readRDS("...")


## merge the peaks by cluster
atac.cluster = unique(x.sp2@cluster)

pcmatrix = matrix(nrow = length(x.sp2@peak),ncol = length(atac.cluster))
rownames(pcmatrix) = x.sp2@peak$name
colnames(pcmatrix) = atac.cluster

for(i in atac.cluster){
  p = t(x.sp2@pmat[(x.sp2@sample == "90028" | x.sp2@sample == "90029") & 
                     x.sp2@cluster == i,])
  pcmatrix[,i] = rowSums(as.matrix(p))
  print(i)
}


pcmatrix2 = matrix(nrow = length(x.sp2@peak),ncol = length(atac.cluster))
rownames(pcmatrix2) = x.sp2@peak$name
colnames(pcmatrix2) = atac.cluster

for(i in atac.cluster){
  p = t(x.sp2@pmat[(x.sp2@sample == "90025" | x.sp2@sample == "90026" | x.sp2@sample == "batch_1" ) & 
                     x.sp2@cluster == i,])
  pcmatrix2[,i] = rowSums(as.matrix(p))
  print(i)
}




table(x.sp2@cluster,x.sp2@sample)

setwd("/home/zhenmiao/kidney_project_4")
# create a gene activity matrix from the peak matrix and GTF, using chromosomes 1:22, X, and Y.
# Peaks that fall within gene bodies, or 2kb upstream of a gene, are considered

activity.matrix <- CreateGeneActivityMatrix(peak.matrix = pcmatrix, annotation.file = "gencode.vM20.annotation.gtf",
                                            seq.levels = paste0("chr",c(1:19, "X", "Y")), upstream = 2000, verbose = TRUE)
activity.matrix2 <- CreateGeneActivityMatrix(peak.matrix = pcmatrix2, annotation.file = "gencode.vM20.annotation.gtf",
                                            seq.levels = paste0("chr",c(1:19, "X", "Y")), upstream = 2000, verbose = TRUE)
activity.matrix2.save = activity.matrix2
activity.matrix2.noXY <- CreateGeneActivityMatrix(peak.matrix = pcmatrix2, annotation.file = "gencode.vM20.annotation.gtf",
                                             seq.levels = paste0("chr",c(1:19)), upstream = 2000, verbose = TRUE)
# 
# 
# saveRDS(activity.matrix, "activity.matrix.P0.rds")
# saveRDS(activity.matrix2, "activity.matrix2.adult.rds")
# saveRDS(activity.matrix2.save, "activity.matrix2.adult updated.rds")

# activity.matrix = readRDS("activity.matrix.P0.rds")
# activity.matrix2 = readRDS("activity.matrix2.adult updated.rds")
# activity.matrix2 = readRDS("activity.matrix2.adult updated no XY.rds")
# 

### get the matrix for scRNA-seq data P0
rna.cluster = unique(pbmc.merge@active.ident)
gcmatrix = matrix(nrow = length(rownames(pbmc.merge[["RNA"]])), ncol = length(rna.cluster))
rownames(gcmatrix) = rownames(pbmc.merge[["RNA"]])
colnames(gcmatrix) = rna.cluster

for(i in rna.cluster){
  p = pbmc.merge[["RNA"]][,pbmc.merge@active.ident == i & pbmc.merge@meta.data$orig.ident == "P0_batch2"]
  gcmatrix[,i] = rowSums(as.matrix(p))
  print(i)
}
  

var.genes = VariableFeatures(pbmc.merge)
var.genes2 = intersect(var.genes,row.names(activity.matrix))



## we can also do normalization and scale, then do the comparison
atac.cluster2 = c("NP","Podo","PT","PT2","LOH","DCT","PC","IC","stroma1","Endo","immune")
# rna.cluster2 = c("6","0","12","2","1","4","3","8","9","7","10")
# rna.cluster3 = c("14","11","3","1","9","8","5","12","13","6","15")

rna.cluster3 = c("14","11","3","1","9","0","5","12","13","6","15")


# var.genes3 = setdiff(var.genes2, c("Isx","Gm2061", "Gm21885", "Otc"))

activity.matrix_p0 = activity.matrix

activity.matrix_p0[,"stroma1"] = activity.matrix_p0[,"stroma1"] + activity.matrix_p0[,"stroma2"] 
activity.matrix_p0 = activity.matrix_p0[,atac.cluster2]

gcmatrix_p0 = gcmatrix
gcmatrix_p0[,"15"] = gcmatrix_p0[,"15"] + gcmatrix_p0[,"10"]
gcmatrix_p0[,"0"] = gcmatrix_p0[,"0"] + gcmatrix_p0[,"16"]

gcmatrix_p0 = gcmatrix_p0[,rna.cluster3]

activity.matrix_p0 = as.matrix(activity.matrix_p0)
gcmatrix_p0 = as.matrix(gcmatrix_p0)

activity.matrix_p0 = apply(activity.matrix_p0, 2, function(x) (x/sum(x))*10000)
gcmatrix_p0 = apply(gcmatrix_p0, 2, function(x) (x/sum(x))*10000)

activity.matrix_p0 = activity.matrix_p0[var.genes2,]
gcmatrix_p0 = gcmatrix_p0[var.genes2,]

for(i in rownames(activity.matrix_p0)){
  activity.matrix_p0[i,] = (activity.matrix_p0[i,] - mean(activity.matrix_p0[i,]) )/sd(activity.matrix_p0[i,])
}

for(i in rownames(gcmatrix_p0)){
  gcmatrix_p0[i,] = ( gcmatrix_p0[i,] - mean(gcmatrix_p0[i,]) ) / sd(gcmatrix_p0[i,]) 
}

i = "NP"
t1 = names(activity.matrix_p0[is.na(activity.matrix_p0[,i]),i])

i = "14"
t2 = names(gcmatrix_p0[is.na(gcmatrix_p0[,i]),i]) 

var.genes3 = setdiff(var.genes2,t1)
var.genes3 = setdiff(var.genes3,t2)

  
  corr_mat2 = matrix(nrow = length(rna.cluster3), ncol = length(atac.cluster2))
  rownames(corr_mat2) = rna.cluster3
  colnames(corr_mat2) = atac.cluster2
  
  
  for(i in rna.cluster3){
    for(j in atac.cluster2){
      corr_mat2[i,j] = cor(gcmatrix_p0[var.genes3,i], activity.matrix_p0[var.genes3,j], method = "pearson")
    }
    print(i)
  }
  
  true_cluster_name = c("NP","Podo","PT S1","PT S3","LOH","DCT","PC","IC","stroma","Endo","immune")
rownames(corr_mat2) <- colnames(corr_mat2) <- true_cluster_name
  

x = corr_mat2

# saveRDS(x, "scRNA-scATAC heatmap P0 updated DCT data.rds")

my_palette <-colorRampPalette(c("blue","white","red"))(50)

# setwd("/Users/zhenmiao/Dropbox/developmental kidney project/data/")

# x = readRDS("scRNA-scATAC heatmap P0 updated DCT data.rds")

# write.csv(x,"scRNA-scATAC correlation values.csv")

library(gplots)

pdf("scRNA-scATAC heatmap P0 updated DCT.pdf")
heatmap.2(x[1:dim(x)[1],], col = my_palette, keysize = 1, density.info="none", trace="none",cexRow=1.5,cexCol = 1.5,
          breaks=c(seq(-0.7,0.7, length=51) ),Rowv=F, Colv=F, dendrogram="none", margins = c(6,6),
          colsep=0:(ncol(x)+1),rowsep=0:(nrow(x)+1),sepcolor="black", sepwidth=c(0.01,0.01))
dev.off()

svg("scRNA-scATAC heatmap P0 updated DCT.svg")

heatmap.2(x[1:dim(x)[1],], col=my_palette, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1.5,cexCol = 1.5,
          breaks=c(seq(-0.7,0.7, length=51) ),Rowv=F, Colv=F, dendrogram="none", margins = c(6,6),
          colsep=0:ncol(x), rowsep=0:nrow(x), sepcolor="black", sepwidth=c(0.01,0.01))
dev.off()

  
  ### get the matrix for scRNA-seq data adult
# activity.matrix2 = readRDS("activity.matrix2.adult.rds")
  rna.cluster2 = unique(pbmc.merge@active.ident)
  gcmatrix2 = matrix(nrow = length(rownames(pbmc.merge[["RNA"]])), ncol = length(rna.cluster2))
  rownames(gcmatrix2) = rownames(pbmc.merge[["RNA"]])
  colnames(gcmatrix2) = rna.cluster2
  
  for(i in rna.cluster2){
    p = pbmc.merge[["RNA"]][,pbmc.merge@active.ident == i & pbmc.merge@meta.data$orig.ident == "Ksp"]
    gcmatrix2[,i] = rowSums(as.matrix(p))
    print(i)
  }

  
  
  var.genes3 = VariableFeatures(pbmc.merge)
  var.genes4 = intersect(var.genes3,row.names(activity.matrix2))
  atac.cluster2 = c("Podo","PT","PT2","LOH","DCT","PC","IC","stroma1","Endo","immune")
  rna.cluster2 = c("11","3","1","9","0","5","12","13","6","15")
  
  true_cluster_name = c("Podo","PT S1","PT S3","LOH","DCT","PC","IC","stroma","Endo","immune")
  # var.genes4 = setdiff(var.genes4, c("Isx","Gm2061","Gm21885","Otc","Gzmk","Ifngas1","Il2ra","H2-Q1","Gm14085","Ceacam18","Cdx2" ))
  
  activity.matrix2[,"stroma1"] = activity.matrix2[,"stroma1"] + activity.matrix2[,"stroma2"] 
  activity.matrix2 = activity.matrix2[,atac.cluster2]
  
  
  gcmatrix2[,"15"] = gcmatrix2[,"15"] + gcmatrix2[,"10"]
  gcmatrix2[,"0"] = gcmatrix2[,"0"] + gcmatrix2[,"16"]
  gcmatrix2 = gcmatrix2[,rna.cluster2]
  
  activity.matrix2 = as.matrix(activity.matrix2)
  gcmatrix2 = as.matrix(gcmatrix2)
  
  activity.matrix2 = apply(activity.matrix2, 2, function(x) (x/sum(x))*10000)
  gcmatrix2 = apply(gcmatrix2, 2, function(x) (x/sum(x))*10000)
  
  activity.matrix2 = activity.matrix2[var.genes4,]
  gcmatrix2 = gcmatrix2[var.genes4,]
  
  sum.gcmatrix2 = rowSums(gcmatrix2)
  genes.to.exclude = names(sum.gcmatrix2[sum.gcmatrix2 == 0])
  
  var.genes4 = setdiff(var.genes4, genes.to.exclude)
  
  activity.matrix2 = activity.matrix2[var.genes4,]
  gcmatrix2 = gcmatrix2[var.genes4,]
  
  
  for(i in rownames(activity.matrix2)){
    activity.matrix2[i,] = (activity.matrix2[i,] - mean(activity.matrix2[i,]) )/sd(activity.matrix2[i,])
  }
  
  for(i in rownames(gcmatrix2)){
    gcmatrix2[i,] = ( gcmatrix2[i,] - mean(gcmatrix2[i,]) ) / sd(gcmatrix2[i,]) 
  }
  
  
  
  corr_mat2 = matrix(nrow = length(rna.cluster2), ncol = length(atac.cluster2))
  rownames(corr_mat2) = rna.cluster2
  colnames(corr_mat2) = atac.cluster2
  

  
  # sum.activity.matrix2 = rowSums(activity.matrix2)
  # sum.activity.matrix2[is.na(sum.activity.matrix2)]
  # 
  # sum.gcmatrix2 = rowSums(gcmatrix2)
  # sum.gcmatrix2[is.na(sum.gcmatrix2)]
  # 
  for(i in rna.cluster2){
    for(j in atac.cluster2){
      corr_mat2[i,j] = cor(gcmatrix2[,i], activity.matrix2[,j], method = "pearson")
    }
    print(i)
  }
  
  colnames(corr_mat2) <- rownames(corr_mat2) <- true_cluster_name
  
  saveRDS(corr_mat2, "scRNA-scATAC heatmap adult updated DCT data noXY.rds")
  saveRDS(corr_mat2, "scRNA-scATAC heatmap adult updated DCT data wXY.rds")
  
  write.csv(corr_mat2,"scRNA-scATAC correlation matrix using all variable genes after scale adult.csv")
  
  x = corr_mat2
  
  pdf("scRNA-scATAC heatmap adult updated DCT.pdf")
  heatmap.2(x[1:dim(x)[1],], col = my_palette, keysize = 1, density.info="none", trace="none",cexRow=1.5,cexCol = 1.5,
            breaks=c(seq(-0.7,0.7, length=51) ),Rowv=F, Colv=F, dendrogram="none", margins = c(6,6),
            colsep=0:(ncol(x)+1),rowsep=0:(nrow(x)+1),sepcolor="black", sepwidth=c(0.01,0.01))
  dev.off()
  
  svg("scRNA-scATAC heatmap adult updated DCT.svg")
  
  heatmap.2(x[1:dim(x)[1],], col=my_palette, keysize = 1, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1.5,cexCol = 1.5,
            breaks=c(seq(-0.7,0.7, length=51) ),Rowv=F, Colv=F, dendrogram="none", margins = c(6,6),
            colsep=0:(ncol(x)+1),rowsep=0:(nrow(x)+1), sepcolor="black", sepwidth=c(0.01,0.01))
  dev.off()
  
  
  
  
  # saveRDS(activity.matrix, "activity.matrix_by_P0_cell_types.rds")
  
  setwd("/Users/zhenmiao/Dropbox/developmental kidney project/data/")
  corr_mat2 = read.csv("scRNA-scATAC correlation matrix using all variable genes after scale adult.csv",header = T)
  rownames(corr_mat2) = corr_mat2$X
  corr_mat2 = corr_mat2[,c("Podo","PT","PT2","LOH","DCT","PC","IC","stroma1","Endo","immune")]
  
  
  
  corr_mat2 = as.matrix(corr_mat2)
  
  corr_mat2 = apply(corr_mat2, 2, as.numeric)
  colnames(corr_mat2) <- rownames(corr_mat2) <- true_cluster_name
  
  pdf("adult scRNA-scATAC heatmap.pdf")
  
  dev.off()
  
  
  x = corr_mat2
  my_palette <-colorRampPalette(c("blue","white","red"))(256)
  
  pdf("scRNA-scATAC heatmap.pdf")
  
  heatmap.2(x[1:dim(x)[1],], col=my_palette, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, 
            breaks=c(seq(-0.7,0.7, length=257) ),Rowv=F, Colv=F, dendrogram="none", 
            colsep=1:ncol(x), rowsep=1:nrow(x), sepcolor="black", sepwidth=c(0.01,0.01))
  dev.off()
  
  
  
  corr_mat1 = read.csv("scRNA-scATAC correlation matrix using all variable genes after scale P0.csv")
  rownames(corr_mat1) = corr_mat1$X
  corr_mat1 = corr_mat1[,c("NP","Podo","PT","PT2","LOH","DCT","PC","IC","stroma1","Endo","immune")]
  
  
  
  corr_mat2 = as.matrix(corr_mat2)
  
  corr_mat2 = apply(corr_mat2, 2, as.numeric)
  colnames(corr_mat2) <- rownames(corr_mat2) <- true_cluster_name
  
  pdf("adult scRNA-scATAC heatmap.pdf")
  
  dev.off()
  
  
  x = corr_mat2
  my_palette <-colorRampPalette(c("blue","white","red"))(256)
  
  pdf("scRNA-scATAC heatmap.pdf")
  
  heatmap.2(x[1:dim(x)[1],], col=my_palette, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, 
            breaks=c(seq(-0.7,0.7, length=257) ),Rowv=F, Colv=F, dendrogram="none", 
            colsep=1:ncol(x), rowsep=1:nrow(x), sepcolor="black", sepwidth=c(0.01,0.01))
  dev.off()
  
  
  
# 
# write.csv(corr_mat,"scRNA-scATAC correlation matrix using all variable genes.csv")
# write.csv(corr_mat2,"scRNA-scATAC correlation matrix using all variable genes after scale.csv")
# 
# 
# saveRDS(activity.matrix, "activity.matrix_by_P0_cell_types.rds")



## get the shared features
share_genes = intersect(rownames(pbmc.atac[["ACTIVITY"]]), rownames(pbmc.merge[["RNA"]]))

## get the highly variable genes only

DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)

var.gene.int = intersect(VariableFeatures(pbmc.atac),VariableFeatures(pbmc.merge))


