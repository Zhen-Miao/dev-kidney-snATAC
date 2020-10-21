## distal element prediction using P0 data
library("Seurat")
library("SnapATAC")
library("GenomicRanges")
library("stringi")

## load marker genes

marker_list_4 = readRDS("...")

## load snapATAC object

x.sp2 = readRDS("...")


## do such for all clusters
cluster_k = c("NP","PT","immune","DCT","LOH","stroma1","Endo","IC","PC","Podo")

csum = data.frame(cell = as.character(x.sp2@barcode),sample =  as.character(x.sp2@sample), cluster = as.character(x.sp2@cluster),stringsAsFactors = F)

## do not study the difference within PT, so merge PT with PT2
csum$cluster[csum$cluster == "PT2"] = "PT"
csum$cluster[csum$cluster == "stroma2"] = "stroma1"

csum$cluster[csum$sample == "90028" | csum$sample == "90029"] = paste0(csum$cluster[csum$sample == "90028" | csum$sample == "90029"],"_p0")
x.sp2@cluster = as.factor(csum$cluster)


## to study DAR for P0
cluster_i = cluster_j = paste0(cluster_k,"_p0")

setwd("/home/zhenmiao/kidney_reg")

idy = rep(list(), length = length(cluster_i))
names(idy) = cluster_i
for(i in cluster_i){
  idy[[i]] = rep(list(), length = length(cluster_i))
  names(idy[[i]]) = cluster_i
}

for(i in cluster_i){
  for(j in cluster_j){
    if(i ==j){
      next
    }
    DARs = findDAR(
      obj=x.sp2,
      input.mat="pmat",
      cluster.pos=i,
      cluster.neg=j,
      # cluster.neg.method="random",
      bcv=0.1,
      test.method="exactTest",
      seed.use=10
    );
    DARs$FDR = p.adjust(DARs$PValue, method="BH");
    idy[[i]][[j]] = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  }
}

# saveRDS(object = idy, file = "differentially accessible peaks between p0 clusters.RDS")

idy = readRDS(file = "differentially accessible peaks between p0 clusters.RDS")

# cluster_k = c("PT","PT2","immune","DCT","LOH","stroma1","stroma2","Endo","IC","PC","Podo")


t = rep(list(),length = length(cluster_i))
names(t) = cluster_i
for(i in cluster_i){
  for(j in cluster_i){
    if(i ==j){
      next
    }

    if(length(t[[i]]) == 0){
      t[[i]] = idy[[i]][[j]]
    }else{
      t[[i]] = intersect(t[[i]],idy[[i]][[j]])
    }
  }
}


## get the seqnames for each peak


peaks = rep(list(),length = length(cluster_i))
names(peaks) = cluster_i
for(i in cluster_i){
  peaks[[i]] = x.sp2@peak$name[t[[i]]]
}


## for each marker genes, find the locations


# install.packages("refGenome")
# library(refGenome)

seq.levels = paste0("chr",c(1:19, "X", "Y"))

setwd("/home/zhenmiao/kidney_project_4")

annotation.file = "gencode.vM20.annotation.gtf"
gtf <- rtracklayer::import(con = annotation.file)


gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')

gtf.genes <- gtf[gtf$type == 'gene']

upstream = 100000
downstream = 100000


gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)



setwd("/home/zhenmiao/kidney_reg")
scicar = read.csv("sci-car distal element summary.csv",stringsAsFactors = F)


##############################
## take Podo as an example
##############################

i = cluster_i[1]


setwd("/home/zhenmiao/kidney_reg")

summary_matrix = matrix(nrow = 6,ncol = length(cluster_i))
colnames(summary_matrix) = cluster_i
rownames(summary_matrix) = c("n_DEG","n_DAR","n_predicted_regulatory","n_predicted_distal_regulatory",
                             "n_sci-car_regulatory","n_overlapped_with_scicar")

marker_list = marker_list_4
# names(marker_list) = paste0(names(marker_list),"_p0",sep = "")

for(i in cluster_i){
  ## get the ith cluster peaks
  peak.df = peaks[[i]]
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = as.data.frame(peak.df))
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
  summary_matrix["n_DAR",i] = length(peaks.gr)
  
  ## remove the promoter peaks
  gtf.promoters <- SummarizedExperiment::promoters(x = gtf.genes, upstream = 1000, downstream = 500)
  peaks.gr_nopromoter = subsetByOverlaps(peaks.gr,gtf.promoters,invert = T)
  
  ## also get the subset of gtf.body_prom
  gene.df = marker_list[[i]]
  gtf.body_prom.sub = gtf.body_prom[gtf.body_prom$gene_name %in% gene.df ,]
  gtf.genes.sub = gtf.genes[gtf.genes$gene_name %in% gene.df,]
  summary_matrix["n_DEG",i] = length(gtf.genes.sub)
  
  print(paste0("successfully get the subset data for cell type ",i))
  
  ## get the overlapped peaks
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom.sub)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes.sub[S4Vectors::subjectHits(x = keep.overlaps)]
  summary_matrix["n_predicted_regulatory",i] = length(gene.ids)
  
  ## get the overlapped non-promoter peaks
  gene.distances_nopromoter <- GenomicRanges::distanceToNearest(x = peaks.gr_nopromoter, subject = gtf.body_prom.sub)
  keep.overlaps_nopromoter <- gene.distances[rtracklayer::mcols(x = gene.distances_nopromoter)$distance == 0]
  peak.ids_nopromoter <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps_nopromoter)]
  gene.ids_nopromoter <- gtf.genes.sub[S4Vectors::subjectHits(x = keep.overlaps_nopromoter)]
  summary_matrix["n_predicted_distal_regulatory",i] = length(gene.ids_nopromoter)
  
  print(paste0("successfully get the peak overlap for cell type ",i))
  
  if (length(gene.ids) == 0) {
    next
  }
  
  ## make the proper dataframe
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids_nopromoter$gene.name <- gene.ids_nopromoter$gene_name
  ## also input the gene range
  granges.ids = as.data.frame(gene.ids)
  peak.ids$gene.range <- paste0(granges.ids$seqnames,":",granges.ids$start,"-",granges.ids$end)
  peak.ids <- as.data.frame(x = peak.ids)
  
  ## do the same for non-promoter peaks
  granges.ids_nopromoter = as.data.frame(gene.ids_nopromoter)
  peak.ids_nopromoter$gene.range <- paste0(granges.ids_nopromoter$seqnames,":",granges.ids_nopromoter$start,"-",granges.ids_nopromoter$end)
  peak.ids_nopromoter <- as.data.frame(x = peak.ids_nopromoter)
  
  ## save the data
  write.csv(peak.ids, paste0("linked peak-gene list in cluster",i,".csv",sep = " "))
  write.csv(peak.ids_nopromoter, paste0("linked peak-gene list in cluster",i,"without promoters.csv",sep = " "))
  
  print(paste0("successfully write the files for cell type ",i))
  
  ## check the overlap with scicar
  # scicar.df = scicar[scicar$gene_name %in% gene.df,]
  scicar.df2 = scicar[scicar$gene_name %in% peak.ids$gene.name,]
  
  if(dim(scicar.df2)[1] == 0){
    summary_matrix["n_sci-car_regulatory",i] = 0
    next
  }
  summary_matrix["n_sci-car_regulatory",i] = dim(scicar.df2)[1]
  
  ## sci-car peaks
  scicarpeaks = do.call(what = rbind, args = strsplit(x = gsub(scicar.df2$peak, pattern = ":", replacement = "-"), split = "-"))
  scicarpeaks = as.data.frame(scicarpeaks)
  colnames(x = scicarpeaks) <- c("chromosome", 'start', 'end')
  scicarpeaks <- GenomicRanges::makeGRangesFromDataFrame(df = as.data.frame(scicarpeaks))
  
  ## predicted peaks
  pdpeaks = peak.ids[,c("seqnames","start","end")]
  colnames(pdpeaks) = c("chromosome", 'start', 'end')
  pdpeaks <- GenomicRanges::makeGRangesFromDataFrame(df = pdpeaks)
  
  summary_matrix["n_overlapped_with_scicar",i] = length(subsetByOverlaps(scicarpeaks, pdpeaks))
}


# write.csv(summary_matrix,"summary information for distal regulatory elements for P0.csv")



## to study that for P0 data
kin = readRDS("P1_obj.RDS")



kin.markers <- FindAllMarkers(kin, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.20)



## process the cell type specific markers in the merged data

setwd("/Users/zhenmiao/Dropbox/developmental kidney project/MS folder/Supplementary Tables/")


mar = read.csv("Supplementary Table 1__Cell type-specific gene expression derived from scRNA-seq.csv")

cluster_k = c("NP","PT","PT2","immune","DCT","LOH","stroma1","stroma2","Endo","IC","PC","Podo")


marker_list_2 = rep(list(), length = length(cluster_k))
names(marker_list_2) = cluster_k
marker_list_2[["NP"]] = mar$NP
marker_list_2[["PT"]] = mar$PT.S1
marker_list_2[["PT2"]] = mar$PT.S3
marker_list_2[["immune"]] = union(mar$Macro,mar$Neutro)
marker_list_2[["LOH"]] = mar$LOH
marker_list_2[["DCT"]] = union(mar$DCT,"Slc12a3")
marker_list_2[["Endo"]] = mar$Endo
marker_list_2[["IC"]] = mar$IC
marker_list_2[["PC"]] = mar$PC
marker_list_2[["Podo"]] = mar$Podo
marker_list_2[["stroma1"]] = union(mar$Stroma, mar$Stroma.like)
marker_list_2[["stroma2"]] = marker_list_2[["stroma1"]]


saveRDS(marker_list_2,file = "marker_list_2.rds")
amarker = read.csv("TableS1_cell_type_markers.csv",stringsAsFactors = F)


cluster_k = c("NP","PT","PT2","immune","DCT","LOH","stroma1","stroma2","Endo","IC","PC","Podo")


marker_list = rep(list(), length = length(cluster_k))
names(marker_list) = cluster_k

marker_list[["PT"]] = union(amarker$PT.S1, amarker$PT.S2)
marker_list[["PT2"]] = union(amarker$PT.S2, amarker$PT.S3)

marker_list[["immune"]] = union(amarker$B1,amarker$B2,amarker$CD4.T, amarker$Granul,amarker$Mono, amarker$DC,amarker$Macro,amarker$pDC, amarker$Baso,
                                amarker$Treg, amarker$Th17, amarker$NKT, amarker$CD8.effector, amarker$CD8.T, amarker$NK)
marker_list[["LOH"]] = union(amarker$ALOH,amarker$DLOH)
marker_list[["DCT"]] = union(amarker$DCT,"Slc12a3")
marker_list[["Endo"]] = amarker$Endo
marker_list[["IC"]] = union(amarker$A.IC,amarker$B.IC)
marker_list[["PC"]] = amarker$PC
marker_list[["Podo"]] = amarker$Podo
marker_list[["stroma1"]] = amarker$GEC
marker_list[["stroma2"]] = amarker$GEC

marker_list_2= readRDS(file = "marker_list_2.rds")

marker_list_3 = rep(list(), length = length(cluster_k))
names(marker_list_3) = cluster_k

for(i in cluster_k){
  marker_list_3[[i]] = union(marker_list[[i]], marker_list_2[[i]])
}

marker_list_3[["PT"]] = union(marker_list_3[["PT"]], marker_list_3[["PT2"]])
marker_list_3[["stroma1"]] = union(marker_list_3[["stroma1"]],marker_list_3[["stroma2"]])


cluster_k2 =  c("NP","PT","immune","DCT","LOH","stroma1","Endo","IC","PC","Podo")

marker_list_4 = rep(list(), length = length(cluster_k2))
names(marker_list_4) = cluster_k2

for(i in cluster_k2){
  marker_list_4[[i]] = marker_list_3[[i]]
}

marker_list_4[["NP"]] = union(marker_list_4[["NP"]], np_additional)

saveRDS(marker_list_4,"marker_list_for_P0.rds")



