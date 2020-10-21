## conduct pair-wise peak comparison

library("SnapATAC");
library("Seurat")
library(viridisLite);
library(ggplot2);


x.sp2 = readRDS("all_comb snapATAC after batch correction with peak info new clu.RDS")


## do such for all clusters

cluster_i = unique(x.sp2@cluster,stringAsFactors=F)
cluster_j = unique(x.sp2@cluster,stringAsFactors=F)



idy = rep(list(), length = 15)
names(idy) = cluster_i
for(i in cluster_i){
  
  idy[[i]] = rep(list(), length = 15)
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
      # test.method="exactTest",
      test.method="",
      seed.use=10
    );
    DARs$FDR = p.adjust(DARs$PValue, method="BH");
    idy[[i]][[j]] = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  }
}

saveRDS(object = idy, file = "differentially accessible peaks.RDS")


idy = readRDS(file = "differentially accessible peaks.RDS")
cluster_k = c("NP","Podo","PT","PT2","LOH","DCT","PC","IC","Endo","immune","stroma1","stroma2")

cluster_j = cluster_k
cluster_i = cluster_k
t = rep(list(),length = length(cluster_j))
names(t) = cluster_j
for(i in cluster_i){
  for(j in cluster_j){
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



## motif detection
sp_clusters = c("PT2","immune","stroma1","Endo","stroma2","Podo","NP")

motifs = rep(list(), length = length(sp_clusters))
names(motifs) = sp_clusters

motifs_summary = lapply(levels(sp_clusters), function(sp_clusters){
  runHomer(
    x.sp2[,t[[sp_clusters]],"pmat"], 
    mat = "pmat",
    path.to.homer = "...",
    result.dir = paste("...",sp_clusters,sep = ""),
    num.cores=5,
    genome = 'mm10',
    motif.length = 10,
    scan.size = 300,
    optimize.count = 2,
    background = 'automatic',
    local.background = FALSE,
    only.known = FALSE,
    only.denovo = FALSE,
    fdr.num = 5,
    cache = 100,
    overwrite = TRUE,
    keep.minimal = FALSE
  )
})


for(i in sp_clusters){
  runHomer(
    x.sp2[,t[[i]],"pmat"], 
    mat = "pmat",
    path.to.homer = "...",
    result.dir = paste("...",i,sep = ""),
    num.cores=8,
    genome = 'mm10',
    motif.length = 10,
    scan.size = 300,
    optimize.count = 2,
    background = 'automatic',
    local.background = FALSE,
    only.known = FALSE,
    only.denovo = FALSE,
    fdr.num = 5,
    cache = 100,
    overwrite = TRUE,
    keep.minimal = FALSE
  )
}
  


names(motifs_summary) = sp_clusters

saveRDS(motifs_summary,file = "motif_summary.RDS")


