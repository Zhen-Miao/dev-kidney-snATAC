## study the peak difference between adult and P0


library("SnapATAC");
library("Seurat")
library(viridisLite);
library(ggplot2);
library("RColorBrewer")
library(reshape2)

x.sp2 = readRDS("...")

# ## run it in parallel
# library(BiocParallel)
# register(MulticoreParam(8))

new_cluster = x.sp2@cluster
bat = x.sp2@sample
library(plyr)
bat = revalue(bat,c('90025' = 'p56', '90026' = 'p21', '90028' = 'p0', '90029' = 'p0', 'batch_1' = 'p56'))
new_cluster = revalue(new_cluster, c("PT"="PCT", "PT2"="PST"))
nn_cluster = paste0(new_cluster,'_',bat)
x.sp2@cluster = as.factor(nn_cluster)



clusters = c("PST","PCT","immune","DCT","LOH","stroma1","Endo","IC","stroma2","PC","Podo","NP")     

comp = c('p0_p21', 'p21_p0','p21_p56', 'p56_p21', 'p0_p56', 'p56_p0')


idy = rep(list(), length = 6)
names(idy) = comp
for(i in comp){
  idy[[i]] = rep(list(), length = 12)
  names(idy[[i]]) = clusters
}

for (t in comp) {
  if(t == 'p0_p21'){
    cluster_i = paste0(clusters,'_p0')
    cluster_j = paste0(clusters,'_p21')
  }else if(t == 'p21_p0'){
    cluster_i = paste0(clusters,'_p21')
    cluster_j = paste0(clusters,'_p0')
  }else if(t == 'p21_p56'){
    cluster_i = paste0(clusters,'_p21')
    cluster_j = paste0(clusters,'_p56')
  }else if(t == 'p56_p21'){
    cluster_i = paste0(clusters,'_p56')
    cluster_j = paste0(clusters,'_p21')
  }else if(t == 'p0_p56'){
    cluster_i = paste0(clusters,'_p0')
    cluster_j = paste0(clusters,'_p56')
  }else if(t == 'p56_p0'){
    cluster_i = paste0(clusters,'_p56')
    cluster_j = paste0(clusters,'_p0')
  }
  for(k in 1:length(clusters)){
    DARs = findDAR(
      obj=x.sp2,
      input.mat="pmat",
      cluster.pos=cluster_i[k],
      cluster.neg=cluster_j[k],
      # cluster.neg.method="random",
      bcv=0.1,
      test.method="exactTest",
      seed.use=10
    );
    DARs$FDR = p.adjust(DARs$PValue, method="BH");
    m = clusters[k]
    idy[[t]][[m]] = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  }
}

## save the comparison results
saveRDS(idy, 'DAR btw stages.rds')


summary_matrix = matrix(nrow = length(clusters), ncol = length(comp))
rownames(summary_matrix) <- clusters
colnames(summary_matrix) <- comp

for (i in clusters) {
  for(j in comp){
    summary_matrix[i,j] <- length(idy[[j]][[i]])
  }
}



