## combining all 5 batches 

library("SnapATAC");
library("Seurat")
library(viridisLite);
library(ggplot2);


file.list = c("90025.snap", "90026.snap","90028.snap","90029.snap","batch1.snap");
sample.list = c("90025", "90026","90028","90029","batch_1");
x.sp.ls = lapply(seq(file.list), function(i){
  x.sp = createSnap(file=file.list[i], sample=sample.list[i], num.cores=8);
  x.sp
})
names(x.sp.ls) = sample.list;
sample.list


## based on barcode information, filter the cells
barcode.file.list = c("90025_cell_barcode_list_for_combining_Nov4.txt", "90026_cell_barcode_list_for_combining_Nov4.txt",
                      "90028_cell_barcode_list_for_combining_Nov4.txt","90029_cell_barcode_list_for_combining_Nov4.txt",
                      "batch_1_cell_barcode_list_for_combining_Nov4.txt");
barcode.list = lapply(barcode.file.list, function(file){
  read.table(file)[,1];
})
x.sp.list = lapply(seq(x.sp.ls), function(i){
  x.sp = x.sp.ls[[i]];
  x.sp  = x.sp[x.sp@barcode %in% barcode.list[[i]],];
})
names(x.sp.list) = sample.list;


## add bmat to the data
x.sp.list = lapply(seq(x.sp.list), function(i){
  x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=5000, num.cores=10);
  x.sp
})
x.sp.list


## combine the two batches
bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) x.sp@feature$name));
x.sp.list <- lapply(x.sp.list, function(x.sp){
  idy = match(bin.shared, x.sp@feature$name);
  x.sp[,idy, mat="bmat"];
})
x.sp = Reduce(snapRbind, x.sp.list);
rm(x.sp.list); # free memory
gc();
table(x.sp@sample);

x.sp = makeBinary(x.sp, mat="bmat");


## remove the blacklist regions
library(GenomicRanges);
black_list = read.table("mm10.blacklist.bed.gz");
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp
## 546103 bins left

## remove the Mitochrondria reads
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp
## 545183 bins left


## remove the highly conserved reads
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp
## 484606 bins left

## run dimension reduction
set.seed(7732)
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
);


## Determine significant components
pdf("All_comb dim reduction scatter plot.pdf",width = 14,height = 14)
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=10000,
  pdf.file.name=NULL, 
  pdf.height=14, 
  pdf.width=14
);
dev.off()

#################################################
################ decide the number of PCs to use
#################################################

## step 9 graph based clustering 

x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
);
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=834
);
x.sp@metaData$cluster = x.sp@cluster;

## step 10 visualization 


x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="umap",
  seed.use=10
);

pdf("All_comb clustering result using 20 PCs with UMAP.pdf",width = 14,height = 14)

plotViz(
  obj=x.sp,
  method="umap", 
  main="All_batch Cluster",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  #down.sample=10000,
  legend.add=FALSE
);

dev.off()


#####################################
#######also try tsne
############################
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="Rtsne",
  seed.use=10
);

pdf("All_comb clustering result using 20 PCs with tsne.pdf",width = 14,height = 14)
plotViz(
  obj=x.sp,
  method="tsne", 
  main="P0 kidney Cluster",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  #down.sample=10000,
  legend.add=FALSE
);

dev.off()



pdf("All_comb snapATAC without batch correction visualized by tsne.pdf",width = 14,height = 14)

plotViz(
  obj=x.sp,
  method="tsne", 
  main="All_comb without batch correction",
  point.color=x.sp@sample, 
  point.size=0.3, 
  text.add= FALSE,
  # down.sample=10000,
  legend.add=TRUE
);

dev.off()



pdf("All_comb snapATAC without batch correction.pdf visualized by UMAP.pdf",width = 14,height = 14)

plotViz(
  obj=x.sp,
  method="umap", 
  main="All_comb without batch correction",
  point.color=x.sp@sample, 
  point.size=0.3, 
  text.add= FALSE,
  # down.sample=10000,
  legend.add=TRUE
);

dev.off()


#########################################
##############correct for batch effect
######################################
library(harmony);
x.after.sp = runHarmony(
  obj=x.sp, 
  eigs.dim=1:20, 
  meta_data=x.sp@sample # sample index
);

x.after.sp = runKNN(
  obj= x.after.sp,
  eigs.dim=1:20,
  k=15
);

x.after.sp = runCluster(
  obj=x.after.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  path.to.snaptools=NULL,
  seed.use=10
);
x.after.sp@metaData$cluster = x.after.sp@cluster;






x.after.sp = runViz(
  obj=x.after.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="umap",
  seed.use=10
);

pdf("All_comb clustering after batch_correction using 20 PCs with UMAP.pdf",width = 14,height = 14)

plotViz(
  obj=x.after.sp,
  method="umap", 
  main="All_batch Cluster after batch_correction",
  point.color=x.after.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  #down.sample=10000,
  legend.add=FALSE
);

dev.off()


#####################################
#######also try tsne
############################
x.after.sp = runViz(
  obj=x.after.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="Rtsne",
  seed.use=10
);

pdf("All_comb clustering after batch_correction using 20 PCs with tsne.pdf",width = 14,height = 14)
plotViz(
  obj=x.after.sp,
  method="tsne", 
  main="All_comb after batch correction Cluster",
  point.color=x.after.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  #down.sample=10000,
  legend.add=FALSE
);

dev.off()


pdf("All_comb snapATAC after batch_correction colored by batch tsne.pdf",width = 14,height = 14)

plotViz(
  obj=x.after.sp,
  method="tsne", 
  main="All_comb after batch correction",
  point.color=x.after.sp@sample, 
  point.size=0.5, 
  text.add= FALSE,
  # down.sample=10000,
  legend.add=TRUE
);

dev.off()



pdf("All_comb snapATAC after batch_correction colored by batch UMAP.pdf",width = 14,height = 14)

plotViz(
  obj=x.after.sp,
  method="umap", 
  main="All_comb after batch correction",
  point.color=x.after.sp@sample, 
  point.size=0.5, 
  text.add= FALSE,
  # down.sample=10000,
  legend.add=TRUE
);

dev.off()



saveRDS(object = x.after.sp, file = "all_comb snapATAC after batch_correction.RDS")

saveRDS(object = x.sp, file = "all_comb snapATAC before batch correction.RDS")




######


## step 11 gene based annotation

# system("wget http://renlab.sdsc.edu/r3fang/share/github/Mouse_Brain_10X/gencode.vM16.gene.bed");
genes = read.table("gencode.vM16.gene.bed");
library(GenomicRanges);
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4]);
## there are 53379 genes in this dataset



marker_genes_new2 = readRDS(file = "marker_genes_new2.RDS")
cell_type <- names(marker_genes_new2)
## re-plot the genes by the clusters

marker_genes_new3 = rep(list(),length = length(cell_type))
names(marker_genes_new3) <- cell_type 

for(i in cell_type){
  marker_gene_cell_type = intersect(marker_genes_new2[[i]], genes.gr$name)
  marker_genes_new3[[i]] = marker_gene_cell_type[1:30]
}

marker.genes1 = unlist(marker_genes_new3)



marker_gene_science_paper = c(
  "Nrp1","Kdr","Nphs1","Npsh2","Slc27a2","Lrp2","Slc12a1","Umod","Slc12a3","Pvalb","Aqp2","Hsd11b2",
  "Atp6v1g3","Atp6v0d2","Insrr","Rhbg","Mki67","Cdca3","Plac8","S100a4","C1qa","C1qb","S100a8",
  "S100a9","Cd79a","Cd79b","Ltb","Cxcr6","Gzma","Nkg7","Stmn1"
)
## 31 genes

marker_gene_development_paper = c(
  "Six2", "Cited1","Pax2","Crym","Meox1","Traf1", "Uncx", "Eya1", "Spock2",                    ## nephron progenitor
  "Penk", "Nts", "Acta2", "Cldn11", "Tagln1", "Alx1",                                          ## Stroma -Collecting duct-associated
  "Ren1", "Fibin", "Mgp", "Hic1", "Igfbp5", "Rgs5", "Fhl2", "Ntn1", "Lhfp", "Foxd1", "Gdnf",   ## Stroma -cortical stroma
  "S100g", "Tmem52b", "Ly6a", "Sostdc1", "Slc12a1", "Wfdc2", "Mal", "Aqp2",                    ## distal nephron
  "Igfbp3", "Fbln5", "Acta2", "Mgp", "Dcn", "Ace2", 
  "Cfh", "Col3a1", "Vegfd", "Col1a1", "Ndufa4l2", "Rgs5",                                      ## Stroma -medullary stroma
  "Sult1d1", "Spink1", "Aldob", "Hdc", "Pdzk1", "Slc34a1", "Fut9",
  "Fxyd2", "Osr2", "Slc39a5", "Keg1", "Cpn1", "Ttc36", "Ly6a",                                 ## early proximal tubule
  "Lhx1", "Pcp4", "Cldn5", "Sfrp2", "Osr2", "Clec18a", "Clu", "Uncx",
  "Npy", "Ccnd1", "Pax8", "Wnt4", "Mafb", "Sox11", "Jag1",                                     ## S-shaped body
  "Fam132a", "Wnt4", "Tmem100", "Bmper", "Pax2", "Eya1", 
  "Fam107a", "Wt1", "Frzb", "Gxylt2", "Kazald1", "Mycn", "Snap91",                             ## renal vesicle
  "Aldob", "Ttc36", "Spp2", "Kap", "Fxyd2", "Slc34a1", "Fbp1",
  "Gsta2", "Sult1d1", "Spink1", "Ass1", "Lrp2", "Gatm", "Pdzk1",                               ## proximal tubule
  "Calb1", "Upk3a", "Rprm", "Aqp2", "Trp63", "Crlf1", "Krt18", 
  "Lcn2", "Gata3", "Wfdc2", "Krt19", "Krt8", "Ret", "Mal", "Mia","Wnt11","Wnt9b",              ## ureteric epithelium
  "Plvap", "Cdh5", "Pecam1", "Cldn5", "Esam", "Cd34", "Flt1", "Kdr", "Tie1", "Ecscr",          ## endothelial
  "Cited1", "Six2", "Pclaf", "Eya1", "Uncx", "Spock2", "Crym", "Wnt4",                         ## committing NP
  "Dlk1", "Dcn", "Igf1", "Meg3", "Col1a1", "Postn", "Col3a1", "Lum", "Tbx18",                  ## stroma-ureter+
  "Cited1", "Acta2", "Six2", "Penk", "Col3a1", "Cfh", 
  "Col14a1", "Crym", "Dcn", "Tpm2", "Gucy1a3",                                                 ## nephron progenitor-stromal
  "R3hdml", "Mafb", "Nphs2", "Magi2", "Podxl", "Cdkn1c", "Cldn5", "Nphs1", 
  "Rasl11a", "Dpp4", "Synpo", "Mapt", "Ptpro", "Wt1",                                          ## podocyte
  "Id1","Id2","Id3","Hey1","Plod2","Sfrp2","C1qtnf12", ## C1qtnf12 is the same as Fam132a      
  "Robo2",                                                                                     ## NPCs primed to diff
  "Emx2","Lhx1","Kdm2b","Dapl1","Hnf1b","Krt8",                                                ## differentiation 1
  "Jag1","Pax8","Hes1","Krt8","Calca","Osr2","Myc","Pcp4",                                     ## differentiation 2
  "Top2a","Ube2c","Cenpf","Cdk1","Mki67","Birc5",                                              ## proliferation 1
  "Pclaf","Zwint","Cenpu","Pcna","Cenpk","Tyms","Dut",                                         ## proliferation 2
  "Dpep1", "Chmp1a", "Cisd2", "Manba","Elf4","Elf3"                                            ## GWAS hits
)

regional_marker_genes = c("Cldn2","Spp2","Lrp2","Aqp1","Sptssb","Slc12a1","Slc12a3","Calb1",   ## Nephron
                          "Hsd11b2","Aqp4","Aqp2","Atp6v1g3",                                  ## ureteric epithelium
                          "Kdr","Cdh5",                                                        ## vascular
                          "Cnn1","Dcn",                                                        ## interstitial cells
                          "Calb1","Pgam2","Slc12a3","Pvalb","Mecom","Umod","Egf","Slc12a1",    ## distal tubules
                          "Clcnka","Prox1","Crlf1","Psca","Sptssb","S100a6","Cldn4",           
                          "Apela","Aqp1","Bcl6","Pax2","Pitx2","Slc39a8","Fst","Jag1",
                          "Slc12a2","Gdf15",                                                   ## thin limb of loop of henle
                          "Cyp7b1","Slc7a13","Cyp2e1","Nat8","Slc22a28","Slc22a7","Cyp51",
                          "Slc7a12","Abcc3","Cyp4a14","Slc22a8","Cyp2d26","Spp2","Etv1",
                          "Hnf4a","Lrp2","Slc34a1",                                            ## proximal tubules
                          "Cp","Msmp","Nphs2","Cldn5",                                         ## RC
                          "Fxyd3","Krt5","Upk1b","Gsmc2"                                      ## deep medullary epithelium
                          )



marker_gene_science_paper_intersetc = marker_gene_science_paper[which(marker_gene_science_paper %in% genes.gr$name)]

develop_marker_intersect = marker_gene_development_paper[which(marker_gene_development_paper %in% genes.gr$name)]
union_marker_genes = union(marker_gene_science_paper_intersetc,marker.genes1)
union_marker_genes = union(union_marker_genes, develop_marker_intersect)
union_marker_genes = union(union_marker_genes, regional_marker_genes)

genes.sel.gr <- genes.gr[which(genes.gr$name %in% union_marker_genes)];


# re-add the cell-by-bin matrix to the snap object;
x.after.sp = addBmatToSnap(x.after.sp,num.cores=10);
x.after.sp = createGmatFromMat(
  obj=x.after.sp, 
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=10
);


# normalize the cell-by-gene matrix
x.after.sp = scaleCountMatrix(
  obj=x.after.sp, 
  cov=Matrix::rowSums(x.after.sp@bmat),
  mat="gmat",
  method = "RPM"
);

pdf("Nrbp1 marker gene without MAGIC log transformed in snapATAC results 20 PCs.pdf")


  plotFeatureSingle(
    obj=x.after.sp,
    feature.value=log(x.after.sp@gmat[,"Nrbp1"]+1),
    method="tsne", 
    main="Nrbp1",
    point.size=0.1, 
    point.shape=19, 
    # down.sample=10000,
    quantiles=c(0, 1)
  )
dev.off()

# smooth the cell-by-gene matrix
x.after.sp = runMagic(
  obj=x.after.sp,
  input.mat="gmat",
  step.size=3
);



# saveRDS(object = x.after.sp, file = "all_comb snapATAC after batch_correction with regional genes.RDS")
saveRDS(object = x.after.sp, file = "all_comb snapATAC after batch_correction with all genes without MAGIC.RDS")

## hierarchical clustering results
ensemble.ls = lapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), function(x){
  SnapATAC::colMeans(x.after.sp[x,], mat="bmat");
})
# cluster using 1-cor as distance  
pdf("All_comb snapATAC heirarchical clustering.pdf")
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
plotViz(
  obj=x.after.sp,
  method="tsne", 
  main="All_comb Cluster",
  point.color=x.after.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  # down.sample=10000,
  legend.add=FALSE
);
plot(hc, hang=-1, xlab="");

dev.off()


pdf("All_comb marker gene without MAGIC log transformed in snapATAC results 20 PCs science paper markers.pdf",height = 35, width = 42)

par(mfrow = c(5,6));
for(i in 1:30){
  plotFeatureSingle(
    obj=x.after.sp,
    feature.value=x.after.sp@gmat[, marker_gene_science_paper_intersetc[i]],
    method="tsne", 
    main=marker_gene_science_paper_intersetc[i],
    point.size=0.1, 
    point.shape=19, 
    # down.sample=10000,
    quantiles=c(0, 1)
  )}
dev.off()



saveRDS(object = x.after.sp, file = "all_comb snapATAC after batch_correction with gene annotation.RDS")



