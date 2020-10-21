
## compare narrow peak files with gene annotations


## only keep 90028 and 90029 batches
# library("SnapATAC");
# library("Seurat")
# library(viridisLite);
library(GenomicRanges)
library(ggplot2);
library(rtracklayer)



intron = read.table(file = "mm10 gene intron.bed",header = F)
exon = read.table(file = "mm10 gene coding exons.bed",header = F)
# promoter = read.table(file = "mm10 gene upstream 2kb.bed",header = F)
promoter = read.table(file = "mm10 gene upstream 5kb.bed",header = F)
utr5 = read.table(file = "mm10 gene 5' UTR.bed",header = F)
utr3 = read.table(file = "mm10 gene 3' UTR exons.bed",header = F)

names(intron) = c("chr","start","end","ID","score","strand")
names(exon) = c("chr","start","end","ID","score","strand")
names(promoter) = c("chr","start","end","ID","score","strand")
names(utr5) = c("chr","start","end","ID","score","strand")
names(utr3) = c("chr","start","end","ID","score","strand")


intron = makeGRangesFromDataFrame(intron)
exon = makeGRangesFromDataFrame(exon)
promoter = makeGRangesFromDataFrame(promoter)
utr5 = makeGRangesFromDataFrame(utr5)
utr3 = makeGRangesFromDataFrame(utr3)

# promoter_spe = setdiff(promoter,exon,ignore.strand=T)
# promoter_spe = setdiff(promoter_spe,intron,ignore.strand=T)
# promoter_spe = setdiff(promoter_spe,utr5,ignore.strand=T)
# promoter_spe = setdiff(promoter_spe,utr3,ignore.strand=T)


# clu_names = system("ls | grep narrowPeak", intern=TRUE)


p56k27 = read.table("P56_H3K27Ac_narrowPeak_ENCFF338WZP.bed",header = F)
p0k27 = read.table("P0_H3K27Ac_narrowPeak_ENCFF872MVE.bed",header = F)

names(p0k27) = c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
names(p56k27) = c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")

p0k27 = p0k27[p0k27$qValue > 2,]
p56k27 = p56k27[p0k27$qValue > 2,]

quantile(p0k27$qValue,c(0.05,0.25,0.5,0.75,0.975))

p0k27 = makeGRangesFromDataFrame(p0k27)
p56k27 = makeGRangesFromDataFrame(p56k27)


## study the overlap with H3K27ac

a_mat = matrix(nrow = 5, ncol = length(clusters))
rownames(a_mat) = c("P0_H3K27Ac","P56_H3K27Ac","both","none","total")
colnames(a_mat) = clusters


## use the same order as the files listed in the directory
# clusters = c("Endo","NP","NP_LOH","Podo","PT","PT2","PT_out","LOH","DCT","PC_DCT","PC","IC","immune","new_stroma1","new_stroma2")
clusters = c("Endo","NP","Podo","PT","PT2","LOH","PC_DCT","DCT","PC","IC","immune","new_stroma1","new_stroma2")


for(i in clusters){
  t = paste("atac_combine_.",i,"_peaks.narrowPeak",sep = "")
  PT = read.table(t,header = F)
  PT = PT[,c("V1","V2","V3","V4","V6")]
  colnames(PT) <- c("chr","start","end","ID","strand")
  PT = makeGRangesFromDataFrame(PT)
  
  if_PT_P0_H3k27ac = (countOverlaps(query = PT, subject = p0k27,type="any",ignore.strand=T) != 0)
  if_PT_P56_H3k27ac = (countOverlaps(query = PT, subject = p56k27,type="any",ignore.strand=T) != 0)
  
  a_mat["P0_H3K27Ac",i] = sum(if_PT_P0_H3k27ac & !(if_PT_P56_H3k27ac))
  a_mat["P56_H3K27Ac",i] = sum(if_PT_P56_H3k27ac & !(if_PT_P0_H3k27ac))
  a_mat["both",i] = sum(if_PT_P56_H3k27ac & if_PT_P0_H3k27ac)
  a_mat["total",i] = length(if_PT_P56_H3k27ac)
  a_mat["none",i] = a_mat["total",i] - a_mat["both",i] - a_mat["P56_H3K27Ac",i] - a_mat["P0_H3K27Ac",i]
}


write.csv(a_mat, "stratified cell type H3K27Ac overlap.csv",quote = F)



## then, study the distribution for DARs with H3K27Ac
tt = readRDS("cell_type specific open chromatin.rds")

cluster_i = c("PT2","PT","immune","DCT","LOH","stroma1","Endo","IC","stroma2","PC","Podo","NP")     
cluster_j = c("PT2","PT","immune","DCT","LOH","stroma1","Endo","IC","stroma2","PC","Podo","NP")    


a_mat2 = matrix(nrow = 5, ncol = length(cluster_i))
rownames(a_mat2) = c("P0_H3K27Ac","P56_H3K27Ac","both","none","total")
colnames(a_mat2) = cluster_i


for(i in cluster_i){
  PT = tt[[i]]
  PT = PT[,c(1,2,3,6,5)]
  colnames(PT) <- c("chr","start","end","ID","strand")
  PT = makeGRangesFromDataFrame(PT)
  
  if_PT_P0_H3k27ac = (countOverlaps(query = PT, subject = p0k27,type="any",ignore.strand=T) != 0)
  if_PT_P56_H3k27ac = (countOverlaps(query = PT, subject = p56k27,type="any",ignore.strand=T) != 0)
  
  a_mat2["P0_H3K27Ac",i] = sum(if_PT_P0_H3k27ac & !(if_PT_P56_H3k27ac))
  a_mat2["P56_H3K27Ac",i] = sum(if_PT_P56_H3k27ac & !(if_PT_P0_H3k27ac))
  a_mat2["both",i] = sum(if_PT_P56_H3k27ac & if_PT_P0_H3k27ac)
  a_mat2["total",i] = length(if_PT_P56_H3k27ac)
  a_mat2["none",i] = a_mat2["total",i] - a_mat2["both",i] - a_mat2["P56_H3K27Ac",i] - a_mat2["P0_H3K27Ac",i]
}

write.csv(a_mat2, "stratified cell type DAP H3K27Ac overlap.csv",quote = F)



clusters = c("Endo","NP","Podo","PT","PT2","LOH","PC_DCT","CNT","PC","IC","immune","new_stroma1","new_stroma2")

s_mat = matrix(nrow = 7, ncol = length(clusters))
rownames(s_mat) = c("exon","3'-UTR","5'-UTR","intron","promoter","distal","total")
colnames(s_mat) = clusters



for(i in clusters){
  t = paste("atac_combine_.",i,"_peaks.narrowPeak",sep = "")
  PT = read.table(t,header = F)
  PT = PT[,c("V1","V2","V3","V4","V6")]
  colnames(PT) <- c("chr","start","end","ID","strand")
  PT = makeGRangesFromDataFrame(PT)
  
  if_PT_promoter = (countOverlaps(query = PT, subject = promoter,type="any",ignore.strand=T) != 0)
  if_PT_exon = (countOverlaps(query = PT, subject = exon,type="any",ignore.strand=T) != 0)
  if_PT_intron = (countOverlaps(query = PT, subject = intron,type="any",ignore.strand=T) != 0)
  if_PT_utr3 = (countOverlaps(query = PT, subject = utr3,type="any",ignore.strand=T) != 0)
  if_PT_utr5 = (countOverlaps(query = PT, subject = utr5,type="any",ignore.strand=T) != 0)
  
  PT_exon = sum(if_PT_exon)
  PT_utr3 = sum(if_PT_utr3 & (!if_PT_exon))
  PT_utr5 = sum(if_PT_utr5 & (!if_PT_exon))
  PT_intron = sum(if_PT_intron & (!if_PT_utr3) & (!if_PT_utr5))
  PT_promoter = sum(if_PT_promoter & (!if_PT_exon) & (!if_PT_intron) & (!if_PT_utr3) & (if_PT_utr5))
  
  PT_total = length(countOverlaps(query = PT, subject = promoter,type="any",ignore.strand=T))
  PT_distal = PT_total - PT_exon - PT_intron - PT_utr3 - PT_utr5 - PT_promoter
  
  s_mat["exon",i] = PT_exon 
  s_mat["3'-UTR",i] = PT_utr3 
  s_mat["5'-UTR",i] = PT_utr5
  s_mat["intron",i] = PT_intron
  s_mat["promoter",i] = PT_promoter
  s_mat["distal",i] = PT_distal
  s_mat["total",i] = PT_total
}

## then, study the distribution for DARs
tt = readRDS("cell_type specific open chromatin.rds")

cluster_i = c("PT2","PT","immune","DCT","LOH","stroma1","Endo","IC","stroma2","PC","Podo","NP")     
cluster_j = c("PT2","PT","immune","DCT","LOH","stroma1","Endo","IC","stroma2","PC","Podo","NP")    


s_mat2 = matrix(nrow = 7, ncol = length(cluster_i))
rownames(s_mat2) = c("exon","3'-UTR","5'-UTR","intron","promoter","distal","total")
colnames(s_mat2) = cluster_i


for(i in cluster_i){
  PT = tt[[i]]
  PT = PT[,c(1,2,3,6,5)]
  colnames(PT) <- c("chr","start","end","ID","strand")
  PT = makeGRangesFromDataFrame(PT)
  
  if_PT_promoter = (countOverlaps(query = PT, subject = promoter,type="any",ignore.strand=T) != 0)
  if_PT_exon = (countOverlaps(query = PT, subject = exon,type="any",ignore.strand=T) != 0)
  if_PT_intron = (countOverlaps(query = PT, subject = intron,type="any",ignore.strand=T) != 0)
  if_PT_utr3 = (countOverlaps(query = PT, subject = utr3,type="any",ignore.strand=T) != 0)
  if_PT_utr5 = (countOverlaps(query = PT, subject = utr5,type="any",ignore.strand=T) != 0)
  
  PT_exon = sum(if_PT_exon)
  PT_utr3 = sum(if_PT_utr3 & (!if_PT_exon))
  PT_utr5 = sum(if_PT_utr5 & (!if_PT_exon))
  PT_intron = sum(if_PT_intron & (!if_PT_utr3) & (!if_PT_utr5))
  PT_promoter = sum(if_PT_promoter & (!if_PT_exon) & (!if_PT_intron) & (!if_PT_utr3) & (if_PT_utr5))
  
  PT_total = length(countOverlaps(query = PT, subject = promoter,type="any",ignore.strand=T))
  PT_distal = PT_total - PT_exon - PT_intron - PT_utr3 - PT_utr5 - PT_promoter
  
  s_mat2["exon",i] = PT_exon 
  s_mat2["3'-UTR",i] = PT_utr3 
  s_mat2["5'-UTR",i] = PT_utr5
  s_mat2["intron",i] = PT_intron
  s_mat2["promoter",i] = PT_promoter
  s_mat2["distal",i] = PT_distal
  s_mat2["total",i] = PT_total
}

write.csv(s_mat2, "stratified cell type specific elements annotations.csv",quote = F)



## then, study the number of combined peaks

peaks.names = system("ls | grep _peaks.narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
PT = reduce(Reduce(c, peak.gr.ls));


if_PT_promoter = (countOverlaps(query = PT, subject = promoter,
                                  type="any",
                                  ignore.strand=T) != 0)
if_PT_exon = (countOverlaps(query = PT, subject = exon,
                                  type="any",
                                  ignore.strand=T) != 0)
if_PT_intron = (countOverlaps(query = PT, subject = intron,
                                  type="any",
                                  ignore.strand=T) != 0)
if_PT_utr3 = (countOverlaps(query = PT, subject = utr3,
                                  type="any",
                                  ignore.strand=T) != 0)
if_PT_utr5 = (countOverlaps(query = PT, subject = utr5,
                                  type="any",
                                  ignore.strand=T) != 0)


PT_exon = sum(if_PT_exon)
PT_utr3 = sum(if_PT_utr3 & (!if_PT_exon))
PT_utr5 = sum(if_PT_utr5 & (!if_PT_exon))
PT_intron = sum(if_PT_intron & (!if_PT_utr3) & (!if_PT_utr5))
PT_promoter = sum(if_PT_promoter & (!if_PT_exon) & (!if_PT_intron) & (!if_PT_utr3) & (if_PT_utr5))

PT_total = length(countOverlaps(query = PT, subject = promoter,
                                type="any",
                                ignore.strand=T))


PT_distal = PT_total - PT_exon - PT_intron - PT_utr3 - PT_utr5 - PT_promoter

tot_stra = c(PT_exon,PT_utr3,PT_utr5,PT_intron,PT_promoter,PT_distal,PT_total)
s_mat2 = cbind(s_mat,tot_stra)


  
  bulk_f = c("P01","P02","P21-1","P21-2","P56-1","P56-2")
  
  b_mat = matrix(nrow = 7, ncol = 3)
  rownames(b_mat) = c("exon","3'-UTR","5'-UTR","intron","promoter","distal","total")
  colnames(b_mat) = c("P01","P21","P56")

  PT_list = rep(list(), length = 3)
  for(i in 2:3){
    m = bulk_f[2*i]
    mm = bulk_f[2*i -1]
    t = paste(m,"-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak",sep = "")
    tt = paste(mm,"-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak",sep = "")
    
    PTf = read.table(t,header = F)
    colnames(PTf) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
    PTf = PTf[PTf$qValue >= 2,]
    PT = makeGRangesFromDataFrame(PTf)
    
    PTf2 = read.table(tt,header = F)
    colnames(PTf2) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
    PTf2 = PTf2[PTf2$qValue >= 2,]
    PT2 = makeGRangesFromDataFrame(PTf2)
    
    o1 = subsetByOverlaps(PT, PT2)
    o2 = subsetByOverlaps(PT2, PT)
    o3 = c(o1,o2)
    o3 = reduce(o3)
    PT_list[[i]] = o3
  }
  
for(i in 1:3){
  # t = paste(i,"-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak",sep = "")
  # PTf = read.table(t,header = F)
  # colnames(PTf) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
  # PTf = PTf[PTf$qValue >=2,]
  # PT = makeGRangesFromDataFrame(PTf)
  
  PT = PT_list[[i]]
  
  if_PT_promoter = (countOverlaps(query = PT, subject = promoter,type="any",ignore.strand=T) != 0)
  if_PT_exon = (countOverlaps(query = PT, subject = exon,type="any",ignore.strand=T) != 0)
  if_PT_intron = (countOverlaps(query = PT, subject = intron,type="any",ignore.strand=T) != 0)
  if_PT_utr3 = (countOverlaps(query = PT, subject = utr3,type="any",ignore.strand=T) != 0)
  if_PT_utr5 = (countOverlaps(query = PT, subject = utr5,type="any",ignore.strand=T) != 0)
  
  PT_exon = sum(if_PT_exon)
  PT_utr3 = sum(if_PT_utr3 & (!if_PT_exon))
  PT_utr5 = sum(if_PT_utr5 & (!if_PT_exon))
  PT_intron = sum(if_PT_intron & (!if_PT_utr3) & (!if_PT_utr5))
  PT_promoter = sum(if_PT_promoter & (!if_PT_exon) & (!if_PT_intron) & (!if_PT_utr3) & (if_PT_utr5))
  
  PT_total = length(countOverlaps(query = PT, subject = promoter,type="any",ignore.strand=T))
  PT_distal = PT_total - PT_exon - PT_intron - PT_utr3 - PT_utr5 - PT_promoter
  
  b_mat["exon",i] = PT_exon 
  b_mat["3'-UTR",i] = PT_utr3 
  b_mat["5'-UTR",i] = PT_utr5
  b_mat["intron",i] = PT_intron
  b_mat["promoter",i] = PT_promoter
  b_mat["distal",i] = PT_distal
  b_mat["total",i] = PT_total
}

f_mat = cbind(s_mat2,b_mat)
f_mat = cbind(tot_stra,b_mat)

write.csv(f_mat,file = "stratified peaks information sc with bulk 5 kb as promoters merged two bulk data.csv",quote = F)

stratify_matrix = matrix(nrow = 6, ncol = 12)
colnames(stratify_matrix) = ID_scs
rownames(stratify_matrix) = ID_bulks

for(i in ID_bulks){
  for(j in ID_scs){
    overlap_count = sum(countOverlaps(query = bulk_gr_narrowPeaks[[i]], subject = gr_narrowPeaks[[j]],
                                      maxgap=-1L, minoverlap=0L,
                                      type="any",
                                      ignore.strand=T) != 0)
    overlap_matrix[i,j] = overlap_count
  }
}



hmmp0 = read.table(file = "P0_kidney_HMM_mm10.bed",stringsAsFactors = F)
hmmp0 = hmmp0[,c("V1","V2","V3","V4")]

# [1] "Heterochromatin" "Weak_Enhancer"   "Transcript"      "Weak_Promoter"  
# [5] "Active_Promoter" "Strong_Enhancer" "Repressed"      
hmmp0$V4[hmmp0$V4 == "Weak_Enhancer"] = "Enhancer"
hmmp0$V4[hmmp0$V4 == "Strong_Enhancer"] = "Enhancer"

hmmp0$V4[hmmp0$V4 == "Weak_Promoter"] = "Promoter"
hmmp0$V4[hmmp0$V4 == "Active_Promoter"] = "Promoter"
hmmp0$V4[hmmp0$V4 == "Heterochromatin"] = "Other"
hmmp0$V4[hmmp0$V4 == "Repressed"] = "Other"

hmm0en = hmmp0[hmmp0$V4 == "Enhancer",]
hmm0pr = hmmp0[hmmp0$V4 == "Promoter",]
hmm0tr = hmmp0[hmmp0$V4 == "Transcript",]
hmm0ot = hmmp0[hmmp0$V4 == "Other",]


hmm0en = hmm0en[,c("V1","V2","V3")]
hmm0pr = hmm0pr[,c("V1","V2","V3")]
hmm0tr = hmm0tr[,c("V1","V2","V3")]
hmm0ot = hmm0ot[,c("V1","V2","V3")]

colnames(hmm0en) = c("chr","start","end")
colnames(hmm0pr) = c("chr","start","end")
colnames(hmm0tr) = c("chr","start","end")
colnames(hmm0ot) = c("chr","start","end")


hmm0en = makeGRangesFromDataFrame(hmm0en)
hmm0pr = makeGRangesFromDataFrame(hmm0pr)
hmm0tr = makeGRangesFromDataFrame(hmm0tr)
hmm0ot = makeGRangesFromDataFrame(hmm0ot)



hmmp8 = read.table(file = "W8_kidney_HMM_mm10.bed",stringsAsFactors = F)
hmmp8 = hmmp8[,c("V1","V2","V3","V4")]

# [1] "Heterochromatin" "Strong_Enhancer" "Transcript"      "Active_Promoter"
# [5] "Insulator"       "Repressed"       "Weak_Promoter"   "Weak_Enhancer"   


hmmp8$V4[hmmp8$V4 == "Weak_Enhancer"] = "Enhancer"
hmmp8$V4[hmmp8$V4 == "Strong_Enhancer"] = "Enhancer"

hmmp8$V4[hmmp8$V4 == "Weak_Promoter"] = "Promoter"
hmmp8$V4[hmmp8$V4 == "Active_Promoter"] = "Promoter"
hmmp8$V4[hmmp8$V4 == "Heterochromatin"] = "Other"
hmmp8$V4[hmmp8$V4 == "Repressed"] = "Other"
hmmp8$V4[hmmp8$V4 == "Insulator"] = "Other"

hmm8en = hmmp8[hmmp8$V4 == "Enhancer",]
hmm8pr = hmmp8[hmmp8$V4 == "Promoter",]
hmm8tr = hmmp8[hmmp8$V4 == "Transcript",]
hmm8ot = hmmp8[hmmp8$V4 == "Other",]


hmm8en = hmm8en[,c("V1","V2","V3")]
hmm8pr = hmm8pr[,c("V1","V2","V3")]
hmm8tr = hmm8tr[,c("V1","V2","V3")]
hmm8ot = hmm8ot[,c("V1","V2","V3")]

colnames(hmm8en) = c("chr","start","end")
colnames(hmm8pr) = c("chr","start","end")
colnames(hmm8tr) = c("chr","start","end")
colnames(hmm8ot) = c("chr","start","end")


hmm8en = makeGRangesFromDataFrame(hmm8en)
hmm8pr = makeGRangesFromDataFrame(hmm8pr)
hmm8tr = makeGRangesFromDataFrame(hmm8tr)
hmm8ot = makeGRangesFromDataFrame(hmm8ot)

en0_en8 = sum(countOverlaps(query = hmm0en, subject = hmm8en, type = "any", ignore.strand = T) != 0) # 47013
en0_en_not8 = sum(countOverlaps(query = hmm0en, subject = hmm8en, type = "any", ignore.strand = T) == 0) # 42191

pr0_pr8 = sum(countOverlaps(query = hmm0pr, subject = hmm8pr, type = "any", ignore.strand = T) != 0) # 38258
pr0_pr_not8 = sum(countOverlaps(query = hmm0pr, subject = hmm8pr, type = "any", ignore.strand = T) == 0) # 10755

en0_pr8 = sum((countOverlaps(query = hmm0en, subject = hmm8pr, type = "any", ignore.strand = T) != 0) &
                (countOverlaps(query = hmm0en, subject = hmm8en, type = "any", ignore.strand = T) == 0)) # 9100
pr0_en8 = sum((countOverlaps(query = hmm0pr, subject = hmm8en, type = "any", ignore.strand = T) != 0) &
                (countOverlaps(query = hmm0pr, subject = hmm8pr, type = "any", ignore.strand = T) == 0)) # 6958







for(i in clusters){
  t = paste("atac_combine_.",i,"_peaks.narrowPeak",sep = "")
  PT = read.table(t,header = F)
  PT = PT[,c("V1","V2","V3","V4","V6")]
  colnames(PT) <- c("chr","start","end","ID","strand")
  PT = makeGRangesFromDataFrame(PT)
  
  if_PT_promoter = (countOverlaps(query = PT, subject = hmm0pr,type="any",ignore.strand=T) != 0)
  if_PT_promoter2 = (countOverlaps(query = PT, subject = hmm8pr,type="any",ignore.strand=T) != 0)
  
  if_PT_enhancer = (countOverlaps(query = PT, subject = hmm0en,type="any",ignore.strand=T) != 0)
  if_PT_enhancer2 = (countOverlaps(query = PT, subject = hmm8en,type="any",ignore.strand=T) != 0)
  
  if_PT_transc = (countOverlaps(query = PT, subject = hmm0tr,type="any",ignore.strand=T) != 0)
  if_PT_transc2 = (countOverlaps(query = PT, subject = hmm8pr,type="any",ignore.strand=T) != 0)
  
  if_PT_other = (countOverlaps(query = PT, subject = hmm0ot,type="any",ignore.strand=T) != 0)
  if_PT_other2 = (countOverlaps(query = PT, subject = hmm8ot,type="any",ignore.strand=T) != 0)
  
  PT_promoter = sum(if_PT_promoter2 & if_PT_promoter)
  PT_promoter0_enhancer8 = sum(if_PT_promoter & if_PT_enhancer2 & (!if_PT_promoter2))
  
  
  PT_utr3 = sum(if_PT_utr3 & (!if_PT_exon))
  PT_utr5 = sum(if_PT_utr5 & (!if_PT_exon))
  PT_intron = sum(if_PT_intron & (!if_PT_utr3) & (!if_PT_utr5))
  PT_promoter = sum(if_PT_promoter & (!if_PT_exon) & (!if_PT_intron) & (!if_PT_utr3) & (if_PT_utr5))
  
  PT_total = length(countOverlaps(query = PT, subject = promoter,type="any",ignore.strand=T))
  PT_distal = PT_total - PT_exon - PT_intron - PT_utr3 - PT_utr5 - PT_promoter
  
  s_mat["exon",i] = PT_exon 
  s_mat["3'-UTR",i] = PT_utr3 
  s_mat["5'-UTR",i] = PT_utr5
  s_mat["intron",i] = PT_intron
  s_mat["promoter",i] = PT_promoter
  s_mat["distal",i] = PT_distal
  s_mat["total",i] = PT_total
}




## study the overlap between bulk and single cell data


peaks.names = system("ls | grep _peaks.narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
comSc = reduce(Reduce(c, peak.gr.ls));




comSc2 = keepStandardChromosomes(comSc,pruning.mode="coarse") 
comSc2 = dropSeqlevels(comSc2, "chrM",pruning.mode="coarse")

comSc3 = dropSeqlevels(comSc, "chrM",pruning.mode="coarse")


peaks.names_bulk = system("ls | grep 300K.narrowPeak", intern=TRUE);
peak.gr.ls_bulk = lapply(peaks.names_bulk, function(x){
  peak.df = read.table(x)
  colnames(peak.df) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
  peak.df = peak.df[peak.df$qValue >=2,]
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
bulk_red = reduce(Reduce(c, peak.gr.ls_bulk));

bulk_red2 = keepStandardChromosomes(bulk_red,pruning.mode="coarse") 












overlap_count_Sc_bulk = sum(countOverlaps(query = comSc2, subject = bulk_red2,
                                               maxgap=-1L, minoverlap=0L,
                                               type="any",
                                               ignore.strand=T) != 0)
total_count_Sc_bulk = length(countOverlaps(query = comSc2, subject = bulk_red2,
                                                maxgap=-1L, minoverlap=0L,
                                                type="any",
                                                ignore.strand=T) != 0)


overlap_count_bulk_Sc = sum(countOverlaps(query = bulk_red2 , subject = comSc2,
                                               maxgap=-1L, minoverlap=0L,
                                               type="any",
                                               ignore.strand=T) != 0)
total_count_bulk_Sc = length(countOverlaps(query = bulk_red2, subject = comSc2,
                                                maxgap=-1L, minoverlap=0L,
                                                type="any",
                                                ignore.strand=T) != 0)


overlap_count_Sc_bulk = sum(countOverlaps(query = comSc3, subject = bulk_red,
                                               maxgap=-1L, minoverlap=0L,
                                               type="any",
                                               ignore.strand=T) != 0)
total_count_Sc_bulk = length(countOverlaps(query = comSc3, subject = bulk_red,
                                                maxgap=-1L, minoverlap=0L,
                                                type="any",
                                                ignore.strand=T) != 0)


overlap_count_bulk_Sc = sum(countOverlaps(query = bulk_red , subject = comSc3,
                                               type="any",
                                               ignore.strand=T) != 0)
total_count_bulk_Sc = length(countOverlaps(query = bulk_red, subject = comSc3,
                                                maxgap=-1L, minoverlap=0L,
                                                type="any",
                                                ignore.strand=T) != 0)





t = "P01-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak"
t2 = "P02-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak"

PTf = read.table(t,header = F)
PTf2 = read.table(t2,header = F)

colnames(PTf) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
colnames(PTf2) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
PTf = PTf[PTf$qValue >=2,]
PTf2 = PTf2[PTf2$qValue >=2,]

bulkp0 = makeGRangesFromDataFrame(PTf)
bulkp02 = makeGRangesFromDataFrame(PTf2)

bulkp0int = GenomicRanges::intersect(bulkp0,bulkp02,ignore.strand = TRUE)

t = "P21-1-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak"
t2 = "P21-2-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak"

PTf = read.table(t,header = F)
PTf2 = read.table(t2,header = F)

colnames(PTf) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
colnames(PTf2) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
PTf = PTf[PTf$qValue >=2,]
PTf2 = PTf2[PTf2$qValue >=2,]

bulkp21 = makeGRangesFromDataFrame(PTf)
bulkp212 = makeGRangesFromDataFrame(PTf2)

bulkp21int = GenomicRanges::intersect(bulkp21,bulkp212,ignore.strand = TRUE)

t = "P56-1-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak"
t2 = "P56-2-50_R1_001.trim.merged.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak"

PTf = read.table(t,header = F)
PTf2 = read.table(t2,header = F)

colnames(PTf) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
colnames(PTf2) <- c("chr","start","end","ID","score","strand","signalValue","pValue","qValue","peak")
PTf = PTf[PTf$qValue >=2,]
PTf2 = PTf2[PTf2$qValue >=2,]

bulkp56 = makeGRangesFromDataFrame(PTf)
bulkp562 = makeGRangesFromDataFrame(PTf2)

bulkp56int = GenomicRanges::intersect(bulkp56,bulkp562,ignore.strand = TRUE)

bulkall_reduce = c(bulkp0int,bulkp21int,bulkp56int)
bulkall_reduce = GenomicRanges::reduce(bulkall_reduce,drop.empty.ranges = T,min.gapwidth = 0,with.revmap = T)

bulkall_reduce2 = keepStandardChromosomes(bulkall_reduce,pruning.mode="coarse") 




overlap_count_comSc_bulkall = sum(countOverlaps(query = comSc2, subject = bulkall_reduce2,
                                  maxgap=-1L, minoverlap=0L,
                                  type="any",
                                  ignore.strand=T) != 0)
total_count_comSc_bulkall = length(countOverlaps(query = comSc2, subject = bulkall_reduce2,
                                             maxgap=-1L, minoverlap=0L,
                                             type="any",
                                             ignore.strand=T) != 0)


overlap_count_bulkall_comSc = sum(countOverlaps(query = bulkall_reduce2 , subject = comSc2,
                                  maxgap=-1L, minoverlap=0L,
                                  type="any",
                                  ignore.strand=T) != 0)
total_count_bulkall_comSc = length(countOverlaps(query = bulkall_reduce2, subject = comSc2,
                                             maxgap=-1L, minoverlap=0L,
                                             type="any",
                                             ignore.strand=T) != 0)

mmm = matrix()


atac_combine_.PT_out_peaks.narrowPeak




## study intersect
grgene <- GRanges(
  seqnames = c('chr19', 'chr8', 'chr20'),
  ranges = IRanges(start = c(100, 1000, 100),
                   end = c(200, 2000, 200)),
  strand = c('-', '+', '-')
) # ignoring the GENEID column
grgene2 <- GRanges(
  seqnames = c('chr19', 'chr8', 'chr20'),
  ranges = IRanges(start = c(200, 1001, 100),
                   end = c(300, 2000, 200)),
  strand = c('-', '+', '-')
) # ignoring the GENEID column

tryint = GenomicRanges::intersect(grgene,grgene2, ignore.strand = TRUE)



ckd_file = read.table("Combined_FullListSNPs_Uniq_Liftover2mm10.bed",header = F)

names(ckd_file) = c("chr","start","end","SNP_ID","strand")
ckd_file = ckd_file[,1:4]

ckd = makeGRangesFromDataFrame(ckd_file)
