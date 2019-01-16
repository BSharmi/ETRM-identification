## TF peak enrichment
rm(list=ls())

library(pheatmap)

## define path
#fpath = '/home/bsharmi6/Methyome_data_analysis/TET1KO_TET2KO/TET1KO/DMS_6week_10week_analysis/redo/All_sites/minmodule_1000/DMS_annotation/ChIPseq/'
fpath = '/home/bsharmi6/Methyome_data_analysis/TET1KO_TET2KO/TET2KO/DMS_6week_10week_analysis/redo/All_sites/minmodule_500/DMS_annotation/ChIPseq/'

## define TF
#TF = 'cluster_1'
#TF = 'Tet2'

## read annotated table. it is annotated with clusters and TF peaks
tmp = read.delim(paste0(fpath, list.files(fpath, pattern = '*annotated.txt')), header = T)

## create a column with cluster annotation
tmp$cluster_annotated = 'Non-TRMs'

## Tet1KO cluster1
tmp$cluster_annotated[! tmp$ap._trm	%in% 0] = 'AP-1-TRM'
tmp$cluster_annotated[! tmp$ato_trm %in% 0] = 'Atoh-TRM'
tmp$cluster_annotated[! tmp$egr_trm %in% 0] = 'Egr1-TRM'
tmp$cluster_annotated[! tmp$mef_trm %in% 0] = 'Mef2-TRM'
tmp$cluster_annotated[! tmp$nf1_trm %in% 0] = 'NF1-TRM'
tmp$cluster_annotated[! tmp$zic_trm %in% 0] = 'Zic3-TRM'
tmp$cluster_annotated[! tmp$sox_trm %in% 0] = 'Sox-TRM'
tmp$cluster_annotated[! tmp$zeb_trm %in% 0] = 'Zeb1-TRM'

## remove no TMRs
tmp = tmp[! tmp$cluster_annotated %in% 'Non-TRMs',]

## order
tmp = tmp[,c(1:3,64,4:55)]

## Tet2KO cluster1
tmp$cluster_annotated[! tmp$atf_trm	%in% 0] = 'Atf3-TRM'
tmp$cluster_annotated[! tmp$mef_trm %in% 0] = 'Mef2-TRM'
tmp$cluster_annotated[! tmp$neu_trm %in% 0] = 'NeuroG2-TRM'
tmp$cluster_annotated[! tmp$nf1_trm %in% 0] = 'NF1-TRM'
tmp$cluster_annotated[! tmp$sox3_trm %in% 0] = 'Sox3-TRM'
tmp$cluster_annotated[! tmp$sox10_trm %in% 0] = 'Sox10-TRM'
tmp$cluster_annotated[! tmp$lhx_trm %in% 0] = 'Lhx-TRM'

## remove no TMRs
tmp = tmp[! tmp$cluster_annotated %in% 'Non-TRMs',]

## order
tmp = tmp[,c(1:3,63,4:55)]

#tmp$cluster_annotated[! tmp$cluster_5_dms %in% 0] = 'cluster_5'
## remove cluster 0
#tmp = tmp[! tmp$cluster_annotated %in% '0',]

#tmp = tmp[,c(1:3,52,4:48)]
# if(TF =='Tet1'){
# 	## order for Tet1
# 	tmp = tmp[,c(1:3,60,4:55)]
# }else{
# 	## order for Tet1
# 	tmp = tmp[,c(1:3,59,4:55)]
# }

## Fishers test
x=tmp
## factor by cluster dms
dat <- c(); total <- c()
for(f in levels(factor(x[,4]))) {
  t <- x[grepl(f, x[,4]),]
  dat <- rbind(dat, colSums(x[grepl(f, x[,4]),5:ncol(x)]))
  total <- rbind(total, sum(grepl(f, x[,4])))
  rm(list=c("t"))
} 
rownames(dat) <- levels(factor(x[,4]))
rownames(total) <- levels(factor(x[,4]))


dat.pval <- c(); dat.est <- c()
for(i in 1:ncol(dat)) for(j in 1:nrow(dat)) {
  tab <- matrix(c(dat[j,i], total[j], sum(dat[,i]), sum(total)), ncol=2) 
  f.et <- fisher.test(tab)
  dat.est <- c(dat.est, f.et$estimate)
  dat.pval <- c(dat.pval, f.et$p.value)
}

dat.est <- matrix(dat.est, ncol=ncol(dat), dimnames=list(rownames(dat), colnames(dat)))
dat.pval <- matrix(dat.pval, ncol=ncol(dat), dimnames=list(rownames(dat), colnames(dat)))
df.dat <- dat.est[,setdiff(colnames(dat.est), c("SRF", "CREB", "Miz1", "Con_enhancers", "Dec_enhancers", "Inc_enhancers", "H3K27ac.KCl","Enhancer"))]

## reorder accroding to cluster
##TET1ko
df.dat = df.dat[c("AP-1-TRM", "Atoh-TRM", "Egr1-TRM", "Mef2-TRM", "NF1-TRM", "Sox-TRM","Zic3-TRM" , "Zeb1-TRM"),]
## TET2KO
df.dat = df.dat[c("Atf3-TRM", "Mef2-TRM" , "NeuroG2-TRM",  "NF1-TRM" , "Sox3-TRM", "Lhx-TRM", "Sox10-TRM"),]
#df.dat <- dat.est[,intersect(colnames(dat.est), c("Ascl1", "Atoh1", "CBP", "Cfos", "CREB", "Egr1", "Fezf2", "Fosb", "Foxa1", "Foxa2", "Foxo3", "Max", "MEF2A", "MEF2C", "NFI", "Npas3", "Npas4", "Oct6", "Olig2", "Zic1_P60.bed", "Zic1_P7.bed", "Sox2.bed", "Sox21.bed", "Sox9.bed", "Tcf3.bed", "Tcf4.bed", "Tbr1.bed"))]

## change cluster names to bring egr1 to front
#if(TF == 'Tet1'){
#rownames(df.dat) = c("cluster_3", "cluster_1", "cluster_2", "cluster_4")
#}else{
#	rownames(df.dat) = c("cluster_2", "cluster_1", "cluster_4", "cluster_3")
#}

## order rows since not clustering on rows
#df.dat = df.dat[order(rownames(df.dat)),]

cols <- colorRampPalette(c('green', 'white', 'red'))(100)
tiff(paste0(fpath, "TF.DMS.corr.tiff"), type="cairo", width = 1000, height = 1000)
#pheatmap(cor(df.dat), color = cols, show_rownames=T, clustering_method="average", fontsize=18)
pheatmap(cor(df.dat), color = cols, show_rownames=T, clustering_method="ward.D2", fontsize=18)
dev.off()

## neuron not needed if we dont see clusters showing difference
#excitatory.neurons <- c("cluster_1")
#inhibitory.neurons <- c("cluster_2")
#neuron <- rep(NA, nrow(dat.pval))
#neuron[rownames(dat.est) %in% excitatory.neurons] <- "Excitatory"
#neuron[rownames(dat.est) %in% inhibitory.neurons] <- "Inhibitory"

#anno <- data.frame(Neuron = factor(neuron))
#anno_color <- list(Neuron = c("orange1", "cyan"))
#names(anno_color$Neuron) <- levels(anno$Neuron)
#rownames(anno) <- rownames(dat.est)

rate <- (max(df.dat)-1)/(max(df.dat)-min(df.dat))
cols <- c(colorRampPalette(c('steelblue', 'white'))(100-round(rate*100)), colorRampPalette(c('white', 'red'))(round(rate*100)))
## tiff
tiff(paste0(fpath, "TF.ctDMRs.FC.tiff"), type="cairo", width = 1200, height = 800)
par(mai=c(1,1,1,1))
#pheatmap(df.dat, color = cols, annotation_row=anno, annotation_colors=anno_color, show_rownames=T, clustering_method="ward.D2", fontsize=18, cellwidth=22, , cellheight=20)
pheatmap(df.dat, color = cols, show_rownames=T, clustering_method="ward.D2", cluster_rows = FALSE, fontsize=18, cellwidth=22, , cellheight=20)
dev.off()

## png
png(paste0(fpath, "TF.ctDMRs.FC.png"), type="cairo", width = 1200, height = 800)
par(mai=c(1,1,1,1))
#pheatmap(df.dat, color = cols, annotation_row=anno, annotation_colors=anno_color, show_rownames=T, clustering_method="ward.D2", fontsize=18, cellwidth=22, , cellheight=20)
pheatmap(df.dat, color = cols, show_rownames=T, clustering_method="ward.D2", cluster_rows = FALSE, fontsize=18, cellwidth=22, , cellheight=20)
dev.off()




