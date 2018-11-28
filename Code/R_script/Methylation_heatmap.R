rm(list=ls())
library(pheatmap)
library(RColorBrewer)
library(WGCNA)

## for heatmap
cols <- colorRampPalette(c('blue', 'yellow'))(100) 

## define TF
TF = 'TET2ko'

###### function for pre-processing and merging replicates
.replicate.merged <- function(x) {
  colnames(x) <- gsub('^X', '', gsub('[.|_]CpG_ML.*', '', colnames(x)))
  colnames(x)[colnames(x) %in% "6wk_F"] <- "6wk"
  colnames(x)[colnames(x) %in% "6wk_neu_F"] <- "6wk NeuN+(F)"
  colnames(x)[colnames(x) %in% "7wk_neu"] <- "7wk NeuN+"
  colnames(x)[colnames(x) %in% "12mo_neu_F"] <- "12mo NeuN+(F)"
  colnames(x)[colnames(x) %in% "6wk_glia_F"] <- "6wk NeuN-(F)"
  colnames(x)[colnames(x) %in% "7wk_glia"] <- "7wk NeuN-"
  colnames(x)[colnames(x) %in% "12mo_glia_F"] <- "12mo NeuN-(F)"
  colnames(x)[colnames(x) %in% "7wk_glia"] <- "7wk NeuN-"
  colnames(x)[colnames(x) %in% "Excitatory.neuron"] <- "Excitatory Neu"
  colnames(x)[colnames(x) %in% "PV.neuron"] <- "PV Neu"
  colnames(x)[colnames(x) %in% "VIP.neuron"] <- "VIP Neu"
  colnames(x)[colnames(x) %in% "Neuron"] <- "E17.5 Neu"
  colnames(x) <-gsub("(m)(hind)|(m)(fore)|(m)(mid)", "\\2\\4\\6", colnames(x))
  
  if(!any(grepl(".rep[1-9]", colnames(x)))) {
    return(x)
  }
  
  y <- x[,!grepl(".rep[1-9]", colnames(x))]
  rep.names <- unique(gsub(".rep[1-9].*", "", colnames(x)[grepl(".rep[1-9]", colnames(x))]))
  
  z <- c()
  for(n in rep.names) {
    tmp <- x[,grepl(n, colnames(x))]
    z <- cbind(z, rowSums(tmp, na.rm=T) / rowSums(!is.na(tmp)))
  }
  colnames(z) <- rep.names
  
  w<-cbind(y, z)
  xlab.sorted = sort(colnames(w),index.return=T)
  xlab <- xlab.sorted$x; w <- w[,xlab.sorted$ix]
  lab.order = c(24:30,23,32:38,31,16:22,15,4,6,7,8,1,5,9,10,2,3,14,40,43,11,12,39,13,41,42) #TET1KO
  #lab.order = c(24:29,23, 31:36, 30, 17:22, 16, 4, 6, 7, 8, 1, 5, 10, 11, 2, 3, 14, 38, 39, 12, 13, 37) # VTCRI
  xlab <- xlab[lab.order]; w <- w[,lab.order]
  return(list(w,xlab))
  
}
 
##### function for plotting
.plot.profile <- function(t,iclust_DMS,xlab,cex=1) {
  ### remove NA for plotting
  x <- t[rowSums(t<0, na.rm=T)<=0 & rowSums(!is.na(t))>=ncol(t),]
  ############### heatmap
  #pdf(paste0('DMS_clust',iclust_DMS,'_heatmap.pdf'), width = 15, height = 13); par(mai=c(3,1,1,1))
  png(paste0('DMS_clust',iclust_DMS,'_heatmap_customized.png'),type="cairo", width = 2000, height = 2000); par(mai=c(3,1,1,1))
  pheatmap(x, color = cols, cluster_col=F, show_rownames=F, clustering_method="ward.D2", 
           fontsize=18, cellwidth=35)
  dev.off()
  ############## methylation profile
  apply(x, 2, mean, na.rm=T) -> avg
  apply(x, 2, sd, na.rm=T) -> sdev
  upper <- approx(1:ncol(x), avg+sdev, n=100, method = "linear")
  lower <- approx(1:ncol(x), avg-sdev, n=100, method = "linear")
  #pdf(paste0('DMS_clust',iclust,'.pdf'),width = 15, height = 13)
  png(paste0('DMS_clust',iclust_DMS,'_profile_customized.png'),type="cairo", width = 1040, height = 580)
  plot(1:ncol(x), avg, type='n', ylim=c(min(lower$y), max(upper$y)), xlab='', ylab='Methylation level', xaxt='n', yaxt='s', cex.axis=cex, cex.lab=cex, cex=cex)
  polygon(c(rev(lower$x), upper$x), c(rev(lower$y), upper$y), col = 'grey80', border = NA)
  ## color points by cluster color
  #points(1:ncol(x), avg, pch=16, col=coltype[iclust_DMS], cex=1.6*cex)
  ## color black for all clusters
  points(1:ncol(x), avg, pch=16, col='black', cex=1.6*cex)
  lines(1:ncol(x), avg, col='black', lty='dashed', lwd=1.8*cex)
  axis(1, at=seq_along(1:ncol(x)), label=xlab, las=2, cex.axis=0.65, cex.lab=cex, cex=cex)
  dev.off()
}


#setwd('/home/bsharmi6/Methyome_data_analysis/TET1KO_TET2KO/TET1KO/DMS_6week_analysis/DMR_up_down_background/Kmeans_WGCNA/Kmeans_cluster_3/')
setwd('/home/bsharmi6/Methyome_data_analysis/TET1KO_TET2KO/TET2KO/DMS_6week_10week_analysis/redo/All_sites/minmodule_500/')
#setwd('/home/bsharmi6/Methyome_data_analysis/TET1KO_TET2KO/EGR1KO/DMS_6week_analysis/')
#setwd('/home/bsharmi6/Methyome_data_analysis/VTCRI/data/DMS_4_list/Math5P6_Math5P23.annotated/')
## load clustering result
load(paste0(getwd(),'/WGCNA_result.RData'))
moduleLabels = net$colors ## membership information of each DMS
mergedColors = labels2colors(net$colors)
unique(mergedColors)

## default color plot
## plot
png(paste0('WGCNA.png'),type="cairo", width = 600, height = 400); par(mai=c(3,1,1,1))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

## change color
mergedColors[mergedColors %in% 'turquoise'] = 'gold'
mergedColors[mergedColors %in% 'yellow'] = 'wheat1'
mergedColors[mergedColors %in% 'blue'] = 'royalblue4'
mergedColors[mergedColors %in% 'green'] = 'aquamarine'
mergedColors[mergedColors %in% 'red'] = 'coral1'
mergedColors[mergedColors %in% 'grey'] = 'red'
mergedColors[mergedColors %in% 'brown'] = 'green'

## customized color plot
png(paste0('WGCNA_customized.png'),type="cairo", width = 600, height = 400); par(mai=c(3,1,1,1))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
### COMMENTED OUT IF YOU WANT TO GENERATE A PDF
#plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
## read meth matrix for assignment for VTCRI
#y <- read.table('VTCRI_meth_matrix.txt', h=T,sep = '\t')
## read meth matrix for assignment for TET1ko/TET2o and different cluster
#y = read.table('/home/bsharmi6/Methyome_data_analysis/TET1KO_TET2KO/TET1KO/DMS_6week_analysis/TET1ko10x_meth_matrix.txt', h=T, sep = '\t')
y <- read.table(paste0(TF, '10x_meth_matrix.txt'), h=T,sep = '\t')
#y <- read.table('cluster_3_mat.txt', h=T,sep = '\t')
## get modules
moduleLabels_list=vector('list', length(unique(moduleLabels)))
names(moduleLabels_list)=unique(moduleLabels)
for(i in 1:length(moduleLabels_list)){
  moduleLabels_list[[i]]=y[moduleLabels %in% names(moduleLabels_list)[i],]
}
## remove unassigned module
#moduleLabels_list[names(moduleLabels_list)=='0'] = NULL
## sort by largest cluster
#moduleLabels_list<-moduleLabels_list[order(sapply(moduleLabels_list,function(x) dim(x)[1]),decreasing=T)]
## use if coloring by clusters
#coltype = c('gold', 'royalblue4', 'green', 'wheat1', 'aquamarine', 'coral1', 'black')
#names(coltype) <-names(moduleLabels_list)
for (iclust in 1:length(moduleLabels_list)){
  #tmp <- y[clara.res$clustering %in% num_clust[iclust],]
  res<-.replicate.merged(moduleLabels_list[[iclust]])
  tmp = res[[1]]
  xlab = res[[2]]
  .plot.profile(tmp,names(moduleLabels_list)[iclust],xlab)
}





