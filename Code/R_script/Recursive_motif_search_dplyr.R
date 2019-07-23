#################################### read files and recursively write ################
args <- commandArgs(trailingOnly=TRUE)
## load libraries
library(gdata)
library(dplyr)
library(data.table)
## file path in args[1]
fpath =  args[1] ## path
## top motif in args[2]
topmotif = args[2] 
## read known motif location file
#known.text.file = read.delim(args[3], header=TRUE, sep='\t')
known.text.file = as.data.frame(fread(args[3], select = c(1:3)))
## read all sequence location file
#all.site.file = read.delim(args[4], header=TRUE, sep='\t')
all.site.file = as.data.frame(fread(args[4]))
## output path
#outpath = args[5]
## split first column
tmp = unname(sapply(known.text.file[,1], function(x) strsplit(as.character(x), '\\:|\\-')))
known.chr = sapply(tmp, '[[', 1)
known.start = as.numeric(sapply(tmp, '[[', 2))
known.end = as.numeric(sapply(tmp, '[[', 3))
#known.chr=unname(sapply(known.text.file[,1], function(x) strsplit(as.character(x), '\\:')[[1]][1]))
#known.start=unname(sapply(known.text.file[,1], function(x) strsplit(as.character(x), '\\:|\\-')[[1]][2]))
#known.end=unname(sapply(known.text.file[,1], function(x) strsplit(as.character(x), '\\:|\\-')[[1]][3]))
known.motif.bed=data.frame(chrom=(known.chr),start=as.numeric(known.start),end = as.numeric(known.end))
## new all site locations
new.all.site.file = anti_join(all.site.file,known.motif.bed,by = c("chrom", "start", "end")) ## with dplyr
## known.text file add the offset column
known.text.file=known.text.file[,c(2:ncol(known.text.file))]
## add chromosome start and end columns
known.text.file=cbind.data.frame(known.motif.bed,known.text.file)
## compute center and add offset
known.text.file$center=known.text.file$start + ceiling((known.text.file$end-known.text.file$start)/2) + known.text.file[,4]
## select unique duplicate rows
known.all.motifs.dat=unique(known.motif.bed[duplicated(known.motif.bed, by=c("chrom","start","end")),])
## get rows setdiff of duplicated
known_unique_sites=anti_join(known.motif.bed,known.all.motifs.dat, by = c("chrom", "start", "end"))
## check if any sites is present after subtracting duplicate candidate TF sites
if(dim(known_unique_sites)[1]>0){
	## get all columns
	temp=inner_join(known_unique_sites,known.text.file,by = c("chrom", "start", "end"))
	## final data frame
	final_dat=data.frame(chrom=as.character(temp$chrom), start = as.numeric(temp$center - 100), end = as.numeric(temp$center + 100))
	## process for duplicated rows
	dup_indx=dim(known.all.motifs.dat)[1]
	for(i_dup in 1:dup_indx){
	  temp = inner_join(known.all.motifs.dat[i_dup,],known.text.file, by = c("chrom", "start", "end"))
	  start_i=min(temp$center)
	  end_i=max(temp$center)
	  temp1=data.frame(chrom=as.character(unique(temp$chrom)),start = as.numeric(start_i-100), end = as.numeric(end_i+100))
	  final_dat=rbind.data.frame(final_dat,temp1)
	}

	## final select unique
	final_dat=unique(final_dat)
	## write candidate motif and subtracted candidate motif file
	write.table(as.data.frame(final_dat), file = paste0(fpath, '/',topmotif,'/',topmotif,'.txt'), sep = "\t", append = F,quote = F,row.names = F,col.names = T,eol = "\n")
	write.table(as.data.frame(new.all.site.file), file = paste0(fpath, '/',topmotif,'_minus/',topmotif,'_minus.txt'), sep = "\t", append = F,quote = F,row.names = F,col.names = T,eol = "\n")
}else{
	## write candidate motif and subtracted candidate motif file
	write.table(as.data.frame(known.all.motifs.dat), file = paste0(fpath, '/',topmotif,'/',topmotif,'.txt'), sep = "\t", append = F,quote = F,row.names = F,col.names = T,eol = "\n")
	write.table(as.data.frame(new.all.site.file), file = paste0(fpath, '/',topmotif,'_minus/',topmotif,'_minus.txt'), sep = "\t", append = F,quote = F,row.names = F,col.names = T,eol = "\n")
}

