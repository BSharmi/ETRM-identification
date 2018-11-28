#################################### read files and recursively plot ################
rm(list=ls())
args <- commandArgs(trailingOnly=TRUE)
## read all known results file
all_known_dat = read.delim(args[1], header=TRUE, sep='\t')
## graph type 
graph_type = args[2]
## outpath
outpath = args[3]
cat (outpath)
## TF name
tfname = args [4]
## add library
library(igraph)	
## preprocess TF name to remove - from last
if(grepl('\\-$', tfname)) tfname = gsub("-", "", tfname) else tfname = tfname
## all known results
if(graph_type == "all"){
	## select top significant motifs
	top_known_all = unlist(lapply(as.character(all_known_dat$Motif.Name[all_known_dat$P.value< 1e-10]), function(x) strsplit(x, split = '\\(')[[1]][1]))
	## check if any exists
	if(!is.null(top_known_all)){
		sigpval_known_all= all_known_dat$P.value[all_known_dat$P.value< 1e-10]
		top_sites = all_known_dat[all_known_dat$P.value< 1e-10,6]
		## plot
		g=graph.empty(length(top_known_all), directed=TRUE)
		g=set.vertex.attribute(g, "name", value=paste(top_known_all,sigpval_known_all,top_sites,sep = '\n'))
		V(g)$label.cex=0.5
		pdf(paste0(outpath,"/knownResults_",tfname,".pdf"), 7, 5)
		plot(g,vertex.size=20,vertex.color="green")
		dev.off()
	}	
}else{
	## read file
	top_known_all = unlist(lapply(as.character(all_known_dat$Motif.Name[all_known_dat$P.value< 1e-10]), function(x) strsplit(x, split = '\\(')[[1]][1]))
	## check if any exists
	if(!is.null(top_known_all)){
		sigpval_known_all= all_known_dat$P.value[all_known_dat$P.value< 1e-10]
		top_sites = all_known_dat[all_known_dat$P.value< 1e-10,6]
		## for fully connected graph
		#g <- make_full_graph(length(top_known_all), directed = F)
		## for star graph
		## get center
		gcenter = which(tolower(top_known_all) %in% tfname)
		## none name matches
		if(length(gcenter)==0) gcenter=1
		g=graph.empty(length(top_known_all), directed=TRUE)
		#g <-make_star(length(top_known_all), mode = "undirected",center=gcenter)
		#g$layout <- layout_in_circle
		g=set.vertex.attribute(g, "name", value=paste(top_known_all,sigpval_known_all,top_sites,sep = '\n'))
		V(g)$label.cex=0.5
		pdf(paste0(outpath,"/knownResults_",tfname,".pdf"), 10, 8)
		plot(g,vertex.size=20,vertex.color="green")
		dev.off()
	}	
}



