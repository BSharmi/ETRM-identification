# TRM-identification
The tool can be used to identify transcriptional regulatory modules (TRMs) enriched within differentially methylated regions (DMRs) by recursively selecting a candidate transcription factor (TF) based on a predefined threshold. Given a set of target DMRs and background regions, it uses HOMER’s known motif enrichment functionality to achieve the goal. The user has the option to cluster the DMS into different subsets. Clustering is based on WGCNA (Weighted Gene Correlation Network Analysis) algorithm. Other clustering algorithms such as K-means can also be used. On the other hand, the user can also directly apply recursive TRM identification on the entire set of DMRs. 

## Required R libraries: 
"dplyr", "WGCNA", "factoextra", "pheatmap", "igraph", "gdata"

## The package contains the following modules-
1.	Code/Shell_script/shell_script_local: to run the commands on a PC
2.	Code/Shell_script/shell_script_server: for users who wish to run the tool on a server
3.	Code/R_script: R scripts needed by the tool
4.	Examples: sample data and output provided as demo
5.	Motif database: Motif libraries compiled from HOMER and other public databases (for details see the paper)

