# TRM-identification
The tool can be used to identify transcriptional regulatory modules (TRMs) enriched within differentially methylated regions (DMRs) by recursively selecting a candidate transcription factor (TF) based on a predefined threshold. Given a set of target DMRs and background regions, it uses HOMER’s known motif enrichment functionality to achieve the goal. The user has the option to cluster the DMS into different subsets. Clustering is based on WGCNA (Weighted Gene Correlation Network Analysis) algorithm. Other clustering algorithms such as K-means can also be used. On the other hand, the user can also directly apply recursive TRM identification on the entire set of DMRs. 

## Required R libraries: 
"dplyr", "WGCNA", "factoextra","pheatmap","igraph","gdata"
