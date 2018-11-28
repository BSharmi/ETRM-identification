# TRM-identification
The tool can be used to identify transcriptional regulatory modules (TRMs) enriched within differentially methylated regions (DMRs) by recursively selecting a candidate transcription factor (TF) based on a predefined threshold. Given a set of target DMRs and background regions, it uses HOMER’s known motif enrichment functionality to achieve the goal. The user has the option to cluster the DMS into different subsets. Clustering is based on WGCNA (Weighted Gene Correlation Network Analysis) algorithm. Other clustering algorithms such as K-means can also be used. On the other hand, the user can also directly apply recursive TRM identification on the entire set of DMRs. 

## Required R libraries: 
"dplyr", "WGCNA", "pheatmap", "igraph", "gdata"

## The package contains the following modules-
1.	Code/Shell_script/shell_script_local: to run the commands on a PC
2.	Code/Shell_script/shell_script_server: for users who wish to run the tool on a server
3.	Code/R_script: R scripts needed by the tool
4.	Examples: sample data and output provided as demo
5.	Motif database: Motif libraries compiled from HOMER and other public databases (for details see the paper)

Please download a genome fa file (e.g. mm10.fa). It is a huge file and has not been provided here.

## Below are some command line arguments to perform different functions:

### Cluster DMS 

./run_DMS_clustering_script.sh –o <path to output directory> -d <3 column BED format DMS file> -x <Methylome matrix> -db <3 column BED format DMS background file> -R <path to R script files>

Example on demo data:
./run_DMS_clustering_script.sh -o Output/ -d DMS.txt -x DMS_methylome.txt -db DMS_background.txt –R R_script/

Explanation of command:
-o: path to store the output
-d: tab delimited three column BED format DMS file with three columns: chromosome, start, end
-x: methylome matrix. One column indicates methylome values for DMS location, one column indicates methylome values for one sample. Number of rows of x must match with the number of rows in the DMS file
-db: tab delimited three column BED format DMS background file to create background regions for each DMS cluster
-R: path to R script files
-h: explains required parameters
The command will cluster DMS, and create a folder for each cluster name containing the tab delimited DMS files and the corresponding background files

### Recursive TRM identification on entire DMS:

./recursive_motif_identification_noclustering.sh –i <path to input directory containing DMS> -H <path to motif database> -r <fasta file> -s <sequence extractor script> -R <path to R script files>

Example on demo data:
./recursive_motif_identification_noclustering.sh -i Data/ -H Motif_db –r mm10.fa –s sequence_extractor.pl –R R_script/

Explanation of command:
-i: path containing the DMS file
-H: path containing motif database 
-s:  sequence extractor Perl script
-R: path to R script files
-r: reference fasta file; examples are mm10.fa, hg19.fa etc.
-h: explains required parameters


### Recursive TRM identification on DMS cluster :

./recursive_motif_identification_clustering.sh –i <path to input directory containing DMS clusters> -H <path to motif database> -r <fasta file> -s <sequence extractor script> -R <path to R script files>

Example on demo data:
./recursive_motif_identification_clustering.sh -i Data/ -H Motif_db –r mm10.fa –s sequence_extractor.pl –R R_script/

Explanation of command:
-i: path containing the DMS clusters
-H: path containing motif database 
-s:  sequence extractor Perl script
-R: path to R script files
-r: reference fasta file; examples are mm10.fa, hg19.fa etc
-h: explains required parameters


