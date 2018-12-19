# Recursive Transcriptional Regulatory Module (TRM) identification using Transcription Factor motifs
The tool can be used to identify transcriptional regulatory modules (TRMs) enriched within differentially methylated regions (DMRs) by recursively selecting a candidate transcription factor (TF) based on a predefined threshold. Given a set of target DMRs and background regions, it uses HOMER’s known motif enrichment functionality to achieve the goal. The user has the option to cluster the DMRs into different subsets. Clustering is based on WGCNA (Weighted Gene Correlation Network Analysis) algorithm. Other clustering algorithms such as K-means can also be used. On the other hand, the user can also directly apply recursive TRM identification on the entire set of DMRs. 

## Methylated motifs
In order to add methylated motif PWMs to an existing motif databse such as provided by HOMER, please add the motif PWMs to an existing library. For example, to add methylated motifs to HOMER's known motif library of vertebrates, please copy the motifs (motifs.zip) to HOMER library/data/knownTFs/motifs and unzip the file. The user can also append the motifs to an existing motif file.

## Required R libraries: 
"dplyr", "WGCNA", "pheatmap", "igraph", "gdata"

## The package contains the following modules-
1.	Code/Shell_script/shell_script_local: to run the commands on a local system
2.	Code/Shell_script/shell_script_server: for users who wish to run the tool on a server
3.	Code/R_script: R scripts needed by the tool
4.	Examples: sample data and output provided as demo
5.	Motif database: Motif libraries compiled from HOMER and other public databases (for details see the paper)

Please download a genome fa file (e.g. mm10.fa). It is a large file and has not been provided here.

## Below are some command line arguments to perform different functions:

### Cluster DMR 

./run_DMR_clustering_script.sh -R <path to R script files> –o <path to output directory> -d <3 column BED format DMS file> -x <Methylome matrix> -db <3 column BED format DMS background file> 

Example on demo data:

```
cd Code/Shell_script/shell_script_local/ 

./run_DMR_clustering_script.sh –R Code/R_script/DMS_clustering_general.R -o Examples/Output/ -d Examples/DMS.txt -x Examples/Methylation_matrix.txt -db Examples/DMS_background.txt 
```

Explanation of command:
-o: path to store the output
-d: tab delimited three column BED format DMS file with three columns: chromosome, start, end
-x: methylome matrix. One column indicates methylome values for DMS location, one column indicates methylome values for one sample. Number of rows of x must match with the number of rows in the DMS file
-db: tab delimited three column BED format DMS background file to create background regions for each DMS cluster
-R: R script used for clustering
-h: explains required parameters
The command will cluster DMRs, and create a folder for each cluster name containing the tab delimited DMRs files and the corresponding background files

#### Running scripts on server

Since WGCNA can be slow on large DMR matrices, it might be efficient to run the scripts by submitting as jobs on a high performance computing system. The script below shows how to run the scripts on server -  

```
sbatch --export=Rpath=Code/R_script/DMS_clustering_general.R,outpath=<Output path>,x=Examples/Methylation_matrix.txt ,dms=Examples/DMS.txt,dms_background=Examples/DMS_background.txt R_clustering_general.sbatch

```

### Recursive TRM identification on entire DMR matrix:

./recursive_motif_identification_noclustering.sh –i <path to input directory containing DMR files> -H <path to motif database> -r <fasta file> -target <DMR matrix> -background <DMR background matrix> -s <sequence extractor script> -R <path to R script files>

Example on demo data:

```
cd Code/Shell_script/shell_script_local/ 

./recursive_motif_identification_noclustering.sh -i Examples/ -H <path to motif database> –r mm10.fa -target Examples/DMR.txt -background Examples/DMR_background.txt –s Code/Perl_script/sequence_extractor.pl –R Code/R_script/
```

Explanation of command:
-i: path containing the DMS file
-H: path containing motif database 
-s:  sequence extractor Perl script
-R: path to R script files
-r: reference fasta file; examples are mm10.fa, hg19.fa etc.
-h: explains required parameters


### Recursive TRM identification on DMR cluster :

./recursive_motif_identification_clustering.sh –i <path to input directory containing DMR clusters> -H <path to motif database> -r <fasta file> -s <sequence extractor script> -R <path to R script files>

Example on demo data:

```
cd Code/Shell_script/shell_script_local/ 

./recursive_motif_identification_clustering.sh -i Examples/ -H <path to motif database> –r mm10.fa –s Code/Perl_script/sequence_extractor.pl –R Code/R_script/
```

Explanation of command:
-i: path containing the DMS clusters
-H: path containing motif database 
-s:  sequence extractor Perl script
-R: path to R script files
-r: reference fasta file; examples are mm10.fa, hg19.fa etc
-h: explains required parameters


