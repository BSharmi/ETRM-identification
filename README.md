# Recursive Transcriptional Regulatory Module (TRM) Identification
The tool is designed for the identification of transcriptional regulatory modules (TRMs) enriched within differentially methylated regions (DMRs). Given a set of DMRs and background regions, it uses HOMER’s known motif enrichment functionality to recursively identify TF motifs significantly enriched. The user has the option to cluster the DMRs into different subsets first or directly apply recursive TRM identification on the entire set of DMRs. Clustering is based on WGCNA (Weighted Gene Correlation Network Analysis) algorithm but other clustering algorithms such as K-means may be used as well. 

## Methylated motifs
The user may use an existing motif database \
OR motifs compiled in this study: copy and unzip the file "motifs.zip" to HOMER library/data/knownTFs/motifs folder \
OR append the HOMER motif file together with the motif files provided: \
For example, to add methylated motifs to HOMER's known motif library of vertebrates, please copy the motifs (motifs.zip) to HOMER library/data/knownTFs/motifs and unzip the file. \
The user can also append the motifs to an existing motif file. For e.g. user can append the methylated motifs to HOMER's known motif files for vertebrates (/home/bsharmi6/HOMER_custom/data/knownTFs/vertebrates/known.motifs). The appended motif files are also provided under Motif_db/all.motifs and Motif_db/known.motifs for convenience. The user can replace the known vertebrates motif files in HOMER with these two motif files. 

An snapshot of the motif database containing methylated motifs in addition to HOMER's known motifs (e.g. path = /home/bsharmi6/HOMER_custom/data/knownTFs/motifs/) is given below - 

![Motif files](https://github.com/BSharmi/TRM-identification/blob/master/Images/motif.png)

## Required R libraries: 
"dplyr", "WGCNA", "pheatmap", "igraph", "gdata"

## The package contains the following modules-
1.	Code/Shell_script/shell_script_local: to run the commands on a local system
2.	Code/Shell_script/shell_script_server: for users who wish to run the tool on a server
3.	Code/R_script: R scripts needed by the tool
4.	Examples: sample data and output provided as demo
5.	Motif database: Motif libraries compiled from HOMER and other public databases (for details see the paper)

Please download a genome fa file (e.g. mm10.fa). It is a large file and has not been provided here.
Clone the Github repository to a local system or a server

## Below are some command line arguments to perform different functions:

### Cluster DMR 

./run_DMR_clustering_script.sh -R <path to R script files> –o <path to output directory> -d <3 column BED format DMR file> -x <Methylome matrix> -db <3 column BED format DMR background file> 

Example on demo data:

```
cd Code/Shell_script/shell_script_local/ 

./run_DMR_clustering_script.sh -R Code/R_script/DMR_clustering_general.R -o Examples/Output/ -d Examples/DMR.txt -x Examples/Methylation_matrix.txt -b Examples/DMR_background.txt 
```

Explanation of command:
-o: path to store the output \
-d: tab delimited three column BED format DMR file with three columns: chromosome, start, end \
-x: methylome matrix. One column indicates methylome values for DMR location, one column indicates methylome values for one sample. Number of rows of x must match with the number of rows in the DMR file \
-b: tab delimited three column BED format DMR background file to create background regions for each DMR cluster \
-R: R script used for clustering \
-h: explains required parameters \
The command will cluster DMRs, and create a folder for each cluster name containing the tab delimited DMRs files and the corresponding background files

#### Running scripts on server

Since WGCNA can be slow on large DMR matrices, it might be efficient to run the scripts by submitting as jobs on a high performance computing system. The script below shows how to run the scripts on server -  

```
sbatch --export=Rpath=Code/R_script/DMR_clustering_general.R,outpath=Examples/Output/,x=Examples/Methylation_matrix.txt,dms=Examples/DMR.txt,dms_background=Examples/DMR_background.txt R_clustering_general.sbatch

```

### Recursive TRM identification on entire DMR matrix:

./recursive_motif_identification_noclustering.sh –i <path to input directory containing DMR files> -H <path to motif database> -r <fasta file> -target <DMR matrix> -background <DMR background matrix> -s <sequence extractor script> -R <path to R script files>

Example on demo data:

```
cd Code/Shell_script/shell_script_local/ 

./recursive_motif_identification_noclustering.sh -i Examples/ -H /home/bsharmi6/HOMER_custom/ -r mm10.fa -t DMR.txt -b DMR_background.txt -s Code/Perl_script/sequence_extractor.pl -R Code/R_script/
```

Explanation of command:
-i: path containing the DMR file \
-H: path to HOMER database \
-s: sequence extractor Perl script \
-R: path to R script files \
-t: DMR matrix \
-b: DMR background matrix \
-r: path containing the reference fasta file; examples are mm10.fa, hg19.fa etc. \
-h: explains required parameters \

#### Running scripts on server

```
sbatch --export=idir=Examples/,Hpath=/home/bsharmi6/HOMER_custom/,refpath=/home/bsharmi6/mm10bowtie2/mm10.fa,seqextractpath=Code/Perl_script/sequence_extractor.pl,Rpath=Code/R_script/,target=DMR.txt,background=DMR_background.txt recursive_motif_identification_noclustering.sbatch 
```


### Recursive TRM identification on DMR cluster :

./recursive_motif_identification_clustering.sh –i <path to input directory containing DMR clusters> -H <path to motif database> -r <fasta file> -s <sequence extractor script> -R <path to R script files>

Example on demo data:

```
cd Code/Shell_script/shell_script_local/ 

./recursive_motif_identification_clustering.sh -i Examples/ -H /home/bsharmi6/HOMER_custom/ -r mm10.fa -s Code/Perl_script/sequence_extractor.pl -R Code/R_script/
```

Explanation of command:
-i: path containing the DMR clusters \
-H: path to HOMER database  \
-s: sequence extractor Perl script \
-R: path to R script files \
-r: path containing the reference fasta file; examples are mm10.fa, hg19.fa etc \
-h: explains required parameters

#### Running scripts on server

```
sbatch --export=idir=Examples/,Hpath=/home/bsharmi6/HOMER_custom/,refpath=/home/bsharmi6/mm10bowtie2/mm10.fa,seqextractpath=Code/Perl_script/sequence_extractor.pl,Rpath=Code/R_script/ recursive_motif_identification_clustering.sbatch
```

Contact: bsharmi6@vt.edu
