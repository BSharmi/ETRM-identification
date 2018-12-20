#!/bin/bash


####################################################
# input parameters
Usage()  
{  
  echo -e "Usage: `basename $0` [-h/--help] <-i/--input-dir-to-DMS string> <-H/--HOMER-motif-path string> <-r/--reference-genome-fasta-path string>  <-t/--DMR-matrix-path string> <-b/--DMR-background-path string> <-s/--sequence-extractor-path string> <-R/--Rscript-path string> <args>\n"; 
  exit 1;  
}  
 
idir=""
Hpath=""
refpath=""
Rpath=""
target=""
background=""
seqextractpath=""
  
# flexible order (Note: the space after ":")
#ARG=`getopt h:n:r:o: $*`  
  
# reset parameters  
#set --$ARG;  
  
# parse options of shell's input 
while getopts :h:i:H:r:t:b:s:R: PARAM_VAL; do  
  case "${PARAM_VAL}" in  
  i|input-dir-to-DMS)  
    idir=$OPTARG;
    ;;  
  R|R-script-path)    
    Rpath=$OPTARG;  
    ;; 
  H|HOMER-motif-path)    
    Hpath=$OPTARG;  
    ;; 
  r|reference-genome-fasta-path)    
    refpath=$OPTARG;  
    ;;
  t|DMR-matrix-path)    
    target=$OPTARG;  
    ;; 
  b|DMR-background-path)    
    background=$OPTARG;  
    ;;
  s|sequence-extractor-path)  
    seqextractpath=$OPTARG;  
    ;;  
  h|help)  
    Usage
    ;;  
  *)  
    Usage
    ;;  
  esac  
done  

shift $(($OPTIND -1))

if [ -z "$idir" ] || [ -z "$Hpath" ] || [ -z "$Rpath" ] || [ -z "$refpath" ] || [ -z "$target" ] || [ -z "$background" ] || [ -z "$seqextractpath" ]; then
    Usage
fi

####################################################
## configure parameters
#target="TET1KO10x_all_dms.txt" #$1 target DMS 
#background="TET1KO10x_all_dms_background.txt" # background DMS
####################################################

################################################### recirsive #################

## add HOMER path
PATH=$PATH:$Hpath/bin/

# cd to HOMER_custom
fpath=$Hpath
cd $fpath/data/knownTFs/motifs/

## set directory
directory=$idir

# read target bed file and convert to fasta
perl $seqextractpath -r $refpath $directory/$target

# read background bed file and convert to fasta
perl $seqextractpath -r $refpath $directory/$background

# run find motifs to find all known motifs
findMotifs.pl $directory/${target%*.txt}.withSquence.fa  mouse $directory/ -mset vertebrates  -fastaBg  $directory/${background%*.txt}.withSquence.fa  -p 10 -nomotif

# delete seqbias if any sed does not work so use ed
ed -s $directory/knownResults.txt <<< $'g/^SeqB/d\nw' 
ed -s $directory/knownResults.txt <<< $'g/^Unkno/d\nw' 

# check the top p-value of known motifs
sigpval=$(awk 'FNR <= 2 {print $3}' $directory/knownResults.txt | awk 'NR>1')

## set parameters
opmotif_minus=""
testpval=1e-10
counter=""
oldmotif_fa=$directory/${target%*.txt}.withSquence.fa

## check recursively
while [ $(perl -e "print $sigpval < $testpval") ]; do

	if [ "$counter" == "" ]; then
		## delete unwanted motifs sed doesnt work
		ed -s $directory/knownResults.txt <<< $'g/^SeqB/d\nw' 
		ed -s $directory/knownResults.txt <<< $'g/^Unkno/d\nw' 
		## TF with the lowest p-value and highest sumber of sites
		topmotif=$(awk 'NR>1 {printf ("%s\t%s\t%d\n", $1, $3, $6) | "sort -gk2,2g -nrk3,3nr"}' $directory/knownResults.txt| awk 'NR==1 {print $1}' | awk -F'(' '{print $1}'| awk '{ print tolower($0) }')
		## update topmotif replace -
		#topmotif=${topmotif//-/}
		## create known motif library
		awk "BEGIN{IGNORECASE=1} FILENAME~ /^${topmotif:0:3}/" *.motif | awk '{print $FILENAME}' > $directory/${topmotif:0:3}.search.motif
		## check if motif match was not found
		if [[ $(wc -l <$directory/${topmotif:0:3}.search.motif) -eq 0 ]]; then
			tmp=${topmotif:0:3}; tmp=${tmp//-/}
			awk "BEGIN{IGNORECASE=1} FILENAME~ /^${tmp}/" *.motif | awk '{print $FILENAME}' > $directory/${topmotif:0:3}.search.motif
		fi
	    findMotifs.pl  $oldmotif_fa  mouse $directory  -mset vertebrates  -fastaBg $directory/${background%*.txt}.withSquence.fa  -find $directory/${topmotif:0:3}.search.motif > $directory/${topmotif:0:3}.motif.txt
		# make directory
		mkdir $directory/${topmotif:0:3}	
		mkdir $directory/${topmotif:0:3}_minus
		##  read r script to parse locations, create two sequence files and 
		Rscript $Rpath/Recursive_motif_search_dplyr.R  $directory ${topmotif:0:3} $directory/${topmotif:0:3}.motif.txt $directory/$target  ### search for similar motif name
	else
		## delete unwanted motifs sed does not work
		ed -s $directory/${topmotif:0:3}_minus/knownResults.txt <<< $'g/^SeqB/d\nw' 
		ed -s $directory/${topmotif:0:3}_minus/knownResults.txt <<< $'g/^Unkno/d\nw'
		topmotif=$(awk 'NR>1 {printf ("%s\t%s\t%d\n", $1, $3, $6) | "sort -gk2,2g -nrk3,3nr"}' $directory/${topmotif:0:3}_minus/knownResults.txt| awk 'NR==1 {print $1}' | awk -F'(' '{print $1}'| awk '{ print tolower($0) }')
		## update topmotif replace -
		#topmotif=${topmotif//-/}
		## create known motif library
		awk "BEGIN{IGNORECASE=1} FILENAME~ /^${topmotif:0:3}/" *.motif | awk '{print $FILENAME}' > $directory/${topmotif:0:3}.search.motif
		## check if motif match was not found
		if [[ $(wc -l <$directory/${topmotif:0:3}.search.motif) -eq 0 ]]; then
			tmp=${topmotif:0:3}; tmp=${tmp//-/}
			awk "BEGIN{IGNORECASE=1} FILENAME~ /^${tmp}/" *.motif | awk '{print $FILENAME}' > $directory/${topmotif:0:3}.search.motif
		fi
		findMotifs.pl $oldmotif_fa  mouse $directory/  -mset vertebrates  -fastaBg $directory/${background%*.txt}.withSquence.fa  -find $directory/${topmotif:0:3}.search.motif > $directory/${topmotif:0:3}.motif.txt
		# make directory
		mkdir $directory/${topmotif:0:3}	
		mkdir $directory/${topmotif:0:3}_minus
		##  read r script to parse locations, create two sequence files and 
		Rscript $Rpath/Recursive_motif_search_dplyr.R  $directory ${topmotif:0:3} $directory/${topmotif:0:3}.motif.txt $oldmotif_txt  ### this is for sing;e motif input to search for similar motif name	
	fi

	## create fasta sequence 
	perl $seqextractpath -r $refpath  $directory/${topmotif:0:3}/${topmotif:0:3}.txt
	## find secondary motifs
	findMotifs.pl $directory/${topmotif:0:3}/${topmotif:0:3}.withSquence.fa  mouse $directory/${topmotif:0:3}/  -mset vertebrates  -fastaBg $directory/${background%*.txt}.withSquence.fa  -nomotif -p 10

	## check if any sites remain
	if [ "$(ls -A $directory/${topmotif:0:3}_minus)" ]; then 
		## create fasta for setdiff
		perl $seqextractpath -r $refpath  $directory/${topmotif:0:3}_minus/${topmotif:0:3}_minus.txt
		# rerun HOMER
		findMotifs.pl $directory/${topmotif:0:3}_minus/${topmotif:0:3}_minus.withSquence.fa  mouse $directory/${topmotif:0:3}_minus/  -mset vertebrates  -fastaBg $directory/${background%*.txt}.withSquence.fa  -nomotif -p 10
		## remove seqbias from output sed does now work
		ed -s $directory/${topmotif:0:3}/knownResults.txt <<< $'g/^SeqB/d\nw' 
		ed -s $directory/${topmotif:0:3}_minus/knownResults.txt  <<< $'g/^SeqB/d\nw' 
		# check the top p-value of remaining known motifs
 		sigpval=$(awk 'FNR <= 2 {print $3}' $directory/${topmotif:0:3}_minus/knownResults.txt| awk 'NR>1')
		## increase counter
		(( counter++ ))
		## update oldmotif
		oldmotif_fa=$directory/${topmotif:0:3}_minus/${topmotif:0:3}_minus.withSquence.fa
		oldmotif_txt=$directory/${topmotif:0:3}_minus/${topmotif:0:3}_minus.txt
	else
		sigpval=1
	fi 	
done

## plot graph
## plot top TFS from all known results
Rscript $Rpath/Recursive_motif_plot.R $directory/knownResults.txt "all" $directory "all"
## get file names
num_TRMs=$(find $directory/ -mindepth 1 -maxdepth 1 -name '*_minus*' -type d  | wc -l)
if [ $((num_TRMs)) -ne 0 ] ; then 
	for dir in $directory/*minus*; do 
		## get folder name
		dname=${dir##/*/}
		tname=${dname%_*}
		knownfile=$directory/$tname/knownResults.txt
		if [ -f "$knownfile" ]; then 
			Rscript $Rpath/Recursive_motif_plot.R $directory/$tname/knownResults.txt "module" $directory $tname			
		fi
	done
fi
