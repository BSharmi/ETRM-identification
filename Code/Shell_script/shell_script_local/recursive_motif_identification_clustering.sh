#!/bin/bash

####################################################
# input parameters
Usage()  
{  
  echo -e "Usage: `basename $0` [-h/--help] <-i/--input-dir-to-DMS string> <-H/--HOMER-motif-path string> <-r/--reference-genome-fasta-path string> <-s/--sequence-extractor-path string> <-R/--Rscript-path string> <args>\n"; 
  exit 1;  
}  
 
idir=""
Hpath=""
refpath=""
Rpath=""
seqextractpath=""
  
# flexible order (Note: the space after ":")
#ARG=`getopt h:n:r:o: $*`  
  
# reset parameters  
#set --$ARG;  
  
# parse options of shell's input 
while getopts :h:i:H:r:s:R: PARAM_VAL; do  
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

if [ -z "$idir" ] || [ -z "$Hpath" ] || [ -z "$Rpath" ] || [ -z "$refpath" ] || [ -z "$seqextractpath" ]; then
    Usage
fi

############################################## recursive #########################
# cd to HOMER_custom
fpath=$Hpath
cd $fpath/data/knownTFs/motifs/

# set directory
directory=$idir


## count number of dms clusters
num_DMS=$(find $directory -mindepth 1 -maxdepth 1 -name 'cluster_*' -type d  | wc -l)

## get names of DMS folder using IFS. If this does not work  the two lines below
IFS=' ' read -r -a name_DMS <<< $(find $directory -mindepth 1  -maxdepth 1 -type d -name 'cluster_*')

## add HOMER path
PATH=$PATH:$Hpath/bin/

# find recursive motif for each DMS. background name used because there are target and input for each dms cluster. Using background avoids counting same dms twice
for ((idms=0;idms<=num_DMS-1;idms++)); do
	
	# read target bed file and convert to fasta
	perl $seqextractpath -r $refpath ${name_DMS[idms]}/*_dms.txt

	# read background bed file and convert to fasta
	perl $seqextractpath -r $refpath  ${name_DMS[idms]}/*_background.txt
	
	# run find motifs to find all known motifs
	findMotifs.pl ${name_DMS[idms]}/*_dms.withSquence.fa  mouse ${name_DMS[idms]}/ -mset vertebrates  -fastaBg  ${name_DMS[idms]}/*_background.withSquence.fa  -nomotif -p 10
	
	# delete seqbias if any sed does not work so use ed
	#sed -i '/^SeqB/d' ${name_DMS[idms]}/knownResults.txt
	ed -s ${name_DMS[idms]}/knownResults.txt <<< $'g/^SeqB/d\nw'
	ed -s ${name_DMS[idms]}/knownResults.txt <<< $'g/^Unkno/d\nw' 
	
	# check the top p-value of known motifs
	sigpval=$(awk 'FNR <= 2 {print $3}' ${name_DMS[idms]}/knownResults.txt | awk 'NR>1')
	
	
	## set parameters
	#topmotif=""
	topmotif_minus=""
	testpval=1e-9
	counter=""
	oldmotif_fa=${name_DMS[idms]}/*_dms.withSquence.fa
	
	## check recursively
	while [ $(perl -e "print $sigpval < $testpval") ]; do
		if [ "$counter" == "" ]; then
			## TF with the lowest p-value and highest number of sites		
			topmotif=$(awk 'NR>1 {printf ("%s\t%s\t%d\n", $1, $3, $6) | "sort -gk2,2g -nrk3,3nr"}' ${name_DMS[idms]}/knownResults.txt| awk 'NR==1 {print $1}' | awk -F'(' '{print $1}'| awk '{ print tolower($0) }')
			## create known motif library
			awk "BEGIN{IGNORECASE=1} FILENAME~ /^${topmotif:0:3}/" *.motif | awk '{print $FILENAME}' > ${name_DMS[idms]}/${topmotif:0:3}.search.motif
			## check if motif match was not found
			if [[ $(wc -l <${name_DMS[idms]}/${topmotif:0:3}.search.motif) -eq 0 ]]; then
				tmp=${topmotif:0:3}; tmp=${tmp//-/}
				awk "BEGIN{IGNORECASE=1} FILENAME~ /^${tmp}/" *.motif | awk '{print $FILENAME}' > ${name_DMS[idms]}/${topmotif:0:3}.search.motif
			fi
			findMotifs.pl  $oldmotif_fa  mouse ${name_DMS[idms]}/  -mset vertebrates  -fastaBg ${name_DMS[idms]}/*_background.withSquence.fa  -find ${name_DMS[idms]}/${topmotif:0:3}.search.motif > ${name_DMS[idms]}/${topmotif:0:3}.motif.txt
			# make directory
			mkdir ${name_DMS[idms]}/${topmotif:0:3}	
			mkdir ${name_DMS[idms]}/${topmotif:0:3}_minus
			##  read r script to parse locations, create two sequence files and 
			Rscript $Rpath/Recursive_motif_search_dplyr.R  ${name_DMS[idms]} ${topmotif:0:3} ${name_DMS[idms]}/${topmotif:0:3}.motif.txt ${name_DMS[idms]}/*_dms.txt   
		else
			topmotif=$(awk 'NR>1 {printf ("%s\t%s\t%d\n", $1, $3, $6) | "sort -gk2,2g -nrk3,3nr"}' ${name_DMS[idms]}/${topmotif:0:3}_minus/knownResults.txt| awk 'NR==1 {print $1}' | awk -F'(' '{print $1}'| awk '{ print tolower($0) }')
			## create known motif library
			awk "BEGIN{IGNORECASE=1} FILENAME~ /^${topmotif:0:3}/" *.motif | awk '{print $FILENAME}' > ${name_DMS[idms]}/${topmotif:0:3}.search.motif
			## check if motif match was not found
			if [[ $(wc -l <${name_DMS[idms]}/${topmotif:0:3}.search.motif) -eq 0 ]]; then
				tmp=${topmotif:0:3}; tmp=${tmp//-/}
				awk "BEGIN{IGNORECASE=1} FILENAME~ /^${tmp}/" *.motif | awk '{print $FILENAME}' > ${name_DMS[idms]}/${topmotif:0:3}.search.motif
			fi
			findMotifs.pl $oldmotif_fa  mouse ${name_DMS[idms]}/  -mset vertebrates  -fastaBg ${name_DMS[idms]}/*_background.withSquence.fa  -find ${name_DMS[idms]}/${topmotif:0:3}.search.motif > ${name_DMS[idms]}/${topmotif:0:3}.motif.txt
			# make directory
			mkdir ${name_DMS[idms]}/${topmotif:0:3}	
			mkdir ${name_DMS[idms]}/${topmotif:0:3}_minus
			##  read r script to parse locations, create two sequence files and 
			Rscript $Rpath/Recursive_motif_search_dplyr.R  ${name_DMS[idms]} ${topmotif:0:3} ${name_DMS[idms]}/${topmotif:0:3}.motif.txt $oldmotif_txt  		
		fi
			
		## create fasta sequence 
	    	perl $seqextractpath -r $refpath  ${name_DMS[idms]}/${topmotif:0:3}/${topmotif:0:3}.txt
        	## find secondary motifs by masking
	    	#findMotifs.pl ${name_DMS[idms]}/${topmotif:0:3}/${topmotif:0:3}.withSquence.fa  mouse ${name_DMS[idms]}/${topmotif:0:3}/  -mset vertebrates  -fastaBg ${name_DMS[idms]}/${name_DMS[idms]##*/}_background.withSquence.fa  -nomotif -p 10 -maskMotif ${name_DMS[idms]}/${topmotif:0:3}.search.motif
	    	findMotifs.pl ${name_DMS[idms]}/${topmotif:0:3}/${topmotif:0:3}.withSquence.fa  mouse ${name_DMS[idms]}/${topmotif:0:3}/  -mset vertebrates  -fastaBg ${name_DMS[idms]}/*_background.withSquence.fa  -nomotif -p 10 
        	## delete unwanted motifs sed does not work
	    	ed -s ${name_DMS[idms]}/${topmotif:0:3}/knownResults.txt <<< $'g/^SeqB/d\nw' 
        	ed -s ${name_DMS[idms]}/${topmotif:0:3}/knownResults.txt  <<< $'g/^Unkno/d\nw' 
            
       		## check if any sites remain
        	if [ "$(ls -A ${name_DMS[idms]}/${topmotif:0:3}_minus)" ]; then 
        		## create fasta for setdiff
          		perl $seqextractpath -r $refpath  ${name_DMS[idms]}/${topmotif:0:3}_minus/${topmotif:0:3}_minus.txt
            		## find new candidate TF motifs
            		findMotifs.pl ${name_DMS[idms]}/${topmotif:0:3}_minus/${topmotif:0:3}_minus.withSquence.fa  mouse ${name_DMS[idms]}/${topmotif:0:3}_minus/  -mset vertebrates  -fastaBg ${name_DMS[idms]}/*_background.withSquence.fa  -nomotif -p 10
            		## delete unwanted motifs sed does not work
			#sed -i '/^SeqB/d' ${name_DMS[idms]}/${topmotif:0:3}_minus/knownResults.txt
			ed -s ${name_DMS[idms]}/${topmotif:0:3}_minus/knownResults.txt <<< $'g/^SeqB/d\nw'
			ed -s ${name_DMS[idms]}/${topmotif:0:3}_minus/knownResults.txt <<< $'g/^Unkno/d\nw'  			
            		## check the top p-value of remaining known motifs
            		sigpval=$(awk 'FNR <= 2 {print $3}' ${name_DMS[idms]}/${topmotif:0:3}_minus/knownResults.txt| awk 'NR>1')
            		## increase counter
            		(( counter++ ))
            		## update oldmotif
            		oldmotif_fa=${name_DMS[idms]}/${topmotif:0:3}_minus/${topmotif:0:3}_minus.withSquence.fa
            		oldmotif_txt=${name_DMS[idms]}/${topmotif:0:3}_minus/${topmotif:0:3}_minus.txt
        	else 
            		sigpval=1 
        	fi    	
	done	
	## plot graph
	## plot top TFS from all known results
	Rscript $Rpath/Recursive_motif_plot.R ${name_DMS[idms]}/knownResults.txt "all" ${name_DMS[idms]}/ "all"
	## get file names
	num_TRMs=$(find ${name_DMS[idms]}/ -mindepth 1 -maxdepth 1 -name '*_minus*' -type d  | wc -l)
	if [ $((num_TRMs)) -ne 0 ] ; then 
		for dir in ${name_DMS[idms]}/*minus*; do 
			## get folder name
			dname=${dir##/*/}
			tname=${dname%_*}
			knownfile=${name_DMS[idms]}/$tname/knownResults.txt
			if [ -f "$knownfile" ]; then 
				Rscript $Rpath/Recursive_motif_plot.R ${name_DMS[idms]}/$tname/knownResults.txt "module" ${name_DMS[idms]}/ $tname
				#head $directory/HOMER_1/$tname/knownResults.txt
			fi
		done
	fi	
done


