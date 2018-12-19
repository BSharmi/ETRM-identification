#!/bin/bash

####################################################
# input parameters
Usage()  
{  
  echo -e "Usage: `basename $0` [-h/--help] <-o/--output-dir string> <-d/--dms-3-columns-BED-file string> <-x/--methylome-matrix string> <-db/--dms-background-3-columnd-BED-file string> <-R/--R-script-path string> <args>\n";  
  exit 1;  
}  
 
dms=""
x=""
outpath=""
Rpath=""
dms_background=""
  
# flexible order (Note: the space after ":")
#ARG=`getopt h:n:r:o: $*`  
  
# reset parameters  
#set --$ARG;  
  
# parse options of shell's input 
while getopts :h:R:d:o:x:db: PARAM_VAL; do  
  case "${PARAM_VAL}" in  
  o|output-dir)  
    outpath=$OPTARG;
    ;;  
  R|R-script-path)    
    Rpath=$OPTARG;  
    ;; 
  x|methylome-matrix)    
    x=$OPTARG;  
    ;; 
  db|dms-background-3-columnd-BED-file)    
    dms_background=$OPTARG;  
    ;;
  d|dms-3-columns-BED-file)  
    dms=$OPTARG;  
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

if [ -z "$x" ] || [ -z "$outpath" ] || [ -z "$Rpath" ] || [ -z "$dms" ] || [ -z "$dms_background" ]; then
    Usage
fi

####################################################
######### run R script
echo $outpath
Rscript $Rpath $outpath $x $dms $dms_background
