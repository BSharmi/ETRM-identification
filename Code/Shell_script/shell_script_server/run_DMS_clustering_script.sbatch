#!/bin/bash

#SBATCH --time=100:00:00


#SBATCH --nodes=1 --ntasks-per-node=1

#SBATCH --mem=122280

#SBATCH -p pegasus_q

#SBATCH -A other_research

# (Uncomment the following line if you wish to use it.)
# module add R/3.4.3-foss-2017b-X11-20171023
module add R/3.4.4-foss-2018a-X11-20180131

# Below here enter the commands to start your job

## check if any input parameter is empty 
if [[ -z "$outpath" ]] || [[ -z "$x" ]] || [[ -z "$dms" ]] || [[ -z "$dms_background" ]] || [[ -z "$Rpath" ]] ; then 
	echo -e "Usage: `basename $0` [-h/--help] <-outpath/--output-file-path string> <-x/--dms-three-column-BED-format string> <-dms_background/--dms-background-three-column-BED-format string> <-Rpath/--Rscript-path string> <args>\n";  
	exit 1; 
fi


######### run R script
Rscript $Rpath $outpath $x $dms $dms_background

