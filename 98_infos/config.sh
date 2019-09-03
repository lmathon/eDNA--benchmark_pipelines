###############################################################################
##
## The config.sh file contains modifiable parameters (and non modifiable too).
## In normal circumstances it should be the only file which may be edited 
## by the end-user. 
##
###############################################################################
##
###############################################################################
## singularity containers path

## mbb
EDNATOOLS_SIMG="/share/reservebenefit/utils/conteneurs/ednatools.simg"
OBITOOLS_SIMG="/share/reservebenefit/utils/conteneurs/obitools.simg"

## Singularity mount point
SING_MNT="/depot"

## Singularity exec command
#SINGULARITY_EXEC_CMD="singularity exec -B .:"${SING_MNT}

COMD="illuminapairedend"
SINGULARITY_EXEC_CMD="singularity exec -B .:${SING_MNT} ${OBITOOLS_SIMG}" 



#OBITOOLS_SIMG="obitools.simg"



###############################################################################
## Input data

## path to folder of reads forward et reverse
DATA_PATH='00_Input_data/forward_reverse_reads'
## path to reference database folder
REFDB_PATH="00_Input_data/reference_database/"


###############################################################################