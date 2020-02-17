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


## Singularity exec command
#SINGULARITY_EXEC_CMD="singularity exec -B .:"${SING_MNT}

SINGULARITY_EXEC_CMD="singularity exec --bind /share:/share" 


### home laetitia
#SINGULARITY_EXEC_CMD="singularity exec"
#OBITOOLS_SIMG="/home/lmathon/obitools.img"



#OBITOOLS_SIMG="obitools.simg"



###############################################################################
## Input data

## path to folder of reads forward et reverse
DATA_PATH="/share/reservebenefit/working/Input_data/Outputs"
## path to data to test
#DATA_PATH="/share/reservebenefit/working/Input_data/Outputs"
## path to tiny data only for debug
#DATA_PATH="/share/reservebenefit/working/Input_data/tiny"

INPUT_DATA=`pwd`"/00_Input_data/"

## path to reference database folder
REFDB_PATH=`pwd`"/00_Input_data/reference_database/"


## home laetitia
#REFDB_PATH="/home/lmathon/reference_database"

###############################################################################

## number of available cores
CORES=1