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

## public
EDNATOOLS_SIMG="/99_utils/containers/ednatools.simg"
OBITOOLS_SIMG="/99_utils/containers/obitools.simg"
SINGULARITY_EXEC_CMD="singularity exec"



###############################################################################
## Input data

## path to folder of reads forward et reverse
DATA_PATH="/benchmark_simulated_dataset/00_Input_data/forward_reverse_reads/"


INPUT_DATA=`pwd`"/benchmark_simulated_dataset/00_Input_data/"

## path to reference database folder
REFDB_PATH=`pwd`"/benchmark_simulated_dataset/00_Input_data/reference_database"


###############################################################################

## number of available cores
CORES=16