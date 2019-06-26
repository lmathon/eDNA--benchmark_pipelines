## define variables
#### number of cores
CORES=16
#### number of samples
NB_SAMPLE=30

seq -w 1 $NB_SAMPLE | xargs -P $CORES simulation_grinder.sh
