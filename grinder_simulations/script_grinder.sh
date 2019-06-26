## define variables
#### number of cores
CORES=16
#### number of samples
NB_SAMPLE=30


for ID_SAMPLE in `seq -w 1 $NB_SAMPLE`; do
	## ((i=i%CORES)); ((i++==0)) && wait
	##  grinder command to run with input files into folder 'Inputs'
	bash grinder_simulations/simulation_grinder.sh "${ID_SAMPLE}" &
done
