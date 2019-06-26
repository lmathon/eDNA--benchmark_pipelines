## define variables
#### number of cores
CORES=16
#### number of samples
NB_SAMPLE=30
#### 
#MONTAGE=""
MONTAGE="/simulations/"

for ID_SAMPLE in `seq -w 1 $NB_SAMPLE`; 
do echo "${ID_SAMPLE}" 
##  grinder command to run with input files into folder 'Inputs'
grinder -rf "${MONTAGE}"grinder_simulations/Inputs/seq_sample"${ID_SAMPLE}".fasta \
    -tr 90000 -nl 12 -di 0 -af "${MONTAGE}"grinder_simulations/Inputs/abund_sample"${ID_SAMPLE}".txt \
    -id 150 -rd 150 -fq 1 -ql 36 30 -mo FR \
    -dc '-' -un 1 -md poly4 3e-3 3.3e-8 -mr 98 2 \
    -hd Balzer -cp 1 -ck 0 -mi "${MONTAGE}"grinder_simulations/Inputs/tag_sample"${ID_SAMPLE}".fasta \
    -bn sample"${ID_SAMPLE}" -od "${MONTAGE}"grinder_simulations/Outputs/sample"${ID_SAMPLE}"
done
