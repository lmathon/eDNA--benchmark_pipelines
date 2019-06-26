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

"""
## commande à mettre en boucle pour les échantillons 01 à 30 : copier les fichiers ranks.txt de chaque échantillon dans un folder ranks
cp grinder_simulations/Outputs/sample"${ID_SAMPLE}"/*ranks.txt grinder_simulations/Outputs/ranks/


bash deinterleave_fastq.sh < interleaved.fastq f.fastq r.fastq compress


    ## commande à mettre en boucle pour les 12 réplicas de chacun des 30 échantillons (sample01-01, sample01-02.. sample30-12)
	grep '/1' grinder_simulations/Outputs/sample01/sample01-01-reads.fastq -A 3 > grinder_simulations/Outputs/grinder_teleo1_R1/sample01-01_R1.fastq
	grep '/2' grinder_simulations/Outputs/sample01/sample01-01-reads.fastq -A 3 > grinder_simulations/Outputs/grinder_teleo1_R2/sample01-01_R2.fastq


## commandes finales -> pas en boucle
cat grinder_simulations/Outputs/grinder_teleo1_R1/* > grinder_simulations/Outputs/grinder_teleo1/grinder_teleo1_R1.fastq
cat grinder_simulations/Outputs/grinder_teleo1_R2/* > grinder_simulations/Outputs/grinder_teleo1/grinder_teleo1_R2.fastq

## quand je concatène ça rajoute des lignes '--' qu'il faut enlever
sed -i -e '/^--/d' grinder_teleo1_R1.fastq
sed -i -e '/^--/d' grinder_teleo1_R2.fastq
"""