#!/bin/bash

## commande grinder à mettre en boucle pour les 30 échantillons. Les fichiers entrée nécessaires sont dans le folder 'Inputs'
grinder -rf grinder_simulations/Inputs/seq_sample01.fasta \
    -tr 90000 -nl 12 -di 0 -af grinder_simulations/Inputs/abund_sample01.txt \
    -id 150 -rd 150 -fq 1 -ql 36 30 -mo FR \
    -dc '-' -un 1 -md poly4 3e-3 3.3e-8 -mr 98 2 \
    -hd Balzer -cp 1 -ck 0 -mi grinder_simulations/Inputs/tag_sample01.fasta \
    -bn sample01 -od grinder_simulations/Outputs/sample01/

## commande à mettre en boucle pour les échantillons 01 à 30 : copier les fichiers ranks.txt de chaque échantillon dans un folder ranks
cp grinder_simulations/Outputs/sample01/*ranks.txt > grinder_simulations/Outputs/ranks/

    ## commande à mettre en boucle pour les 12 réplicas de chacun des 30 échantillons (sample01-01, sample01-02.. sample30-12)
	grep '/1' grinder_simulations/Outputs/sample01/sample01-01-reads.fastq -A 3 > grinder_simulations/Outputs/grinder_teleo1_R1/sample01-01_R1.fastq
	grep '/2' grinder_simulations/Outputs/sample01/sample01-01-reads.fastq -A 3 > grinder_simulations/Outputs/grinder_teleo1_R2/sample01-01_R2.fastq


## commandes finales -> pas en boucle
cat grinder_simulations/Outputs/grinder_teleo1_R1/* > grinder_simulations/Outputs/grinder_teleo1/grinder_teleo1_R1.fastq
cat grinder_simulations/Outputs/grinder_teleo1_R2/* > grinder_simulations/Outputs/grinder_teleo1/grinder_teleo1_R2.fastq

## quand je concatène ça rajoute des lignes '--' qu'il faut enlever
sed -i -e '/^--/d' grinder_teleo1_R1.fastq
sed -i -e '/^--/d' grinder_teleo1_R2.fastq