#!/bin/bash
###############################################################################
## Codes for the paper:
##   ..............
##
## Authors : GUERIN Pierre-Edouard, MATHON Laetitia
## Montpellier 2019-2020
## 
###############################################################################
## Usage:
##    bash 07_assignation/Assign_vsearch.sh
##
## Description:
##  ..............    
##
##
##
###############################################################################
## load config global variables
source 98_infos/config.sh


## Obitools
illuminapairedend=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" illuminapairedend"
obigrep=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obigrep"
ngsfilter=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" ngsfilter"
obisplit=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obisplit"
obiuniq=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obiuniq"
obiannotate=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obiannotate"
obiclean=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obiclean"
ecotag=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" ecotag"
obisort=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obisort"
obitab=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obitab"

container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"



## EDNAtools
usearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" usearch"

## Prefix for all generated files
pref="grinder_teleo1"
## Prefix of the final table, including the step and the program tested (ie: merging_obitools)
step=assign_sintax
## Path to the directory containing forward and reverse reads
R1_fastq="$DATA_PATH"/"$pref"/"$pref"_R1.fastq
R2_fastq="$DATA_PATH"/"$pref"/"$pref"_R2.fastq
## path to 'sample_description_file.txt'
sample_description_file=${INPUT_DATA}"/sample_description_file.txt"
# Chemin vers le fichier 'db_sim_teleo1.fasta'
refdb_dir=`pwd`"/07_assignation/db_teleo_sintax.fasta"
## path to outputs final and temporary (main)
main_dir=`pwd`"/07_assignation/Outputs/02_sintax/main"
fin_dir=`pwd`"/07_assignation/Outputs/02_sintax/final"


################################################################################################

## forward and reverse reads assembly
#assembly=${main_dir}"/"${pref}".fastq"
#$illuminapairedend -r ${R2_fastq} ${R1_fastq} > ${assembly}
## Remove non-aligned reads
#assembly_ali="${assembly/.fastq/.ali.fastq}"
#$obigrep -p 'mode!="joined"' ${main_dir}"/"${pref}".fastq" > ${assembly_ali}
## Assign each sequence to a sample
#identified="${assembly_ali/.ali.fastq/.ali.assigned.fasta}"
#unidentified="${assembly_ali/.ali.fastq/_unidentified.fastq}"
#$ngsfilter -t ${sample_description_file} -u ${unidentified} ${assembly_ali} --fasta-output > ${identified}
## Split big file into one file per sample 
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta ${identified}


all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
## Dereplication of reads in unique sequences
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo "/usr/bin/time $obiuniq -m sample "$sample" > "$dereplicated_sample > $sample_sh;
## Keep only sequences longer than 20pb with no ambiguous bases
good_sequence_sample="${dereplicated_sample/.fasta/.l20.fasta}"
echo "/usr/bin/time $obigrep -s '^[ACGT]+$' -l 20 "$dereplicated_sample" > "$good_sequence_sample >> $sample_sh
## Removal of PCR and sequencing errors (variants)
clean_sequence_sample="${good_sequence_sample/.fasta/.r005.clean.fasta}"
echo "/usr/bin/time $obiclean -r 0.05 -H "$good_sequence_sample" > "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
## Concatenation of all samples in one file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.l20.r005.clean.fasta > $all_sample_sequences_clean
## Dereplication in unique sequences
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
/usr/bin/time $obiuniq -m sample $all_sample_sequences_clean > $all_sample_sequences_uniq
## Removal of useless attributes in sequences headers
all_sample_sequences_ann="${all_sample_sequences_uniq/.fasta/.ann.fasta}"
$obiannotate --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
 --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount --delete-tag=obiclean_head  --delete-tag=obiclean_headcount \
 --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
 --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag --delete-tag=forward_tag --delete-tag=forward_score \
 --delete-tag=forward_primer --delete-tag=forward_match --delete-tag=tail_quality --delete-tag=mode --delete-tag=seq_a_single --delete-tag=status \
 --delete-tag=direction --delete-tag=seq_a_insertion --delete-tag=seq_b_insertion --delete-tag=seq_a_deletion --delete-tag=seq_b_deletion \
 --delete-tag=ali_length --delete-tag=head_quality --delete-tag=seq_b_single $all_sample_sequences_uniq > $all_sample_sequences_ann
## Sort sequences by 'count'
all_sample_sequences_sort="${all_sample_sequences_ann/.fasta/.sort.fasta}"
$obisort -k count --uppercase -r $all_sample_sequences_ann > $all_sample_sequences_sort
## converting lowercase letters to uppercase sequence ATGC letters in fasta file (usearch compatibility)
all_sample_sequences_sort_uppercase="${all_sample_sequences_sort/.fasta/.uppercase.fasta}"
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $0}}' $all_sample_sequences_sort > $all_sample_sequences_sort_uppercase
## Taxonomic assignation
all_sample_sequences_sintax_ann="${all_sample_sequences_sort/.fasta/.sintax_ann.csv}"
$usearch -sintax $all_sample_sequences_sort_uppercase -db $refdb_dir -sintax_cutoff 0.98 -strand plus -tabbedout $all_sample_sequences_sintax_ann
## convert usearch output to obifasta
all_sample_sequences_sintax_ann_fas="${all_sample_sequences_sintax_ann/.csv/.fasta}"
$container_python2 07_assignation/convert_assign_sintax_2_obifasta.py -f $all_sample_sequences_sort -s $all_sample_sequences_sintax_ann -o $all_sample_sequences_sintax_ann_fas
## Create final table
$obitab -o $all_sample_sequences_sintax_ann_fas > $fin_dir/"$step".csv

################################################################################################