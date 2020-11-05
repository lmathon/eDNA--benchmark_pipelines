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
##    bash obitools_reference/total_obitools.sh
##
## Description:
##  ..............    
##
##
##
###############################################################################
## load config global variables
source benchmark_real_dataset/98_infos/config.sh


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
## EDNAtools
fastp=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" fastp"
container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"


## Prefix for all generated files
pref="Carry"
## Prefix of the final table, including the step and the program tested (ie: merging_obitools)
step="filter_fastp"
## Path to the directory containing forward and reverse reads
R1_fastq="${DATA_PATH}"/"$pref"/"$pref"_R1.fastq.gz
R2_fastq="${DATA_PATH}"/"$pref"/"$pref"_R2.fastq.gz
## path to 'sample_description_file.txt'
sample_description_file=${INPUT_DATA}"/sample_description_file.txt"
## path to the file 'db_sim_teleo1.fasta'
refdb_dir=${REFDB_PATH}"/db_carry_all.fasta"
## Path to embl files of the reference database
base_dir=${REFDB_PATH}
### remove '.' and  '_' from the prefix files
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//'g | awk -F/ '{print $NF}' | uniq`
## path to outputs final and temporary (main)
main_dir=`pwd`"/benchmark_real_dataset/04_filtering/Outputs/05_fastp/main"
fin_dir=`pwd`"/benchmark_real_dataset/04_filtering/Outputs/05_fastp/final"


################################################################################################
## forward and reverse reads assembly
assembly=${main_dir}"/"${pref}".fastq"
$illuminapairedend -r ${R2_fastq} ${R1_fastq} > ${assembly}
## Remove non-aligned reads
assembly_ali="${assembly/.fastq/.ali.fastq}"
$obigrep -p 'mode!="joined"' ${main_dir}"/"${pref}".fastq" > ${assembly_ali}
# Keep only sequences longer than 20pb with no N bases
filtered="${assembly_ali/.ali.fastq/.ali.l20.fastq}"
/usr/bin/time $fastp -i ${assembly_ali} -n 0 -l 20 -q 0 -u 100 -w 16 -o ${filtered}
## Assign each sequence to a sample
identified="${filtered/.l20.fastq/.l20.assigned.fasta}"
unidentified="${filtered/.l20.fastq/_unidentified.fastq}"
$ngsfilter -t ${sample_description_file} -u ${unidentified} ${filtered} --fasta-output > ${identified}
## Split big file into one file per sample 
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta ${identified}

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
# Dereplicate reads in unique sequences
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo "$obiuniq -m sample "$sample" > "$dereplicated_sample > $sample_sh;
# Removal of PCR and sequencing errors (variants)
clean_sequence_sample="${dereplicated_sample/.fasta/.r005.clean.fasta}"
echo "$obiclean -r 0.05 -H "$dereplicated_sample" > "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
# Concatenate all files into one main file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.r005.clean.fasta > $all_sample_sequences_clean
# Dereplicate in unique sequences
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
$obiuniq -m sample $all_sample_sequences_clean > $all_sample_sequences_uniq
# Taxonomic assignation
all_sample_sequences_tag="${all_sample_sequences_uniq/.fasta/.tag.fasta}"
$ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $all_sample_sequences_uniq -m 0.98 > $all_sample_sequences_tag
# Removal of unnecessary attributes in sequence headers
all_sample_sequences_ann="${all_sample_sequences_tag/.fasta/.ann.fasta}"
$obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
 --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
 --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
 --delete-tag=obiclean_head  --delete-tag=obiclean_headcount \
 --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status \
 --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
 --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag \
 --delete-tag=forward_tag --delete-tag=forward_score --delete-tag=forward_primer --delete-tag=forward_match \
 --delete-tag=tail_quality --delete-tag=mode --delete-tag=seq_a_single $all_sample_sequences_tag > $all_sample_sequences_ann
# Sort sequences by 'count'
all_sample_sequences_sort="${all_sample_sequences_ann/.fasta/.sort.fasta}"
$obisort -k count -r $all_sample_sequences_ann > $all_sample_sequences_sort
# Create the final table
$obitab -o $all_sample_sequences_sort > $fin_dir/"$step".csv

gzip $main_dir/*