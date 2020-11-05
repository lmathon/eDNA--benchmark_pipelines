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
fastp=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" fastp"

# Prefix for all generated files
pref=Carry
# Prefix of the final table 
step=merging_fastp
# Path to forward and reverse fastq files
R1_fastq="${DATA_PATH}"/"$pref"/"$pref"_R1.fastq.gz
R2_fastq="${DATA_PATH}"/"$pref"/"$pref"_R2.fastq.gz
# Path to file 'sample_description_file.txt'
sample_description_file=${INPUT_DATA}"/sample_description_file.txt"
# Path to file 'db_sim_teleo1.fasta'
refdb_dir=${REFDB_PATH}"/db_carry_all.fasta"
# Path to all files 'embl' of the reference database
base_dir=${REFDB_PATH}
### Prefix of the ref database files must not contain "." ou "_"
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//g' | awk -F/ '{print $NF}' | uniq`
# Path to intermediate and final folders
main_dir=`pwd`"/benchmark_real_dataset/01_merging/Outputs/06_fastp/main"
fin_dir=`pwd`"/benchmark_real_dataset/01_merging/Outputs/06_fastp/final"


###################################################################################################################

## forward and reverse reads assembly
assembly=${main_dir}"/"${pref}".fastq"
/usr/bin/time $fastp -i ${R1_fastq} -I ${R2_fastq} --merge --overlap_len_require 10 --overlap_diff_limit 15 --overlap_diff_percent_limit 25 --disable_length_filtering -w 16 --merged_out ${assembly}

## Assign each sequence to a sample
identified="${assembly/.fastq/.ali.assigned.fasta}"
unidentified="${assembly/.fastq/_unidentified.fastq}"
$ngsfilter -t ${sample_description_file} -u ${unidentified} ${assembly} --fasta-output > ${identified}
## Split big file into one file per sample
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta ${identified}

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
# Dereplicate reads into unique sequences
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo "/usr/bin/time $obiuniq -m sample "$sample" > "$dereplicated_sample > $sample_sh;
# Keep sequences longuer than 20bp without ambiguous bases
good_sequence_sample="${dereplicated_sample/.fasta/.l20.fasta}"
echo "/usr/bin/time $obigrep -s '^[ACGT]+$' -l 20 "$dereplicated_sample" > "$good_sequence_sample >> $sample_sh
# Removal of PCR and sequencing errors
clean_sequence_sample="${good_sequence_sample/.fasta/.r005.clean.fasta}"
echo "/usr/bin/time $obiclean -r 0.05 -H "$good_sequence_sample" > "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
# Concatenation of all samples into one file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.l20.r005.clean.fasta > $all_sample_sequences_clean
# Dereplication into unique sequences
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
/usr/bin/time $obiuniq -m sample $all_sample_sequences_clean > $all_sample_sequences_uniq
# Taxonomic assignation
all_sample_sequences_tag="${all_sample_sequences_uniq/.fasta/.tag.fasta}"
/usr/bin/time $ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $all_sample_sequences_uniq -m 0.98 > $all_sample_sequences_tag
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