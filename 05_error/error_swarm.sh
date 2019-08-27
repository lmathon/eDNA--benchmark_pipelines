#!/bin/bash
##Obitools

illuminapairedend='singularity exec 99_utils/images/obitools.img illuminapairedend'
obigrep='singularity exec 99_utils/images/obitools.img obigrep'
ngsfilter='singularity exec 99_utils/images/obitools.img ngsfilter'
obisplit='singularity exec 99_utils/images/obitools.img obisplit'
obiuniq='singularity exec 99_utils/images/obitools.img obiuniq'
obiannotate='singularity exec 99_utils/images/obitools.img obiannotate'
obiclean='singularity exec 99_utils/images/obitools.img obiclean'
ecotag='singularity exec 99_utils/images/obitools.img ecotag'
obisort='singularity exec 99_utils/images/obitools.img obisort'
obitab='singularity exec 99_utils/images/obitools.img obitab'
swarm='singularity exec 99_utils/images/ednatools.img swarm'


# Path to directory containing forward et reverse reads
DATA_PATH='00_Input_data/forward_reverse_reads'
# Prefix for all generated files
pref=grinder_teleo1
# Prefix for the final table, containing the step and the program used (ex: merging_obitools)
step=error_swarm
# Files containing forward et reverse reads
R1_fastq="$DATA_PATH"/"$pref"_R1.fastq
R2_fastq="$DATA_PATH"/"$pref"_R2.fastq
# Path to the file 'sample_description_file.txt'
sample_description_file='00_Input_data/sample_description_file.txt'
# Path to the file 'db_sim_teleo1.fasta'
refdb_file='00_Input_data/reference_database/db_sim_teleo1.fasta'
# Path to the files 'embl' of the reference database
base_dir='00_Input_data/reference_database'
### Prefix of the reference database files must not contain "." or "_"
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//g' | awk -F/ '{print $NF}' | uniq`
# Path to intermediate and final output directories
main_dir='05_error/Outputs/01_swarm/main/'
fin_dir='05_error/Outputs/01_swarm/final/'


################################################################################################

# Merging of forward et reverse reads
$illuminapairedend -r $R2_fastq $R1_fastq > $main_dir/"$pref".fastq
# Removal of non-aligned reads
$obigrep -p 'mode!="joined"' $main_dir/"$pref".fastq > $main_dir/"$pref".ali.fastq
# Assign each sequence to its sample
$ngsfilter -t $sample_description_file -u $main_dir/"$pref"_unidentified.fastq $main_dir/"$pref".ali.fastq --fasta-output > $main_dir/"$pref".ali.assigned.fasta
# Split the file in smaller sample files
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta $main_dir/"$pref".ali.assigned.fasta.gz

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
  sample_sh="${sample/.fasta/_cmd.sh}"
  echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
  # Dereplication of reads in unique sequences
  dereplicated_sample="${sample/.fasta/.uniq.fasta}"
  echo "$obiuniq -m sample "$sample" > "$dereplicated_sample > $sample_sh;
  # Keep only sequences longer than 20pb with no ambiguous bases
  good_sequence_sample="${dereplicated_sample/.fasta/.l20.fasta}"
  echo "$obigrep -s '^[ACGT]+$' -l 20 "$dereplicated_sample" > "$good_sequence_sample >> $sample_sh
  # Format fasta file to process sequence with swarm
  formated_sequence_sample="${good_sequence_sample/.fasta/.formated.fasta}"
  echo "$obiannotate -R 'count:size'  "$good_sequence_sample" | obiannotate -k size -k merged_sample > "$formated_sequence_sample >> $sample_sh
  # Removal of PCR and sequencing errors (variants) with swarm
  clean_sequence_sample="${formated_sequence_sample/.fasta/.clean.fasta}"
  echo " /usr/bin/time $swarm -z -f -w "$clean_sequence_sample" "$formated_sequence_sample >> $sample_sh
  # Format swarm fasta file to continue the pipeline process
  formated_clean_sequence_sample="${clean_sequence_sample/.fasta/.formated.fasta}"
  echo "$obiannotate -R 'size:count' "$clean_sequence_sample" > "$formated_clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
# Concatenation of all samples in one file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.l20.formated.clean.formated.fasta > $all_sample_sequences_clean
# Dereplication in unique sequences
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
$obiuniq -m sample $all_sample_sequences_clean > $all_sample_sequences_uniq
# Taxonomic assignation
all_sample_sequences_tag="${all_sample_sequences_uniq/.fasta/.tag.fasta}"
$ecotag -d $base_dir/"${base_pref}" -R $refdb_file $all_sample_sequences_uniq > $all_sample_sequences_tag
# Removal of useless attributes in sequences headers
all_sample_sequences_ann="${all_sample_sequences_tag/.fasta/.ann.fasta}"
$obiannotate --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
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
# Create final table
$obitab -o $all_sample_sequences_sort > $fin_dir/"$step".csv

gzip $main_dir/*