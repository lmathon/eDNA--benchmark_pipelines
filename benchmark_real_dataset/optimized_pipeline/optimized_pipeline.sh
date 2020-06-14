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
##    bash optimized_pipeline/optimized_pipeline.sh
##
## Description:
##  ..............    
##
##
##
###############################################################################
## load config global variables
source benchmark_real_dataset/98_infos/config.sh


obisplit=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obisplit"
obiannotate=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obiannotate"
obisort=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obisort"
vsearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" vsearch"
cutadapt=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" cutadapt"
flexbar=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" flexbar"
container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"

# Prefix for all generated files
pref="Banyuls"
# Path to forward and reverse fastq files
R1_fastq="${DATA_PATH}"/"$pref"/"$pref"_R1.fastq.gz
R2_fastq="${DATA_PATH}"/"$pref"/"$pref"_R2.fastq.gz
## path to 'tags.fasta'
Tags_F=`pwd`"/benchmark_real_dataset/optimized_pipeline/Tags_F.fasta"
Tags_R=`pwd`"/benchmark_real_dataset/optimized_pipeline/Tags_R.fasta"
# Path to file 'db_sim_teleo1.fasta'
refdb_dir=${REFDB_PATH}"/db_banyuls_vsearch.fasta"
# Path to all files 'embl' of the reference database
base_dir=${REFDB_PATH}
### Prefix of the ref database files must not contain "." ou "_"
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//g' | awk -F/ '{print $NF}' | uniq`
# Path to intermediate and final folders
main_dir=`pwd`"/benchmark_real_dataset/optimized_pipeline/Outputs/main"
fin_dir=`pwd`"/benchmark_real_dataset/optimized_pipeline/Outputs/final"

###################################################################################################################


## assign each sequence to a sample
$cutadapt --pair-adapters --pair-filter=both -g file:$Tags_F -G file:$Tags_R -y '; sample={name};' -e 0 -j 16 -o $main_dir/R1.assigned.fastq -p $main_dir/R2.assigned.fastq --untrimmed-paired-output $main_dir/unassigned_R2.fastq --untrimmed-output $main_dir/unassigned_R1.fastq $R1_fastq $R2_fastq

## Remove primers
$cutadapt --pair-adapters --pair-filter=both -g assigned=^ACACCGCCCGTCACTCT -G assigned=^CTTCCGGTACACTTACCATG -e 0.12 -j 16 -o $main_dir/R1.assigned2.fastq -p $main_dir/R2.assigned2.fastq --untrimmed-paired-output $main_dir/untrimmed_R2.fastq --untrimmed-output $main_dir/untrimmed_R1.fastq $main_dir/R1.assigned.fastq $main_dir/R2.assigned.fastq

##Format file post cutadapt for obitools
$obiannotate $main_dir/R1.assigned2.fastq -k sample > $main_dir/R1.assigned3.fastq
$obiannotate $main_dir/R2.assigned2.fastq -k sample > $main_dir/R2.assigned3.fastq

## forward and reverse reads assembly
assembly=${main_dir}"/"${pref}".fasta"
/usr/bin/time $vsearch --fastq_mergepairs $main_dir/R1.assigned3.fastq --reverse $main_dir/R2.assigned3.fastq --fastq_allowmergestagger  --fastaout ${assembly}

## Format file for obitools
sed  -i -e "s/_sample/ sample/g" ${assembly}

## Split big file into one file per sample
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta ${assembly}

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
# Dereplicate reads into unique sequences
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo /usr/bin/time $vsearch" --derep_fulllength "$sample" --sizeout --fasta_width 0 --notrunclabels --relabel_keep --minseqlength 20 --output "$dereplicated_sample > $sample_sh
# Formate vsearch output to obifasta
formated_sample="${dereplicated_sample/.fasta/.formated.fasta}"
echo "$container_python2 optimized_pipeline/vsearch_to_obifasta.py -f "$dereplicated_sample" -o "$formated_sample >> $sample_sh
# Keep sequences longuer than 20bp without ambiguous bases
good_sequence_sample="${formated_sample/.fasta/.l20.fasta}"
echo "/usr/bin/time $flexbar --reads "$formated_sample" --max-uncalled 0 --min-read-length 20 -n 16 -o -t "$good_sequence_sample >> $sample_sh
# Format fasta file to process sequence with vsearch
formated_sequence_sample="${good_sequence_sample/.fasta.fasta/.formated.fasta}"
echo "$obiannotate -R 'count:size'  "$good_sequence_sample".fasta | $obiannotate -k size -k merged_sample > "$formated_sequence_sample >> $sample_sh
echo "sed -i 's/; size/;size/g' "$formated_sequence_sample >> $sample_sh
# Removal of PCR and sequencing errors (variants) with swarm
clean_sequence_sample="${formated_sequence_sample/.fasta/.clean.fasta}"
echo "/usr/bin/time $vsearch --cluster_unoise "$formated_sequence_sample" --sizein --sizeout --minsize 1 --unoise_alpha 2 --notrunclabels --minseqlength 20 --threads 16 --relabel_keep --centroids  "$clean_sequence_sample >> $sample_sh
# Format vsearch fasta file to continue the pipeline process
echo "sed -i 's/;size/count/g' "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
# Concatenation of all samples into one file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.formated.l20.clean.fasta > $all_sample_sequences_clean
# Removal of unnecessary attributes in sequence headers
all_sample_sequences_ann="${all_sample_sequences_clean/.fasta/.ann.fasta}"
$obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
 --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
 --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
 --delete-tag=obiclean_head  --delete-tag=obiclean_headcount \
 --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status \
 --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
 --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag \
 --delete-tag=forward_tag --delete-tag=forward_score --delete-tag=forward_primer --delete-tag=forward_match \
 --delete-tag=tail_quality --delete-tag=mode --delete-tag=seq_a_single $all_sample_sequences_clean > $all_sample_sequences_ann
# Sort sequences by 'count'
all_sample_sequences_sort="${all_sample_sequences_ann/.fasta/.sort.fasta}"
$obisort -k count -r $all_sample_sequences_ann > $all_sample_sequences_sort
# unique ID for each sequence
all_sample_sequences_uniqid="${all_sample_sequences_sort/.fasta/.uniqid.fasta}"
python3 optimized_pipeline/unique_id_obifasta.py $all_sample_sequences_sort > $all_sample_sequences_uniqid
# Taxonomic assignation
all_sample_sequences_vsearch_tag="${all_sample_sequences_uniqid/.fasta/.tag.fasta}"
$vsearch --usearch_global $all_sample_sequences_uniqid --db $refdb_dir --qmask none --dbmask none --notrunclabels --id 0.98 --top_hits_only --threads 16 --fasta_width 0 --maxaccepts 20 --maxrejects 20 --minseqlength 20 --maxhits 20 --query_cov 0.6 --blast6out $all_sample_sequences_vsearch_tag --dbmatched $main_dir/db_matched.fasta --matched $main_dir/query_matched.fasta
## Create final table
### preformat
all_sample_sequences_vsearch_preformat="${all_sample_sequences_vsearch_tag/.fasta/.preformat.fasta}"
tr "\t" ";" < $all_sample_sequences_vsearch_tag | sed 's/ merged_sample=/; merged_sample=/g' > $all_sample_sequences_vsearch_preformat
python3 optimized_pipeline/vsearch2obitab.py -a $all_sample_sequences_vsearch_preformat -o $fin_dir/opt_pipeline.csv


#gzip $main_dir/*
