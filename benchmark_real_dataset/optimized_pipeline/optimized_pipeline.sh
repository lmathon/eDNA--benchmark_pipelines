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
obicomplement=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obicomplement"
obisort=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obisort"
vsearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" vsearch"
cutadapt=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" cutadapt"
swarm=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" swarm"
container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"

# Prefix for all generated files
pref="Carry"
# Path to forward and reverse fastq files
R1_fastq="${DATA_PATH}"/"$pref"/"$pref"_R1.fastq
R2_fastq="${DATA_PATH}"/"$pref"/"$pref"_R2.fastq
## path to 'tags.fasta'
Tags=`pwd`"/benchmark_real_dataset/optimized_pipeline/Tags.fasta"
# Path to file 'db_sim_teleo1.fasta'
refdb_dir=${REFDB_PATH}"/db_carry_vsearch_all.fasta"
# Path to all files 'embl' of the reference database
base_dir=${REFDB_PATH}
### Prefix of the ref database files must not contain "." ou "_"
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//g' | awk -F/ '{print $NF}' | uniq`
# Path to intermediate and final folders
main_dir=`pwd`"/benchmark_real_dataset/optimized_pipeline/Outputs/main"
fin_dir=`pwd`"/benchmark_real_dataset/optimized_pipeline/Outputs/final"

###################################################################################################################

## forward and reverse reads assembly
assembly=${main_dir}"/"${pref}".fasta"
/usr/bin/time $vsearch --fastq_mergepairs $R1_fastq --fast_maxdiffpct 25--reverse $R2_fastq --fastq_allowmergestagger  --fastaout ${assembly}


## assign each sequence to a sample
identified="${assembly/.fasta/.assigned.fasta}"
unidentified="${assembly/.fasta/_unidentified.fasta}"
/usr/bin/time $cutadapt -g file:$Tags -y ' sample={name};' -e 0 -j 16 -O 8 --revcomp -o ${identified} \
--untrimmed-output ${unidentified} ${assembly}
## Remove primers
trimmed="${identified/.assigned.fasta/.assigned.trimmed.fasta}"
untrimmed="${identified/.assigned.fasta/_untrimmed.fasta}"
/usr/bin/time $cutadapt -g "cttccggtacacttaccatg...agagtgacgggcggtgt" \
-e 0.12 -j 16 -O 15 --revcomp -o ${trimmed} --untrimmed-output ${untrimmed} \
${identified}


annotate="${trimmed/.fasta/.ann.fasta}"
$obiannotate -k sample --fasta-output ${trimmed} > ${annotate}
sed -i 's/  1:N/_1:N/g' ${annotate}

## Split big file into one file per sample
/usr/bin/time $obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta ${annotate}


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
echo "$container_python2 /benchmark_real_dataset/optimized_pipeline/vsearch_to_obifasta.py -f "$dereplicated_sample" -o "$formated_sample >> $sample_sh
# Keep sequences longuer than 20bp without ambiguous bases
good_sequence_sample="${formated_sample/.fasta/.l20.fasta}"
echo "/usr/bin/time $vsearch --fastx_filter "$formated_sample" --notrunclabels --threads 16 --fastq_maxns 0 --fastq_minlen 20 --fastaout "$good_sequence_sample >> $sample_sh
# Format fasta file to process sequence with swarm
formated_sequence_sample="${good_sequence_sample/.fasta/.formated.fasta}"
echo "$obiannotate -S 'size:count' "$good_sequence_sample" | $container_python2 /benchmark_real_dataset/optimized_pipeline/formate_header.py > "$formated_sequence_sample >> $sample_sh
# Removal of PCR and sequencing errors (variants) with swarm
clean_sequence_sample="${formated_sequence_sample/.fasta/.clean.fasta}"
echo " /usr/bin/time $swarm -z -f -t 16 -w "$clean_sequence_sample" "$formated_sequence_sample >> $sample_sh
echo "sed -i 's/;/; /g' "$clean_sequence_sample >> $sample_sh
echo "sed -i 's/}/}; /g' "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
# Concatenation of all samples into one file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.formated.l20.formated.clean.fasta > $all_sample_sequences_clean
# Removal of unnecessary attributes in sequence headers
all_sample_sequences_ann="${all_sample_sequences_clean/.fasta/.ann.fasta}"
$obiannotate  -k count -k size -k merged_sample $all_sample_sequences_clean > $all_sample_sequences_ann
# Sort sequences by 'count'
all_sample_sequences_sort="${all_sample_sequences_ann/.fasta/.sort.fasta}"
$obisort -k count -r $all_sample_sequences_ann > $all_sample_sequences_sort
# unique ID for each sequence
all_sample_sequences_uniqid="${all_sample_sequences_sort/.fasta/.uniqid.fasta}"
python3 /benchmark_real_dataset/optimized_pipeline/unique_id_obifasta.py $all_sample_sequences_sort > $all_sample_sequences_uniqid
# Reverse-complement all sequences
all_sample_sequence_RC="${all_sample_sequences_uniqid/.fasta/.RC.fasta}"
/usr/bin/time $obicomplement $all_sample_sequences_uniqid > $all_sample_sequence_RC



# Taxonomic assignation
all_sample_sequences_vsearch_tag="${all_sample_sequences_RC/.fasta/.tag.fasta}"
$vsearch --usearch_global $all_sample_sequences_RC --db $refdb_dir --qmask none --dbmask none --notrunclabels --id 0.98 --top_hits_only --threads 16 --fasta_width 0 --maxaccepts 20 --maxrejects 20 --minseqlength 20 --maxhits 20 --query_cov 0.6 --blast6out $all_sample_sequences_vsearch_tag --dbmatched $main_dir/db_matched.fasta --matched $main_dir/query_matched.fasta
## Create final table
### preformat
all_sample_sequences_vsearch_preformat="${all_sample_sequences_vsearch_tag/.fasta/.preformat.fasta}"
tr "\t" ";" < $all_sample_sequences_vsearch_tag > $all_sample_sequences_vsearch_preformat
python3 /benchmark_real_dataset/optimized_pipeline/vsearch2obitab.py -a $all_sample_sequences_vsearch_preformat -o $fin_dir/opt_pipeline.csv


#gzip $main_dir/*