#!/bin/bash
##Obitools

illuminapairedend='singularity exec /home/lmathon/obitools.img illuminapairedend'
obigrep='singularity exec /home/lmathon/obitools.img obigrep'
cutadapt='singularity exec /home/lmathon/ednatools.img cutadapt'
obisplit='singularity exec /home/lmathon/obitools.img obisplit'
obiuniq='singularity exec /home/lmathon/obitools.img obiuniq'
obiannotate='singularity exec /home/lmathon/obitools.img obiannotate'
obiclean='singularity exec /home/lmathon/obitools.img obiclean'
ecotag='singularity exec /home/lmathon/obitools.img ecotag'
obisort='singularity exec /home/lmathon/obitools.img obisort'
obitab='singularity exec /home/lmathon/obitools.img obitab'

# Path to the directory containing forward and reverse reads
DATA_PATH='/home/lmathon/Comparaison_pipelines/01_In_silico/00_Inputs'
# Prefix for all generated files
pref=grinder_teleo1
# Prefix of the final table, including the step and the program tested (ie: merging_obitools) 
step=demultiplex_cutadapt
# Files containing forward and reverse reads
R1_fastq="$DATA_PATH"/"$pref"_R1.fastq.gz
R2_fastq="$DATA_PATH"/"$pref"_R2.fastq.gz
# Path to file 'tags.fasta'
Tags_F='/home/lmathon/Comparaison_pipelines/01_In_silico/02_demultiplex/Tags_F.fasta'
Tags_R='/home/lmathon/Comparaison_pipelines/01_In_silico/02_demultiplex/Tags_R.fasta'
# Path to file 'db_sim_teleo1.fasta'
refdb_dir='/home/lmathon/Comparaison_pipelines/01_In_silico/db_sim_teleo1.fasta'
# Path to files 'embl' of the reference database
base_dir='/home/lmathon/reference_database'
### Reference database files should not contain "." or "_"
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//'g | awk -F/ '{print $NF}' | uniq`
# Path to intermediate and final input directory
main_dir='/home/lmathon/Comparaison_pipelines/01_In_silico/02_demultiplex/Outputs/01_cutadapt/main'
fin_dir='/home/lmathon/Comparaison_pipelines/01_In_silico/02_demultiplex/Outputs/01_cutadapt/final'


################################################################################################
## Assignation of each read to its sample
/usr/bin/time singularity exec /home/lmathon/ednatools.img bash -c "export LC_ALL=C.UTF-8 ; cutadapt --pair-adapters --pair-filter=both -g file:$Tags_F -G file:$Tags_R -y '; sample={name};' -e 0 -o $main_dir/R1.assigned.fastq -p $main_dir/R2.assigned.fastq --untrimmed-paired-output $main_dir/unassigned_R2.fastq --untrimmed-output $main_dir/unassigned_R1.fastq $R1_fastq $R2_fastq"
/usr/bin/time singularity exec /home/lmathon/ednatools.img bash -c "export LC_ALL=C.UTF-8 ; cutadapt --pair-adapters --pair-filter=both -g assigned=^ACACCGCCCGTCACTCT -G assigned=^CTTCCGGTACACTTACCATG -e 0.12 -o $main_dir/R1.assigned2.fastq -p $main_dir/R2.assigned2.fastq --untrimmed-paired-output $main_dir/untrimmed_R2.fastq --untrimmed-output $main_dir/untrimmed_R1.fastq $main_dir/R1.assigned.fastq $main_dir/R2.assigned.fastq"

##Format file post cutadapt for obitools
$obiannotate $main_dir/R1.assigned2.fastq -k sample > $main_dir/R1.assigned3.fastq
$obiannotate $main_dir/R2.assigned2.fastq -k sample > $main_dir/R2.assigned3.fastq
sed  -i -e "s/ sample/_sample/g" $main_dir/R1.assigned3.fastq
sed  -i -e "s/ sample/_sample/g" $main_dir/R2.assigned3.fastq

## Assembly of forward and reverse reads
/usr/bin/time $illuminapairedend -r $main_dir/R2.assigned3.fastq $main_dir/R1.assigned3.fastq > $main_dir/"$pref".assigned.fastq
# Format header for obitools
sed  -i -e "s/_CONS//g" $main_dir/"$pref".assigned.fastq
sed  -i -e "s/_sample/_CONS sample/g" $main_dir/"$pref".assigned.fastq

## Discard non aligned reads
/usr/bin/time $obigrep -p 'mode!="joined"' --fasta-output $main_dir/"$pref".assigned.fastq > $main_dir/"$pref".ali.assigned.fasta
## Separate main file in one file per sample 
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta $main_dir/"$pref".ali.assigned.fasta

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
## Dereplicate reads in unique sequences
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo "/usr/bin/time $obiuniq -m sample "$sample" > "$dereplicated_sample > $sample_sh;
## Keep only sequences longer than 20pb with no N bases
good_sequence_sample="${dereplicated_sample/.fasta/.l20.fasta}"
echo "/usr/bin/time $obigrep -s '^[ACGT]+$' -l 20 "$dereplicated_sample" > "$good_sequence_sample >> $sample_sh
## Removal of PCR and sequencing errors (variants)
clean_sequence_sample="${good_sequence_sample/.fasta/.r005.clean.fasta}"
echo "/usr/bin/time $obiclean -r 0.05 -H "$good_sequence_sample" > "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
## Concatenate all files into one main file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.l20.r005.clean.fasta > $all_sample_sequences_clean
## Dereplicate in unique sequences
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
/usr/bin/time $obiuniq -m sample $all_sample_sequences_clean > $all_sample_sequences_uniq
## Taxonomic assignation
all_sample_sequences_tag="${all_sample_sequences_uniq/.fasta/.tag.fasta}"
/usr/bin/time $ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $all_sample_sequences_uniq > $all_sample_sequences_tag
## Removal of useless attributes in header
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
## Sort sequences by 'count'
all_sample_sequences_sort="${all_sample_sequences_ann/.fasta/.sort.fasta}"
$obisort -k count -r $all_sample_sequences_ann > $all_sample_sequences_sort
## Create final tab
$obitab -o $all_sample_sequences_sort > $fin_dir/"$step".csv

gzip $main_dir/*.fasta