#!/bin/bash
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
## EDNAtools
cutadapt=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" cutadapt"

## Prefix of all generated files
pref="grinder_teleo1"
## Prefix of final table files
step="demultiplex_cutadapt"
## Path to forward and reverse fastq files
R1_fastq="${DATA_PATH}"/"$pref"/"$pref"_R1.fastq.gz
R2_fastq="${DATA_PATH}"/"$pref"/"$pref"_R2.fastq.gz
## path to 'tags.fasta'
Tags_F=`pwd`"/02_demultiplex/Tags_F.fasta"
Tags_R=`pwd`"/02_demultiplex/Tags_R.fasta"
## path to the file 'db_sim_teleo1.fasta'
refdb_dir=${REFDB_PATH}"/db_sim_teleo1.fasta"
## Path to embl files of the reference database
base_dir=${REFDB_PATH}
### remove '.' and  '_' from the prefix files
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//'g | awk -F/ '{print $NF}' | uniq`
## path to outputs final and temporary (main)
main_dir=`pwd`"/02_demultiplex/Outputs/01_cutadapt/main"
fin_dir=`pwd`"/02_demultiplex/Outputs/01_cutadapt/final"


################################################################################################

## assign each sequence to a sample
$cutadapt --pair-adapters --pair-filter=both -g file:$Tags_F -G file:$Tags_R \
-y sample="; sample={name};" -e 0 -o $main_dir/R1.assigned.fastq -p $main_dir/R2.assigned.fastq \
--untrimmed-paired-output $main_dir/unassigned_R2.fastq \
--untrimmed-output $main_dir/unassigned_R1.fastq \
$R1_fastq $R2_fastq

## ???
$cutadapt --pair-adapters --pair-filter=both \
-g assigned=^ACACCGCCCGTCACTCT -G assigned=^CTTCCGGTACACTTACCATG \
-e 0.12 -o $main_dir/R1.assigned2.fastq -p $main_dir/R2.assigned2.fastq \
--untrimmed-paired-output $main_dir/untrimmed_R2.fastq \
--untrimmed-output $main_dir/untrimmed_R1.fastq \
$main_dir/R1.assigned.fastq $main_dir/R2.assigned.fastq

##Format file post cutadapt for obitools
$obiannotate $main_dir/R1.assigned2.fastq -k sample > $main_dir/R1.assigned3.fastq
$obiannotate $main_dir/R2.assigned2.fastq -k sample > $main_dir/R2.assigned3.fastq
sed  -i -e "s/ sample/_sample/g" $main_dir/R1.assigned3.fastq
sed  -i -e "s/ sample/_sample/g" $main_dir/R2.assigned3.fastq


## forward and reverse reads assembly
assembly=${main_dir}"/"${pref}".assigned.fastq"
$illuminapairedend -r $main_dir/R2.assigned3.fastq $main_dir/R1.assigned3.fastq > ${assembly}

## Remove non-aligned reads
assembly_ali="${assembly/.fastq/.ali.fastq}"
$obigrep -p 'mode!="joined"' ${main_dir}"/"${pref}".fastq" > ${assembly_ali}

sed -i -e "s/sample/NN; sample/g" ${assembly_ali}
# split global file into sample files
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta ${assembly_ali}

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
# Déréplication des reads en séquences uniques
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo "/usr/bin/time $obiuniq -m sample "$sample" > "$dereplicated_sample > $sample_sh;
# On garde les séquences de plus de 20pb sans bases ambigues
good_sequence_sample="${dereplicated_sample/.fasta/.l20.fasta}"
echo "/usr/bin/time $obigrep -s '^[ACGT]+$' -l 20 "$dereplicated_sample" > "$good_sequence_sample >> $sample_sh
# Supression des erreurs de PCR et séquençage (variants)
clean_sequence_sample="${good_sequence_sample/.fasta/.r005.clean.fasta}"
echo "/usr/bin/time $obiclean -r 0.05 -H "$good_sequence_sample" > "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
# Concatenation de tous les échantillons en un fichier
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.l20.r005.clean.fasta > $all_sample_sequences_clean
# Déréplication en séquences uniques
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
/usr/bin/time $obiuniq -m sample $all_sample_sequences_clean > $all_sample_sequences_uniq
# Assignation taxonomique
all_sample_sequences_tag="${all_sample_sequences_uniq/.fasta/.tag.fasta}"
/usr/bin/time $ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $all_sample_sequences_uniq > $all_sample_sequences_tag
# Supression des attributs inutiles dans l'entête des séquences
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
# Tri des séquences par 'count'
all_sample_sequences_sort="${all_sample_sequences_ann/.fasta/.sort.fasta}"
$obisort -k count -r $all_sample_sequences_ann > $all_sample_sequences_sort
# Création d'un tableau final
$obitab -o $all_sample_sequences_sort > $fin_dir/"$step".csv

gzip $main_dir/*.fasta