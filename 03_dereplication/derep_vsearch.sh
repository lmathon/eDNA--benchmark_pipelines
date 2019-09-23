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
vsearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" vsearch"
container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"


## Prefix for all generated files
pref="grinder_teleo1"
## Prefix of the final table, including the step and the program tested (ie: merging_obitools)
step="demultiplex_cutadapt"
## Path to the directory containing forward and reverse reads
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
main_dir=`pwd`"/03_dereplication/Outputs/01_vsearch/main"
fin_dir=`pwd`"/03_dereplication/Outputs/01_vsearch/final"
## sample description file
sample_description_file='00_Input_data/sample_description_file.txt'


################################################################################################

## forward and reverse reads assembly
assembly=${main_dir}"/"${pref}".fastq"
$illuminapairedend -r ${R2_fastq} ${R1_fastq} > ${assembly}
## Remove non-aligned reads
assembly_ali="${assembly/.fastq/.ali.fastq}"
$obigrep -p 'mode!="joined"' ${main_dir}"/"${pref}".fastq" > ${assembly_ali}
## Assign each sequence to a sample
identified="${assembly_ali/.ali.fastq/.ali.assigned.fasta}"
unidentified="${assembly_ali/.ali.fastq/_unidentified.fastq}"
$ngsfilter -t ${sample_description_file} -u ${unidentified} ${assembly_ali} --fasta-output > ${identified}
## Séparation du fichier global en un fichier par échantillon 
$obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta ${identified}

################################################################################################

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
## Dereplicate reads in unique sequences
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo $vsearch" --derep_fulllength "$sample" --sizeout --fasta_width 0 --notrunclabels --output "$dereplicated_sample > $sample_sh
# Formatage des sorties vsearch en obifasta
formated_sample="${dereplicated_sample/.fasta/.formated.fasta}"
echo $container_python2" 03_dereplication/vsearch_to_obifasta.py -f "$dereplicated_sample" -o "$formated_sample >> $sample_sh
## Keep only sequences longer than 20pb with no N bases
good_sequence_sample="${formated_sample/.fasta/.l20.fasta}"
echo $obigrep" -s '^[ACGT]+$' -l 20 "$formated_sample" > "$good_sequence_sample >> $sample_sh
## Removal of PCR and sequencing errors (variants)
clean_sequence_sample="${good_sequence_sample/.fasta/.r005.clean.fasta}"
echo $obiclean" -r 0.05 -H "$good_sequence_sample" > "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
## Concatenate all files into one main file
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.formated.l20.r005.clean.fasta > $all_sample_sequences_clean
## Dereplicate in unique sequences
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
$vsearch --derep_fulllength $all_sample_sequences_clean --sizeout --uc $main_dir/info_seq --fasta_width 0 --notrunclabels --output $all_sample_sequences_uniq
# formatage des sorties vsearch pour obitools
all_sample_sequences_uniq_formated="${all_sample_sequences_uniq/.fasta/.formated.fasta}"
$container_python2 03_dereplication/allvsearch_into_obifasta.py -i $main_dir/info_seq -f $all_sample_sequences_uniq -o $all_sample_sequences_uniq_formated
# Assignation taxonomique
all_sample_sequences_tag="${all_sample_sequences_uniq_formated/.fasta/.tag.fasta}"
$ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $all_sample_sequences_uniq_formated > $all_sample_sequences_tag
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

#gzip $main_dir/*




#derep="/share/reservebenefit/working/pierre/eDNA--benchmark_pipelines/03_dereplication/Outputs/01_vsearch/main/grinder_teleo1_all_sample_clean.uniq.fasta"
#dinfo="/share/reservebenefit/working/pierre/eDNA--benchmark_pipelines/03_dereplication/Outputs/01_vsearch/main/info_seq"

#$container_python2 03_dereplication/allvsearch_into_obifasta.py -i $dinfo -f $derep -o test
