#!/bin/bash
##Obitools

# Chemin vers r�pertoire contenant les reads forward et reverse
DATA_PATH=
# Prefixe pour tous les fichiers g�n�r�s
pref=grinder_teleo1
# Prefixe du tableau final, contenant l'�tape et le programme test� (ex: merging_obitools) 
step=
# Fichiers contenant les reads forward et reverse
R1_fastq="$DATA_PATH"/"$pref"_R1.fastq
R2_fastq="$DATA_PATH"/"$pref"_R2.fastq
# Chemin vers le fichier 'tags.txt'
sample_description_file=
# Chemin vers le fichier 'db_sim_teleo1.fasta'
refdb_dir=
# Chemin vers les fichiers 'embl' de la base de r�f�rence
base_dir=
### Les pr�fixes des fichiers de la base de ref ne doivent pas contenir "." ou "_"
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//'g | awk -F/ '{print $NF}' | uniq`
# Chemin vers les r�pertoires de sorties interm�diaires et finales
main_dir=
fin_dir=


################################################################################################

# Assemblage des reads forward et reverse
illuminapairedend -r $R2_fastq $R1_fastq > $main_dir/"$pref".fastq
# Supression des reads non align�s
obigrep -p 'mode!="joined"' $main_dir/"$pref".fastq > $main_dir/"$pref".ali.fastq
# Assignation de chaque s�quence � son �chantillon
ngsfilter -t $sample_description_file -u $main_dir/"$pref"_unidentified.fastq $main_dir/"$pref".ali.fastq --fasta-output > $main_dir/"$pref".ali.assigned.fasta
# S�paration du fichier global en un fichier par �chantillon 
obisplit -p $main_dir/"$pref"_sample_ -t sample --fasta $main_dir/"$pref".ali.assigned.fasta

all_samples_parallel_cmd_sh=$main_dir/"$pref"_sample_parallel_cmd.sh
echo "" > $all_samples_parallel_cmd_sh
for sample in `ls $main_dir/"$pref"_sample_*.fasta`;
do
sample_sh="${sample/.fasta/_cmd.sh}"
echo "bash "$sample_sh >> $all_samples_parallel_cmd_sh
# D�r�plication des reads en s�quences uniques
dereplicated_sample="${sample/.fasta/.uniq.fasta}"
echo "$obiuniq -m sample "$sample" > "$dereplicated_sample > $sample_sh;
# On garde les s�quences de plus de 20pb sans bases ambigues
good_sequence_sample="${dereplicated_sample/.fasta/.l20.fasta}"
echo "$obigrep -s '^[ACGT]+$' -l 20 "$dereplicated_sample" > "$good_sequence_sample >> $sample_sh
# Supression des erreurs de PCR et s�quen�age (variants)
clean_sequence_sample="${good_sequence_sample/.fasta/.r005.clean.fasta}"
echo "$obiclean -r 0.05 -H "$good_sequence_sample" > "$clean_sequence_sample >> $sample_sh
done
parallel < $all_samples_parallel_cmd_sh
# Concatenation de tous les �chantillons en un fichier
all_sample_sequences_clean=$main_dir/"$pref"_all_sample_clean.fasta
cat $main_dir/"$pref"_sample_*.uniq.l20.r005.clean.fasta > $all_sample_sequences_clean
# D�r�plication en s�quences uniques
all_sample_sequences_uniq="${all_sample_sequences_clean/.fasta/.uniq.fasta}"
obiuniq -m sample $all_sample_sequences_clean > $all_sample_sequences_uniq
# Assignation taxonomique
all_sample_sequences_tag="${all_sample_sequences_uniq/.fasta/.tag.fasta}"
ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $all_sample_sequences_uniq > $all_sample_sequences_tag
# Supression des attributs inutiles dans l'ent�te des s�quences
all_sample_sequences_ann="${all_sample_sequences_tag/.fasta/.ann.fasta}"
obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
 --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
 --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
 --delete-tag=obiclean_head  --delete-tag=obiclean_headcount \
 --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status \
 --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
 --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag \
 --delete-tag=forward_tag --delete-tag=forward_score --delete-tag=forward_primer --delete-tag=forward_match \
 --delete-tag=tail_quality --delete-tag=mode --delete-tag=seq_a_single $all_sample_sequences_tag > $all_sample_sequences_ann
# Tri des s�quences par 'count'
all_sample_sequences_sort="${all_sample_sequences_ann/.fasta/.sort.fasta}"
obisort -k count -r $all_sample_sequences_ann > $all_sample_sequences_sort
# Cr�ation d'un tableau final
obitab -o $all_sample_sequences_sort > $fin_dir/"$step".csv
