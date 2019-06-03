#!/bin/bash
##Obitools

# Chemin vers répertoire contenant les reads forward et reverse
DATA_PATH=
# Prefixe pour tous les fichiers générés
pref=grinder_teleo1
# Prefixe du tableau final, contenant l'étape et le programme testé (ex: merging_obitools) 
step=
# Fichiers contenant les reads forward et reverse
R1_fastq="$DATA_PATH"/"$pref"_R1.fastq
R2_fastq="$DATA_PATH"/"$pref"_R2.fastq
# Chemin vers le fichier 'tags.txt'
sample_description_file=
# Chemin vers le fichier 'db_sim_teleo1.fasta'
refdb_dir=
# Chemin vers les fichiers 'embl' de la base de référence
base_dir=
### Les préfixes des fichiers de la base de ref ne doivent pas contenir "." ou "_"
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//'g | awk -F/ '{print $NF}' | uniq`
# Chemin vers les répertoires de sorties intermédiaires et finales
main_dir=
fin_dir=


################################################################################################

# Assemblage des reads forward et reverse
illuminapairedend -r $R2_fastq $R1_fastq > $main_dir/"$pref".fastq
# Supression des reads non alignés
obigrep -p 'mode!="joined"' $main_dir/"$pref".fastq > $main_dir/"$pref".ali.fastq
# Assignation de chaque séquence à son échantillon
ngsfilter -t $sample_description_file -u $main_dir/"$pref"_unidentified.fastq $main_dir/"$pref".ali.fastq --fasta-output > $main_dir/"$pref".ali.assigned.fasta
# Déréplication des reads en séquences uniques
obiuniq -m sample $main_dir/"$pref".ali.assigned.fasta > $main_dir/"$pref".ali.assigned.uniq.fasta
# On garde les séquences de plus de 20pb sans bases ambigues, et on enlève les singletons
obigrep -p 'count>=2' -s '^[ACGT]+$' -l 20 $main_dir/"$pref".ali.assigned.uniq.fasta > $main_dir/"$pref".ali.assigned.uniq.l20.fasta
# Supression des erreurs de PCR et séquençage (variants)
obiclean -r 0.05 -H $main_dir/"$pref".ali.assigned.uniq.l20.fasta > $main_dir/"$pref".ali.assigned.uniq.l20.clean.fasta
# Assignation taxonomique
ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $main_dir/"$pref".ali.assigned.uniq.l20.clean.fasta > $main_dir/"$pref".ali.assigned.uniq.l20.clean.tag.fasta
# Supression des attributs inutiles dans l'entête des séquences
obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
 --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
 --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
 --delete-tag=obiclean_head  --delete-tag=obiclean_headcount \
 --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status \
 --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
 --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag \
 --delete-tag=forward_tag --delete-tag=forward_score --delete-tag=forward_primer --delete-tag=forward_match \
 --delete-tag=tail_quality --delete-tag=mode --delete-tag=seq_a_single $main_dir/"$pref".ali.assigned.uniq.l20.clean.tag.fasta > $main_dir/"$pref".ali.assigned.uniq.l20.clean.tag.ann.fasta
 # Tri des séquences par 'count'
obisort -k count -r $main_dir/"$pref".ali.assigned.uniq.l20.clean.tag.ann.fasta > $main_dir/"$pref".ali.assigned.uniq.l20.clean.tag.ann.sort.fasta
# Création d'un tableau final
obitab -o $main_dir/"$pref".ali.assigned.uniq.l20.clean.tag.ann.sort.fasta > $fin_dir/"$step".csv