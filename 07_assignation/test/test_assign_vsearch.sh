source 98_infos/config.sh
vsearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" vsearch"
container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"
obitab=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obitab"


refdb_dir=`pwd`"/07_assignation/db_teleo_vsearch.fasta"
all_sample_sequences_sort=`pwd`"/07_assignation/test/01_vsearch/grinder_teleo1_all_sample_clean.uniq.ann.sort.fasta"

#all_sample_sequences_vsearch_tag=`pwd`"/07_assignation/Outputs/01_vsearch/final/assign_vsearch"



all_sample_sequences_sort="07_assignation/Outputs/01_vsearch/main/grinder_teleo1_all_sample_clean.uniq.ann.sort.fasta"


ASSIGNED=`pwd`"/07_assignation/test/01_vsearch/ASSIGNED"

PREFINI=`pwd`"/07_assignation/test/01_vsearch/PREFINI"


## test assignation vsearch
$vsearch --usearch_global $all_sample_sequences_sort --db $refdb_dir --notrunclabels --id 0.8 --fasta_width 0 --top_hits_only --blast6out $ASSIGNED

$container_python2 07_assignation/convert_assign_vsearch_2_obifasta.py -f $all_sample_sequences_sort -a $ASSIGNED -o $PREFINI

$obitab -o $PREFINI > test.csv


