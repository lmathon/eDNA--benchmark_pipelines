source 98_infos/config.sh
vsearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" vsearch"
container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"
obitab=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obitab"
usearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" usearch"

refdb_dir=`pwd`"/07_assignation/db_teleo_sintax.fasta"


all_sample_sequences_sort="07_assignation/Outputs/01_vsearch/main/grinder_teleo1_all_sample_clean.uniq.ann.sort.fasta"


ASSIGNED=`pwd`"/07_assignation/test/02_sintax/assigned_sintax.csv"


#all_sample_sequences_sort="/share/reservebenefit/working/lmathon/eDNA--benchmark_pipelines/07_assignation/Outputs/02_sintax/main/grinder_teleo1.ali.assigned.fasta"


all_sample_sequences_sort_upper=`pwd`"/07_assignation/test/02_sintax/grinder_teleo1_all_sample_clean.uniq.ann.sort.up.fasta"

awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $0}}' $all_sample_sequences_sort > $all_sample_sequences_sort_upper


PREFINI=`pwd`"/07_assignation/test/02_sintax/PREFINI_sintax"

$usearch -sintax $all_sample_sequences_sort_upper -db $refdb_dir -sintax_cutoff 0.8 -strand plus -tabbedout $ASSIGNED

$container_python2 07_assignation/convert_assign_sintax_2_obifasta.py -f $all_sample_sequences_sort -a $ASSIGNED -o $PREFINI

$obitab -o $PREFINI > sintax_test.csv