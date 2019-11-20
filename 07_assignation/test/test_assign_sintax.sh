source 98_infos/config.sh
vsearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" vsearch"
container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"
obitab=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obitab"
usearch=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" usearch"


all_sample_sequences_sort="07_assignation/Outputs/01_vsearch/main/grinder_teleo1_all_sample_clean.uniq.ann.sort.fasta"


ASSIGNED=`pwd`"07_assignation/test/02_sintax/assigned_sintax.csv"

PREFINI=`pwd`"/07_assignation/test/02_sintax/PREFINI_sintax"

$usearch -sintax $all_sample_sequences_sort -db $refdb_dir -sintax_cutoff 0.8 -strand plus -tabbedout $ASSIGNED

