source 98_infos/config.sh

container_python2=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" python2"


### test allvsearch_into_obitfasta.py
INFO_SEQ=$(pwd)/03_dereplication/test/01_vsearch/grinder_teleo1_all_sample_clean.uniq.fasta
ALL_SAMPLE_SEQ_UNIQ=$(pwd)/03_dereplication/test/01_vsearch/grinder_teleo1_all_sample_clean.uniq.fasta
FORMATED=$(pwd)/03_dereplication/test/01_vsearch/FORMATED.fasta

$container_python2 03_dereplication/allvsearch_into_obifasta.py -i $INFO_SEQ -f $ALL_SAMPLE_SEQ_UNIQ -o $FORMATED
