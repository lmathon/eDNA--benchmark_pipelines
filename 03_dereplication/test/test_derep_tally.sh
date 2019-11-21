source 98_infos/config.sh
tally=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" tally"
ALL_SAMPLE_SEQ_UNIQ=$(pwd)/03_dereplication/test/01_vsearch/grinder_teleo1_all_sample_clean.uniq.fasta


$tally --fasta-in --fasta-out -i ${ALL_SAMPLE_SEQ_UNIQ} -o test