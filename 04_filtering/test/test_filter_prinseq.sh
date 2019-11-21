source 98_infos/config.sh
prinseq=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" perl prinseq-lite-0.20.4/prinseq-lite.pl"
obiclean=${SINGULARITY_EXEC_CMD}" "${OBITOOLS_SIMG}" obiclean"


dereplicated_sample=`pwd`"/04_filtering/test/04_prinseq/grinder_teleo1_sample_S15-11.uniq.l20.fasta"
good_sequence_sample=`pwd`"/04_filtering/test/04_prinseq/OUT_PRINTSEQ"

$prinseq -fasta "${dereplicated_sample}" -min_len 20 -ns_max_n 0 -noniupac --out_format 1 -out_good "${good_sequence_sample}"

clean_sequence_sample="${good_sequence_sample}".r005.clean.fasta
$obiclean -r 0.05 -H "${good_sequence_sample}".fasta > "${clean_sequence_sample}"