source 98_infos/config.sh
flexbar=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" flexbar"

## Prefix for all generated files

pref="grinder_teleo1"
## Path to the directory containing forward and reverse reads
R1_fastq="${DATA_PATH}"/"$pref"/"$pref"_R1.fastq
R2_fastq="${DATA_PATH}"/"$pref"/"$pref"_R2.fastq
## path to 'tags.fasta'
Tags=`pwd`"/02_demultiplex/Tags.fasta"

$flexbar -r $R1_fastq -p $R2_fastq --barcodes2 $Tags -bt 0 -t 02_demultiplex/test/02_flexbar/"$pref" -n 4


$flexbar -r $R1_fastq --barcodes $Tags -bt 0 -t 02_demultiplex/test/02_flexbar/"$pref"


## clean
rm 02_demultiplex/test/02_flexbar/*