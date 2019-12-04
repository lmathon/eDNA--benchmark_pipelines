source 98_infos/config.sh
flexbar=${SINGULARITY_EXEC_CMD}" "${EDNATOOLS_SIMG}" flexbar"

## Prefix for all generated files

pref="grinder_teleo1"
## Path to the directory containing forward and reverse reads
R1_fastq="${DATA_PATH}"/"$pref"/"$pref"_R1.fastq
R2_fastq="${DATA_PATH}"/"$pref"/"$pref"_R2.fastq
## path to 'tags.fasta'
Tags=`pwd`"/02_demultiplex/Tags.fasta"

Tags_F=`pwd`"/02_demultiplex/Tags_F.fasta"
Tags_R=`pwd`"/02_demultiplex/Tags_R.fasta"
Primer_F=`pwd`"/02_demultiplex/Primer_F.fasta"
Primer_R=`pwd`"/02_demultiplex/Primer_R.fasta"

## flexbar - flexible barcode and adapter removal


## demultiplexing
$flexbar -r $R1_fastq -p $R2_fastq --barcodes $Tags -t 02_demultiplex/test/02_flexbar/"$pref" -n 4





######## cimetiere des commandes non fonctionnelles 
#$flexbar -r $R1_fastq -p $R2_fastq -b $Tags -b2 $Tags -a $Primer_F -a2 $Primer_R -at 0.1 -t 02_demultiplex/test/02_flexbar/"$pref" -n 4
#$flexbar -r $R1_fastq -p $R2_fastq --barcodes $Tags_F --barcodes2 $Tags_R -bt 0 -t 02_demultiplex/test/02_flexbar/"$pref" -n 4
#$flexbar -r $R1_fastq -p $R2_fastq --barcodes2 $Tags -bt 0 -t 02_demultiplex/test/02_flexbar/"$pref" -n 4
#$flexbar -p $R2_fastq --barcodes2 $Tags -t 02_demultiplex/test2/02_flexbar/"$pref" -n 4

## remove adapters
$flexbar -r $R1_fastq -p $R2_fastq --adapters $Primer_F --adapters2 $Primer_R -bt 0 -t 02_demultiplex/test/02_flexbar/"$pref" -n 4

## clean
rm 02_demultiplex/test/02_flexbar/*