R1_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R1.fastq.gz'
R2_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R2.fastq.gz'
output= #output dir and output pref

pear -f $R1_fastq -r $R2_fastq -v 10 -c 0 -n 0 -o $output

# output with *.assembled.fastq