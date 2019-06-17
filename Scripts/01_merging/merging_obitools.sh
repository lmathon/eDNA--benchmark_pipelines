
R1_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R1.fastq.gz'
R2_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R2.fastq.gz'
output=

illuminapairedend -r $R2_fastq $R1_fastq | obigrep -p 'mode!="joined"' > $output