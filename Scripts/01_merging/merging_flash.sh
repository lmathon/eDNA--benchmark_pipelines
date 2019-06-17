R1_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R1.fastq.gz'
R2_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R2.fastq.gz'
output_dir=


flash $R1_fastq $R2_fastq -m 10 -M 150 -x 0.40 -d $output_dir 