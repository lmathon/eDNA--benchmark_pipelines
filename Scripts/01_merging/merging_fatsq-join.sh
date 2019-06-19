R1_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R1.fastq.gz'
R2_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R2.fastq.gz'
output_dir= #output_dir and output_pref

fastqjoin $R1_fastq $R2_fastq -m 10 -p 25 -o "$output"_%.fastq

# output avec *_join.fastq