R1_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R1.fastq.gz'
R2_fastq='00_Input_data/forward_reverse_reads/grinder_teleo1_R2.fastq.gz'
output=

vsearch --fastq_mergepairs $R1_fastq --reverse $R2_fastq --fastq_allowmergestagger --fastqout $output