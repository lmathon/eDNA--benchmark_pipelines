
### Paired end alignment then keep reads with quality > 40
rule illuminapairedend:
    input:
        R1=config["input_data"]["reads"]+"/{run}_R1.fastq.gz",
        R2=config["input_data"]["reads"]+"/{run}_R2.fastq.gz"
    output:
        fq=config["output_data"]["demultiplex"]+"/{run}.fastq"
    singularity:
        config["container"]["obitools"]
    log:
        'log/illuminapairedend/{run}.log'
    params:
        s_min=config["illuminapairedend"]["s_min"]
    shell:
        '''illuminapairedend -r {input.R2} {input.R1} --score-min={params.s_min} > {output.fq} 2> {log}'''
