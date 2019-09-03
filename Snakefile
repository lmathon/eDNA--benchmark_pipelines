configfile: "98_infos/obitools_reference_config.yaml"


rule all:
	input:
		expand(config["output_data"]["demultiplex"]+"/{run}.fastq", run=config["filename"]),
		expand('log/illuminapairedend/{run}.log',run=RUNS)



include: "rules/obitools_reference.smk"
