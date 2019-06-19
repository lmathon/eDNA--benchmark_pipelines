 [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2878)

# benchmark_pipelines

Assessment of metabarcoding data pipelines :
The goal of this study is to compare the performance of several bioinformatic tools to analyze fish eDNA data, in order to build the most efficient pipeline.
Seven steps of the analysis were identified, from the assembly of paired-end reads to the taxonomic assignation.
For each of these steps, different programs to compare were identified, as shown below :

![pipeline_schema](schema_protocole.PNG)

For each step, all the programs are compared, while the start and the end of the pipeline are standardized with a reference pipeline ([Script_obitools_reference.sh](Script_obitools_reference.sh)). This pipeline is based on [ObiTools](https://git.metabarcoding.org/obitools/obitools/wikis/home), a set of python programs designed to analyse Next Generation Sequencer outputs (illumina) in the context of DNA Metabarcoding.

The optimal pipeline obtained will be again compared to existant complete pipelines (QIIME2, Mothur, BARQUE, DADA2 and SLIM).

# Dependencies

* [BASH](https://www.gnu.org/software/bash/)
	- [openssl](https://www.openssl.org/)
	- [wget](https://www.gnu.org/software/wget/)
	- [curl](https://curl.haxx.se/)

* [Singularity](https://www.sylabs.io/docs/)

# Installation

## Install from source code

To install all the programs used in this study, please follow the instructions on their installation pages : [ObiTools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools), [VSEARCH](https://github.com/torognes/vsearch), [USEARCH](https://drive5.com/usearch/download.html), [PEAR](http://www.exelixis-lab.org/web/software/pear), [FLASH](https://sourceforge.net/projects/flashpage/files), [PANDAseq](https://github.com/neufeld/pandaseq), [CASPER](http://best.snu.ac.kr/casper/index.php?name=manual), [fastq-join](https://github.com/brwnj/fastq-join), [fatsq-multx](https://github.com/brwnj/fastq-multx), [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html), [Flexbar](https://github.com/seqan/flexbar), [TagCleaner](https://sourceforge.net/projects/tagcleaner/files), [deML](https://github.com/grenaud/deml), [Tally](https://www.ebi.ac.uk/research/enright/software/kraken), [DADA2](https://benjjneb.github.io/dada2/dada-installation.html), [Prinseq](https://sourceforge.net/projects/prinseq/files/), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic), [SWARM](https://github.com/torognes/swarm), [Perseus](https://code.google.com/archive/p/ampliconnoise/), [CATCh](https://github.com/M-Mysara/CATCh), and [PROTAX](https://www.helsinki.fi/en/researchgroups/statistical-ecology/software#section-49869).

The installation guidelines for the complete pipelines can be found here : [QIIME2](https://docs.qiime2.org/2019.4/install), [MOTHUR](https://github.com/mothur/mothur), [BARQUE](https://github.com/enormandeau/barque) and [SLIM](https://github.com/yoann-dufresne/SLIM).

## Singularity containers

Other possibility : all the programs have been installed in a [Singularity container](https://www.sylabs.io/docs/), which you can access here :

First you need to install [Singularity](https://github.com/sylabs/singularity/blob/master/INSTALL.md).

We provide ready to run versions of [Singularity containers](https://www.sylabs.io/).

Our complete collection of singularity recipes is available [here](https://github.com/Grelot/bioinfo_singularity_recipes).

To download [ObiTools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools)'s container :
```
singularity pull --name 99_utils/images/obitools.img shub://Grelot/bioinfo_singularity_recipes:obitools
```

To download container with vsearch, usearch, PEAR, FLASh, CASPER, cutadapt, fastq-multx, fastq-join, TAGcleaner, Reaper, SWARM, Flexbar, deML :
```
singularity pull --name ednatools.img shub://Grelot/bioinfo_singularity_recipes:ednatools
```


# Data requirement

## Reference database

Download, uncrypt and unzip `reference_database` from our [MEGA cloud](https://mega.nz/). The `reference_database` folder will be created at [00_Input_data/reference_database](00_Input_data) :
```
bash 99_utils/mega_download/download_input_data.sh
```

## Input data

The dataset used for this study, containing forward and reverse reads from 12S mitochondrial gene fragment of fish, has been simulated with [Grinder](https://sourceforge.net/projects/biogrinder/). For the full simulation protocole, please read the [README](grinder_simulations/README.md) in [grinder_simulations](grinder_simulations).

The dataset is stored on MEGA. To dowload, uncrypt and unzip `forward_reverse_reads` from our [MEGA cloud](https://mega.nz/), run :
```
--
```

The `forward_reverse_reads` folder will be created at [00_Input_data/forward_reverse_reads](00_Input_data).

## Sample description file

The sample description file, containing the primers and the tags associated to each sample, is available [here](00_Input_data/sample_description_file.txt).


# Performance measures

To assess the efficiency of each program, we measure the time, % of CPU used and the memory used (among other metrics).

Each time you test a different program for a given analysis step, you can record the memory usage of this command by running [record_memory_usage.sh](99_utils/record_memory_usage.sh) :

```
bash 99_utils/record_memory_usage.sh PID NUMBER_ITER TABLE_MEM
```


# Analysis steps

For simplicity, each command of the pipeline are separated in folders corresponding to the steps, in [Scripts](Scripts). To build the pipeline, pick the scripts corresponding to the commands you want in the rigth folder. 

## 1 - Merging paired-end reads

The first step consists in assembling forward and reverse reads of the sequences. We tested several assemblers, with no specific parameters.

[01_merging](Scripts/01_merging) contains the scripts to run each of these programs.

## 2 - Demultiplexing

Once the reads assembled, the primers are removed (max. 2 mismatches allowed by primers). The tags are also removed (no mismatch allowed) and each read is assigned to the sample it comes from.

[02_demultiplex](Scripts/02_demultiplex) contains the scripts to run each programs used at this step.

## 3 - Dereplicating

The reads are then dereplicated: identical reads are gathered in a unique read and the count is saved.

All the scripts to run the different programs are in [03_dereplication](Scripts/03_dereplication).

## 4 - Quality filtering

Reads are then checked for their quality : sequences longer than 20bp and with no ambiguous bases are kept.

The scripts to run the different programs are in [04_quality-filter](Scripts/04_quality-filter)

## 5 - PCR / Sequencing error removal

Each program or pipeline offers different tools to remove PCR or sequencing errors. For ObiTools, the program obiclean keeps only the head sequences for each PCR replicate.

The scripts for the different programs are in [05_error-removal](Scripts/05_error-removal).

## 6 - Chimera removal

Several programs are specialized in identifying and removing chimeras. 

[06_chimera-removal](Scripts/06_chimera-removal) contains the scripts to run these programs.

## 7 - Taxonomic assignation

The last step of the analysis is to assign every sequence to a taxa. In our case, we use a homemade reference database. To be assigned at the species level, the query sequence must be at least similar at 98% to the reference sequence.

[07_taxonomic-assignation](Scripts/07_taxonomic-assignation) contains the scripts to run the different assigning programs.


The script [last_step](Scripts/last_step.sh) organize the data in a table, and can be run at the end of the pipeline.

# Outputs

## Taxa/sample tables

The outputs of each pipeline tested can be found in [Outputs](Outputs), in the folder corresponding to the step tested, under the name of the program tested.
For example, to find the results of the pipeline testing the program flash for merging reads : Outputs/01_merging/merging_flash.csv

## Time and memory reports

Time and memory reports for each program compared are stored in [Reports](Reports), also sorted by step.