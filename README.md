 [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2878)

# benchmark_pipelines

Assessment of metabarcoding data pipelines :
The goal of this study is to compare the performance of several bioinformatic tools to analyze fish eDNA data, in order to build the most efficient pipeline.
Seven steps of the analysis were identified, from the assembly of paired-end reads to the taxonomic assignation.
For each of these steps, different programs to compare were identified, as shown below :

![pipeline_schema](schema_protocole.PNG)

For each step, all the programs are compared, while the start and the end of the pipeline are standardized with a reference pipeline ([obitools_reference](obitools_reference/total_obitools.sh)). This pipeline is based on [Obitools](https://git.metabarcoding.org/obitools/obitools/wikis/home), a set of python programs designed to analyse Next Generation Sequencer outputs (illumina) in the context of DNA Metabarcoding.

The optimal pipeline obtained will be again compared to existant complete pipelines (QIIME2 and BARQUE).

# Dependencies

* [BASH](https://www.gnu.org/software/bash/)
	- [openssl](https://www.openssl.org/)
	- [wget](https://www.gnu.org/software/wget/)
	- [curl](https://curl.haxx.se/)

* [Singularity](https://www.sylabs.io/docs/)

# Installation

## Install from source code

To install all the programs used in this study, please follow the instructions on their installation pages : [ObiTools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools), [VSEARCH](https://github.com/torognes/vsearch), [PEAR](http://www.exelixis-lab.org/web/software/pear), [FLASH](https://sourceforge.net/projects/flashpage/files), [CASPER](http://best.snu.ac.kr/casper/index.php?name=manual), [fastq-join](https://github.com/brwnj/fastq-join), [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html), [Prinseq](https://sourceforge.net/projects/prinseq/files/)), [SWARM](https://github.com/torognes/swarm) and [SINTAX](https://www.drive5.com/usearch/manual/cmd_sintax.html).

The installation guidelines for the complete pipelines can be found here : [QIIME2](https://docs.qiime2.org/2019.4/install) and [BARQUE](https://github.com/enormandeau/barque).

## Singularity containers

Other possibility : all the programs have been installed in a [Singularity container](https://www.sylabs.io/docs/), which you can access here :

First you need to install [Singularity](https://github.com/sylabs/singularity/blob/master/INSTALL.md).

We provide ready to run versions of [Singularity containers](https://www.sylabs.io/).

Our complete collection of singularity recipes is available [here](https://github.com/Grelot/bioinfo_singularity_recipes).

To download [ObiTools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools)'s container :
```
singularity pull --name obitools.simg shub://Grelot/bioinfo_singularity_recipes:obitools
```

To download container with vsearch, PEAR, FLASh, CASPER, cutadapt, fastq-join, Tally, Prinseq, usearch (sintax), SWARM, Flexbar :
```
singularity pull --name ednatools.simg shub://Grelot/bioinfo_singularity_recipes:ednatools
```


# Data requirement

## Input FASTQ files and Reference database

The dataset used for this study, containing forward and reverse reads from 12S mitochondrial gene fragment of fish, has been simulated with [Grinder](https://sourceforge.net/projects/biogrinder/). For the full simulation protocole, please visit https://github.com/lmathon/metabarcoding_data_simulation.

The dataset is stored on MEGA. To dowload, uncrypt and unzip `forward_reverse_reads` from our [MEGA cloud](https://mega.nz/), run :

```
bash 99_utils/mega_download/download_input_data.sh
```
The `forward_reverse_reads` folder will be created at [00_Input_data/forward_reverse_reads](00_Input_data).

The reference database needed for the taxonomic assignment step is also stored on MEGA and will be downloaded at the same time as the input FASTQ files. The `reference_database` folder will be created at [00_Input_data/reference_database](00_Input_data) :


## Sample description file

The sample description file, containing the primers and the tags associated to each sample, is available [here](00_Input_data/sample_description_file.txt).

## Abundance data 

Thanks to the simulated dataset, we know exactly the relative abundance of each species in each sample and replicate. These data can be found here [species_abundance_per_sample](grinder_simulations/Outputs/species_abundance_per_sample) and will be compared to the output of each pipeline tested to assess their efficiency.

(Note that the input FASTQ files and the abundance data will change each time you run a grinder simulation. The files given here correspond to the grinder simulation made to obtain the data for our program comparison)
  
# Performance measures

To assess the efficiency of each program, we measure the execution time (among other metrics).

Each time you test a different program for a given analysis step, you can record the time, memory usage, and CPU usage of this command by running `time` in front of the command :

```
/usr/bin/time command [-options]
```
This will give this output in the standard error file :

![image_time](image_time.PNG)

where :

%elapsed = time in hours:minutes:seconds


Other performance metrics will be calculated for each pipeline tested : sensitivity and F-measure index will be calculated from the number of true positive, false positive and false negative outputted by each pipeline. 
Relative abundances outputted by each pipeline are compared to the expected abundances ([species_abundance_per_sample](00_Input_data/species_abundance_per_sample)).

To run any script, run this command line :

```
bash 99_utils/submitjob_sge_cluster/bash2sge.sh SCRIPT.sh
qsub
```


# Analysis steps

For simplicity, each pipelines are separated in folders corresponding to the steps, in [01_merging](01_merging), [02_demultiplex](02_demultiplex), [03_dereplication](03_dereplication), [04_quality-filter](04_quality-filter), [05_error-removal](05_error-removal), [06_chimera-removal](06_chimera-removal) and [07_taxonomic-assignation](07_taxonomic-assignation). 

## 1 - Merging paired-end reads

The first step consists in assembling forward and reverse reads of the sequences. We tested several assemblers, with no specific parameters.

[01_merging](01_merging) contains the scripts to run each of these programs.

## 2 - Demultiplexing

Once the reads assembled, the primers are removed (max. 2 mismatches allowed by primers). The tags are also removed (no mismatch allowed) and each read is assigned to the sample it comes from.

[02_demultiplex](02_demultiplex) contains the scripts to run each programs used at this step.

## 3 - Dereplicating

The reads are then dereplicated: identical reads are gathered in a unique read and the count is saved.

All the scripts to run the different programs are in [03_dereplication](03_dereplication).

## 4 - Quality filtering

Reads are then checked for their quality : sequences longer than 20bp and with no ambiguous bases are kept.

The scripts to run the different programs are in [04_quality-filter](04_quality-filter)

## 5 - PCR / Sequencing error removal

Each program or pipeline offers different tools to remove PCR or sequencing errors. For ObiTools, the program obiclean keeps only the head sequences for each PCR replicate.

The scripts for the different programs are in [05_error-removal](05_error-removal).

## 6 - Chimera removal

Several programs are specialized in identifying and removing chimeras. 

[06_chimera-removal](06_chimera-removal) contains the scripts to run these programs.

## 7 - Taxonomic assignation

The last step of the analysis is to assign every sequence to a taxa. In our case, we use a homemade reference database. To be assigned at the species level, the query sequence must be at least similar at 98% to the reference sequence.

[07_taxonomic-assignation](07_taxonomic-assignation) contains the scripts to run the different assigning programs.


# Outputs

## Taxa/sample tables

The outputs of each pipeline tested can be found in the `Outputs` folder, in the folder corresponding to the step tested, under the name of the program tested.
The `main` folder contains all the intermediate files produced by the pipeline. The `final` folder contains the taxonomic table.
For example, to find the results of the pipeline testing the program flash for merging reads : 01_merging/Outputs/02_flash/final/merging_flash.csv

## Time and memory reports

Time reports for each program compared are stored [here](99_utils/submitjob_sge_cluster/qsub_outputs).

# Optimal pipeline

The optimal pipeline built with the most performant program from each step is [here](optimal_pipeline/optimal_pipeline.sh)