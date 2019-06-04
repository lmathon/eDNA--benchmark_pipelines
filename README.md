# benchmark_pipelines

Assessment of metabarcoding data pipelines :
The goal of this study is to compare the performance of several bioinformatic tools to analyse fish eDNA data, in order to build the most efficient pipeline.
Seven steps of the analysis were identified, from the assembly of paired-end reads to the taxonomic assignation.
For each of these steps, different programs to compare were identified, as shown below :

![pipeline_schema](schema_protocole.PNG)

For each step, all the programs are compared, while the start and the end of the pipeline are standardized with a reference pipeline ([Script_obitools_reference.sh](Script_obitools_reference.sh)). This pipeline is based on [OBItools](https://git.metabarcoding.org/obitools/obitools/wikis/home), a set of python programs designed to analyse Next Generation Sequencer outputs (illumina) in the context of DNA Metabarcoding.

The optimal pipeline obtained will be again compared to existant complete pipelines (QIIME2, Mothur, BARQUE and SLIM).

# Dependencies

# Installation

To install all the programs used in this stidy, please follow the instructions on their installation pages : [Obitools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools), [VSEARCH](https://github.com/torognes/vsearch), [USEARCH](https://drive5.com/usearch/download.html), [PEAR](http://www.exelixis-lab.org/web/software/pear), [FLASH](https://sourceforge.net/projects/flashpage/files), [PANDAseq](https://github.com/neufeld/pandaseq), [CASPER](http://best.snu.ac.kr/casper/index.php?name=manual), [EA-util](https://github.com/ExpressionAnalysis/ea-utils/tree/master), [QIIME](http://qiime.org/install/install.html), [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html), [Flexbar](https://github.com/seqan/flexbar), [TagCleaner](https://sourceforge.net/projects/tagcleaner/files), [deML](https://github.com/grenaud/deml), [Tally](https://www.ebi.ac.uk/research/enright/software/kraken), [DADA2](https://benjjneb.github.io/dada2/dada-installation.html), [Prinseq](https://sourceforge.net/projects/prinseq/files/), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic), [SWARM](https://github.com/torognes/swarm), [Perseus](https://code.google.com/archive/p/ampliconnoise/), [CATCh](https://github.com/M-Mysara/CATCh), and [PROTAX](https://www.helsinki.fi/en/researchgroups/statistical-ecology/software#section-49869).
The installation guidelines for the complete pipelines can be found here : [QIIME2](https://docs.qiime2.org/2019.4/install), [MOTHUR](https://mothur.org/wiki/Download_mothur), [BARQUE](https://github.com/enormandeau/barque) and [SLIM](https://github.com/yoann-dufresne/SLIM).

Other possibility : all the programs have been installed in a Singularity container, which you can access here :
But first you need to install [Singularity](https://github.com/sylabs/singularity/blob/master/INSTALL.md).

# Analysis steps

## 1 - Merging paired-end reads
## 2 - Demultiplexing
## 3 - Dereplicating
## 4 - Quality filtering
## 5 - PCR / Sequencing error removal
## 6 - Chimera removal
## 7 - Taxonomic assignation


