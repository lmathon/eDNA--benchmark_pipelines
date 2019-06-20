# Grinder simulations

The dataset used for this comparative study has been simulated with [Grinder](https://sourceforge.net/projects/biogrinder/) (Angly et al 2012), to reproduce real dataset from Illumina sequencing with the highest fidelity.

This simulator allows the user to specify the species (=sequence) list for each sample, the species relative abundance in each sample, to simulate sample replicates and to simulate PCR and sequencing errors similar to reality. 

Here we simulated 30 samples, with 12 replicates, containing between 18 and 83 species of fish. The sequences used are fragments of the mitochondrial gene 12S, amplified by the primers Teleo_F and Teleo_R used in Valentini et al 2016. 

Grinder [inputs](grinder_simulations/Inputs) :

- one fasta file per sample containing the species names and their sequences : species have been selected randomly from a dataset of 2076 sequences of Actinopteri, Chondrichthyes, Cladistia, Cyclostomata, Myxini and Sarcopterygii extracted from GenBank.

- one text file per sample containing the relative abundances of each species (identical between the 12 replicates) : abondances are different for each sample and have been chosen to represent natural samples from marine and freshwater ecosystems (some species very abundant and some very rare). In samples 10 and 26, all species have equal abundances. 

- one fasta file per sample containing the tags to be added to the sequences (12 different tags for the 12 replicates) : tags are sequences of 8 nucleotids, all tags are differing by 3 nucleotids. Tags are unique for each replicate of each samples, in order to attribute each sequence to the sample and the replicate it comes from.
