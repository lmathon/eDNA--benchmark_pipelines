#===============================================================================
#INFORMATIONS
#===============================================================================
"""
CEFE - EPHE - BENCHMARK eDNA  2019
mathon laetitia, guerin pierre-edouard
creer un fichier OBIFASTA a partir des sequences SPY.fas generes par VSEARCH

description:

input:

output:


usage:

python2 vsearch_to_obifasta.py -f input_vsearch.fasta -o output_obifasta.fasta


"""
#===============================================================================
#MODULES
#===============================================================================
import argparse
import os
import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#===============================================================================
#ARGUMENTS
#===============================================================================

parser = argparse.ArgumentParser('convert fasta to obiconvert format')
parser.add_argument("-o","--output",type=str)
parser.add_argument("-f","--vsearch_fasta",type=str)


#===============================================================================
#MAIN
#===============================================================================

args = parser.parse_args()
vsearchFile = args.vsearch_fasta
outputFile = args.output

#vsearch_fasta="Outputs/01_vsearch/main/grinder_teleo1_sample_S1-01.fasta"
mes_records=[]
for seq_record in SeqIO.parse(vsearchFile, "fasta", alphabet=IUPAC.unambiguous_dna):
    #seq_record_idSplit=seq_record.id.split(";")
    seq_record_DescriptionSplit=seq_record.description.split(";")
    local_id=""
    for elem in seq_record_DescriptionSplit:
        if "sample" in elem:
            sampleName=elem.split('=')[1]
            local_id=local_id+";"+str(elem)
        elif "size" in elem:
            sizeSeq=elem.split('=')[1]
        else:
            local_id=local_id+";"+str(elem)
    local_id=local_id[1:]+" count="+sizeSeq+"; merged_sample={'"+sampleName+"': "+sizeSeq+"};"
    local_seq=str(seq_record.seq.lower()).replace("'","")
    local_record=SeqRecord(Seq(local_seq,IUPAC.unambiguous_dna), id=local_id,description="")
    mes_records.append(local_record)
SeqIO.write(mes_records, outputFile, "fasta")
