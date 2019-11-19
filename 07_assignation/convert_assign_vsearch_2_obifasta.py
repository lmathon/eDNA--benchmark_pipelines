#===============================================================================
#INFORMATIONS
#===============================================================================
"""
CEFE - EPHE - BENCHMARK eDNA  2019
mathon laetitia, guerin pierre-edouard

convert output from vsearch assignation to obifasta

description:

input:

output:


usage:

python2 convert_assign_vsearch_2_obifasta.py -a assignvsearchfile -f obisortoutput -o output.obifasta

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
#CLASS
#===============================================================================

## nada

#===============================================================================
#ARGUMENTS
#===============================================================================

parser = argparse.ArgumentParser('convert fasta to obiconvert format')
parser.add_argument("-o","--output",type=str)
parser.add_argument("-a","--vsearch_assignation",type=str)
parser.add_argument("-f","--obifasta", type=str)


#===============================================================================
#MAIN
#===============================================================================

args = parser.parse_args()
vsearchFile = args.vsearch_assignation
outputFile = args.output
obitFastaFile = args.obifasta

dicOfseq={}

##  read vsearch assignement --blast6out
with open(vsearchFile,'r') as readFile:
    for ligne in readFile.readlines():     
        ligneSplit=ligne.split(";")
        #print ligneSplit
        infoline=""
        for elem in ligneSplit[1:]:           
            elemSplit=elem.split("=")
            if len(elemSplit) > 1:
                infoTag=elemSplit[0]
                infoVal=elemSplit[1]
                if infoTag == "merged_sample":
                    infoline+=infoTag+"="+infoVal+"; "
                elif infoTag == "family_name":
                    infoline+=infoTag+"="+infoVal+"; "
                elif infoTag == "genus_name":
                    infoline+=infoTag+"="+infoVal+"; "
                elif infoTag == "species_name":
                    infoline+=infoTag+"="+infoVal+"; "
                elif infoTag == "rank":
                    infoline+=infoTag+"="+infoVal+"; "
                else:
                    infoline=infoline
        infoline=ligneSplit[0]+"; "+infoline
        #print infoline
        dicOfseq[ligneSplit[0].split(" ")[0]]=infoline



mes_records=[]
for seq_record in SeqIO.parse(obitFastaFile, "fasta", alphabet=IUPAC.unambiguous_dna):
    seq_record_DescriptionSplit=seq_record.description.split(";")
    vSeqId=str(seq_record_DescriptionSplit[0].split(" ")[0])    
    if vSeqId in dicOfseq:
        print vSeqId
        local_id=dicOfseq[vSeqId]
        local_seq=str(repr(str(seq_record.seq.lower()))).replace("'","")
        local_record=SeqRecord(Seq(local_seq,IUPAC.unambiguous_dna), id=local_id,description="")
        mes_records.append(local_record)
        ## remove the key (we need only one sequence by key: the "seed" sequence)
        del dicOfseq[vSeqId]
SeqIO.write(mes_records, outputFile, "fasta")
