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

##  read vsearch assignement --blast6out
dicOfseq={}
with open(vsearchFile,'r') as readFile:
    pair=0;
    for ligne in readFile.readlines():
        if pair !=1:
            ligneSplit=ligne.split(";")
            ## edit ID tax best match
            bestmatch=ligneSplit[3].replace(" ","")
            ligneSplit[3]="best_match={'db_teleo1': '"+bestmatch+"'}"
            ## edit best_identity
            best_identity=ligneSplit[-1].split()[0]
            ligneSplit[-1]="{'db_teleo1': "+best_identity+"}"
            lignePropre="; ".join(str(elem) for elem in ligneSplit[0:])+";"
            dicOfseq[ligneSplit[0].split(" ")[0]]=lignePropre
            pair=1
        else:
            pair=0

print dicOfseq

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
