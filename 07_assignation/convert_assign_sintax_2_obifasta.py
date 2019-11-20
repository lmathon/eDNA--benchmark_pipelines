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

python2 convert_assign_vsearch_2_obifasta.py -s assignsintaxfile -f obisortoutput -o output.obifasta

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
parser.add_argument("-s","--sintax_assignation",type=str)
parser.add_argument("-f","--obifasta", type=str)


#===============================================================================
#MAIN
#===============================================================================

args = parser.parse_args()
sintaxFile = args.sintax_assignation
outputFile = args.output
obitFastaFile = args.obifasta


#sintaxFile="07_assignation/test/02_sintax/assigned_sintax.csv"
#obiFastaFile="07_assignation/test/02_sintax/grinder_teleo1_all_sample_clean.uniq.ann.sort.fasta"
#outputFile="sintax_test.fasta"

dicOfseq={}

sintax_regex = re.compile('f:.*,g:.*,s:.*')


##  read vsearch assignement --blast6out
with open(sintaxFile,'r') as readFile:
    for ligne in readFile.readlines():     
        ligneSplit=ligne.split(";")
        #print ligneSplit
        infoline=""
        for elem in ligneSplit[1:]:           
            elemSplit=elem.split("=")
            if len(elemSplit) > 1:
                infoTag=elemSplit[0].replace(" ","")
                infoVal=elemSplit[1]
                if infoTag == "merged_sample":
                    infoline+=infoTag+"="+infoVal+"; "
                else:
                    infoline=infoline
            else:
            	if re.search(sintax_regex, elem) is not None:
            		raw_assign= elem.replace("\t","").replace("+","")
            		for lv in raw_assign.split(","):
            			lvSplit=lv.split("(")[0].split(":")
            			if len(lvSplit) == 2:
            				assignTag=lvSplit[0]
            				assignVal=lvSplit[1]
            				if assignTag == "f" and "family_name" not in infoline:
            					infoline+="family_name="+assignVal+"; "
            				elif assignTag == "g" and "species_name" not in infoline:
            					infoline+="genus_name="+assignVal+"; "
            				elif assignTag == "s" and "species_name" not in infoline:
            					infoline+="species_name="+assignVal+"; "
            				else:
            					infoline=infoline
        infoline=ligneSplit[0]+"; "+infoline
        #print infoline
        dicOfseq[ligneSplit[0].split(" ")[0]]=infoline



mes_records=[]
for seq_record in SeqIO.parse(obiFastaFile, "fasta", alphabet=IUPAC.unambiguous_dna):
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
