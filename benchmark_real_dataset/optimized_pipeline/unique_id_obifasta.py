#===============================================================================
#INFORMATIONS
#===============================================================================
"""
CEFE - EPHE - BENCHMARK eDNA  2019
mathon laetitia, guerin pierre-edouard


description:

attribute an unique ID to each sequence of a fasta file

usage:

python3 unique_id_obifasta.py input.fasta > output.fasta

"""
#===============================================================================
#MODULES
#===============================================================================
import sys


#===============================================================================
#MAIN
#===============================================================================
obifasta = sys.argv[1]


#obifasta="07_assignation/Outputs/01_vsearch/main/grinder_teleo1_all_sample_clean.uniq.ann.sort.fasta"

counterid=1
with open(obifasta,'r') as readFile:
    for ligne in readFile.readlines():       
        if ligne[0] == ">":
            ligneSplit=ligne.split(";")
            current_id=ligneSplit[0]
            new_id=">"+str(counterid)+" "+" ".join(current_id.split(" ")[1:])
            ligneSplit[0]=new_id
            print(*[elem for elem in ligneSplit],sep=";",end='')
            counterid+=1
        else:
            print(ligne,end='')
