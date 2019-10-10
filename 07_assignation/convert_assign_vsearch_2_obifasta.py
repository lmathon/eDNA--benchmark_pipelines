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

class Obinfo:
     def __init__(self,id,merged_sample,count):
         self.id=id
         self.merged_sample=merged_sample
         self.count=count

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

#vsearchFile="cattest_derep"
#outputFile="cattest.obifasta"
#infoSeqFile="info_seq"



with open(vsearchFile,'r') as readFile:
    for ligne in readFile.readlines():
        ligneSplit=ligne.split(";")
        print ligneSplit


