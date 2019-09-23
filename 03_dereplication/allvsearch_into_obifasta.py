#===============================================================================
#INFORMATIONS
#===============================================================================
"""
CEFE - EPHE - BENCHMARK eDNA  2019
mathon laetitia, guerin pierre-edouard

convert concatenated vsearch outputs into obifasta format
count merged sample unique sequences occurence

description:

input:

output:


usage:

python2 allvsearch_into_obifasta.py -i info_seq -f input_vsearch.fasta -o output_obifasta.fasta


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
parser.add_argument("-f","--vsearch_fasta",type=str)
parser.add_argument("-i","--info_seq", type=str)



#===============================================================================
#MAIN
#===============================================================================

args = parser.parse_args()
vsearchFile = args.vsearch_fasta
outputFile = args.output
infoSeqFile = args.info_seq

#vsearchFile="cattest_derep"
#outputFile="cattest.obifasta"
#infoSeqFile="info_seq"


# from info_seq, calculate total number of sequences and by sample
dicOfObinfo={}
with open(infoSeqFile,'r') as readFile:
    for ligne in readFile.readlines():
        if ligne.startswith("S"):
            idSeq=str(ligne.split("\t")[8].split(";")[0].split(" ")[0])
        if ligne.startswith("H"):
            seqSampleDic={}
            for seqSample in ligne.split("\t")[8:]:
                seqSampleSplit=seqSample.split(";")
                for elem in seqSampleSplit:
                    if "merged_sample" in elem:
                        myCode=elem.replace(" ","")
                        exec(myCode)
                        seqSampleKey=merged_sample.items()[0][0]
                        seqSampleCount=merged_sample.items()[0][1]
                        seqSampleDic[seqSampleKey]=seqSampleCount
            seqMergedSample=str(seqSampleDic)
            seqCountSum=0
            for i in seqSampleDic.items():
                iCount=int(i[1])
                seqCountSum+=iCount
            dicOfObinfo[idSeq] = Obinfo(idSeq,seqMergedSample,seqCountSum)


# write again vsearch derep fasta files with new description including number of sequences and by sample for each unique sequence
mes_records=[]
for seq_record in SeqIO.parse(vsearchFile, "fasta", alphabet=IUPAC.unambiguous_dna):
    seq_record_DescriptionSplit=seq_record.description.split(";")
    #print seq_record_DescriptionSplit    
    vSeqId=str(seq_record_DescriptionSplit[0].split(" ")[0])
    if vSeqId in dicOfObinfo:
        diSeqId= dicOfObinfo[vSeqId].id
        diSeqSamples= dicOfObinfo[vSeqId].merged_sample
        diSeqCount= str(dicOfObinfo[vSeqId].count)
        local_id=diSeqId+"; count="+diSeqCount+"; merged_sample="+diSeqSamples
        print local_id
        local_seq=str(repr(str(seq_record.seq.lower()))).replace("'","")
        local_record=SeqRecord(Seq(local_seq,IUPAC.unambiguous_dna), id=local_id,description="")
        mes_records.append(local_record)
SeqIO.write(mes_records, outputFile, "fasta")



