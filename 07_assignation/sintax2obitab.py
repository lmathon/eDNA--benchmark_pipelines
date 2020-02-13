#===============================================================================
#INFORMATIONS
#===============================================================================
"""
CEFE - EPHE - BENCHMARK eDNA  2019
mathon laetitia, guerin pierre-edouard


description:

convert output from SINTAX assignation to obitab


input:

- output from sintax assignation command


output:

- output.tsv: a tsv file with this header:
"id definition count family_name genus_name species_name sample*"
where sample* is a list of all known samples labels
Each line is a line from grinder_teleo1_all_sample_clean.uniq.ann.sort.uppercase.tag.fasta


usage:

python3 sintax2obitab.py -s 07_assignation/test/grinder_teleo1_all_sample_clean.uniq.ann.sort.uppercase.tag.fasta -o output.tsv

"""
#===============================================================================
#MODULES
#===============================================================================
import argparse
import os



#===============================================================================
#CLASS
#===============================================================================

class Ligne:
    def __init__(self,id_ligne,merged_sample,family_name,genus_name,species_name,definition,count):
        self.id_ligne = id_ligne
        self.merged_sample= merged_sample
        self.family_name = family_name
        self.genus_name = genus_name
        self.species_name = species_name
        self.definition = definition
        self.count = count
    def fill_missing_keys(self, all_keys):
        for key in all_keys:
            if key not in self.merged_sample.keys():
                self.merged_sample[key]=0
        for key in all_keys:
            self.merged_sample[key] = self.merged_sample.pop(key)



#===============================================================================
#ARGUMENTS
#===============================================================================

parser = argparse.ArgumentParser('convert fasta to obiconvert format')
parser.add_argument("-o","--output",type=str)
parser.add_argument("-s","--sintax_assignation",type=str)



#===============================================================================
#MAIN
#===============================================================================

args = parser.parse_args()
sintaxFile = args.sintax_assignation
outputFile = args.output


#sintaxFile="07_assignation/test/grinder_teleo1_all_sample_clean.uniq.ann.sort.uppercase.tag.fasta"
#outputFile="07_assignation/test/tabfin.tsv"


listOfLignes = []

with open(sintaxFile,'r') as readFile:
    for ligne in readFile.readlines():     
        ligneSplit=ligne.split(";")
        thisLigne= Ligne("NA","NA","NA","NA","NA","NA","NA")
        thisLigne.id_ligne=ligneSplit[0].split(" ")[0]
        thisLigne.count=ligneSplit[0].split(" ")[1].split("=")[1]
        #print(ligneSplit)
        for elem in ligneSplit[1:]:
            thisLigne.id_ligne=str(ligneSplit[0])
            if "=" in elem:
                elemSplit=elem.split("=")
                if len(elemSplit) > 1:
                    infoTag=elemSplit[0].replace(" ","")
                    infoVal=elemSplit[1]
                    if infoTag == "merged_sample":
                        exec(elem.replace("\t"," ").replace(" ",""))
                        thisLigne.merged_sample=merged_sample
                    else:
                        thisLigne.definition="."
            else:
                elemFormat=elem.replace("\t","").replace("\n","").lstrip().split(",")
                for elemf in elemFormat:
                    if elemf[0] =='f':
                        thisLigne.family_name=elemf.split(':')[1].split('(')[0] #family
                    elif elemf[0] == 'g':
                        thisLigne.genus_name=elemf.split(':')[1].split('(')[0] #genus
                    else:
                        thisLigne.species_name=elemf.split(':')[1].split('(')[0] #species
        listOfLignes.append(thisLigne)

## get all keys of sample from merged_sample dic
all_keys = list(set().union(*(d.merged_sample.keys() for d in listOfLignes)))


## add zero value to each missing keys into merged_sample dic for each line
for ligne in listOfLignes:
    ligne.fill_missing_keys(all_keys) 


## write the final table
sourceFile= open(outputFile,"w+")
### write header
print("id","definition","count","family_name","genus_name","species_name",*[key for key in all_keys],sep="\t",file=sourceFile)
### write line content
for ligne in listOfLignes:
    print(ligne.id_ligne,ligne.definition,ligne.count,ligne.family_name,ligne.genus_name,ligne.species_name,*[ligne.merged_sample[key] for key in all_keys ],sep="\t",file=sourceFile)

sourceFile.close()