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

python3 sintax2obitab.py -s 06_assignation/Outputs/02_sintax/main/grinder_teleo1_all_sample_clean.uniq.ann.sort.uppercase.tag.fasta -o output.tsv

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


listOfLignes = []
listOfIds = []

with open(sintaxFile,'r') as readFile:
    for ligne in readFile.readlines():     
        ligneSplit=ligne.split(";")
        thisLigne= Ligne("NA","NA","NA","NA","NA","NA","NA")
        thisLigne.id_ligne=ligneSplit[0].split(" ")[0]        
        if thisLigne.id_ligne not in listOfIds:
            listOfIds.append(thisLigne.id_ligne)
            thisLigne.count=ligneSplit[0].split(" ")[1].split("=")[1]        
            for elem in ligneSplit[1:]:                
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
                    ## check if taxon field is not empty
                    if 'f' in elem:
                        elemFormatAll=elem.replace("\t","").replace("\n","").lstrip().split("+")                       
                        if len(elemFormatAll) > 1:
                            if len(elemFormatAll[1]) > 1:
                                elemFormat = elemFormatAll[1].split(",")                            
                                for elemf in elemFormat:
                                    if elemf[0] =='f':
                                        thisLigne.family_name=elemf.split(':')[1].split('(')[0] #family
                                    elif elemf[0] == 'g':
                                        thisLigne.genus_name=elemf.split(':')[1].split('(')[0] #genus
                                    elif elemf[0] == 's':
                                        thisLigne.species_name=elemf.split(':')[1].split('(')[0] #species
                                    else:
                                        print("Error no taxon information at line %s.", thisLigne.id_ligne)
            listOfLignes.append(thisLigne)
        else:
            print("non")
            for i in range(len(listOfLignes)):
                if listOfLignes[i].id_ligne == thisLigne.id_ligne:
                    thisLigne.definition="r"
                    for elem in ligneSplit[1:]:
                        if "=" not in elem:
                            elemFormatAll=elem.replace("\t","").replace("\n","").lstrip().split("+")
                            if len(elemFormatAll) > 1:
                                if len(elemFormatAll[1]) > 1:
                                    elemFormat = elemFormatAll[1].split(",")       
                                    for elemf in elemFormat:
                                        if elemf[0] =='f':
                                            listOfLignes[i].family_name=thisLigne.family_name+","+elemf.split(':')[1].split('(')[0] #family
                                        elif elemf[0] == 'g':
                                            listOfLignes[i].genus_name=thisLigne.genus_name+","+elemf.split(':')[1].split('(')[0] #genus
                                        else:
                                            listOfLignes[i].species_name=thisLigne.species_name+","+elemf.split(':')[1].split('(')[0] #species


## get all keys of sample from merged_sample dic
all_keys = list(set().union(*(d.merged_sample.keys() for d in listOfLignes)))


## add zero value to each missing keys into merged_sample dic for each line
for ligne in listOfLignes:
    ligne.fill_missing_keys(all_keys) 


## write the final table
sourceFile= open(outputFile,"w+")
### write header
print("id","definition","count","family_name","genus_name","species_name", *[key for key in all_keys],sep="\t",file=sourceFile)
### write line content
for ligne in listOfLignes:
    print(ligne.id_ligne,ligne.definition,ligne.count,ligne.family_name,ligne.genus_name,ligne.species_name,*[ligne.merged_sample[key] for key in all_keys ],sep="\t",file=sourceFile)

sourceFile.close()
