input=
output=
unidentified=
sample_description_file='00_Input_data/sample_description_file.txt'

ngsfilter -t $sample_description_file -u $unidentified $input --fasta-output > $output