input=
output=
base_dir='00_Input_data/reference_database/'
base_pref=`ls $base_dir/*sdx | sed 's/_[0-9][0-9][0-9].sdx//'g | awk -F/ '{print $NF}' | uniq`
refdb_dir='00_Input_data/reference_database/db_sim_teleo1.fasta'


ecotag -d $base_dir/"${base_pref}" -R $refdb_dir $input > $output