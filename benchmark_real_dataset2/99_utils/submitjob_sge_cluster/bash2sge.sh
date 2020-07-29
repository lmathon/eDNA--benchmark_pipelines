#!/bin/bash
##=============================================================================
## HEADER
##=============================================================================
## USAGE
##    subbash monscript.sh
##
## DESCRIPTION
##    write a SGE script to submit a given bash script with `qsub`
##    This SGE file can be submitted with `qsub` on MBB SGE cluster
##
##    We use a template file `template_qsub.sh` to write a base of SGE script
##    with default parameters and values. We add bash command running the bash
##    script defined as first argument `monscript.sh`
##
##    SGE script will be written into folder 99_utils/submitjob_sge_cluster/qsub_scripts/
##    SGE script name contains name of users, name of bash script to submit
##    and date of creation.
##
##    By default this SGE script will submit the job to `cemeb20.q` queue 
##    with multithread60 16 cores
##    SGE's ouputs or error are defined to be written into folder
##    99_utils/submitjob_sge_cluster/qsub_outputs/
##    
##
##
##=============================================================================
## END_OF_HEADER
##=============================================================================

template_script="99_utils/submitjob_sge_cluster/template_qsub.sh"
qsub_script="99_utils/submitjob_sge_cluster/qsub_scripts/"
run_output="99_utils/submitjob_sge_cluster/qsub_outputs/"
nom_bash_script=`echo $1 | cut -d "." -f 1`
location_bash_script=`pwd`
nom_qsub_script=`echo $nom_bash_script"_"$USER"_"$(date +%Y%m%d-%H%M%S) | sed -e "s|/|_|g"`
sed -e "s/JOB_NAME/$nom_qsub_script/g" \
-e "s|OUT_FILE|$run_output/$nom_qsub_script\.out|g" \
-e "s|ERR_FILE|$run_output/$nom_qsub_script\.err|g" \
$template_script  > $qsub_script/temp_script.sh
echo -e "cd "$(pwd) >> $qsub_script/temp_script.sh
echo -e "bash $location_bash_script/$nom_bash_script.sh\n" >> $qsub_script/temp_script.sh
mv $qsub_script/temp_script.sh $qsub_script/$nom_qsub_script".sh"
echo $qsub_script/$nom_qsub_script".sh is ready to be qsubmitted."

