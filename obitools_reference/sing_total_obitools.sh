source 98_infos/config.sh

singularity exec -B .:"${SING_MNT}" $OBITOOLS_SIMG bash -c "cd $SING_MNT ; bash obitools_reference/total_obitools.sh"