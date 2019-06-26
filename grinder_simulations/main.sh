## run grinder 30 simulations 12 replicats
singularity exec -B .:/simulations grinder.img bash -c "cd /simulations ; bash grinder_simulations/script_grinder.sh" 2> grinder_simulations/Logs/
