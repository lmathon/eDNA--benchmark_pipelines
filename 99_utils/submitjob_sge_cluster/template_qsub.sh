# Job name
#$ -N job_JOB_NAME
# Using current working directory (otherwise, you will have to use '#$ wd /path/to/run')
#$ -cwd
# job time limits (h_rt is required [s_rt == software time limit / h_rt == hardware time limit])
#$ -l s_rt=999:55:00
#$ -l h_rt=920:00:00
# choose to run on a specific queue
# (qconf -sql (to list queues) qconf -sq queue_name (to print informations on this queue))
#$ -q cemeb20.q
# Get a mail when the job begins, ends or is suspended
#$ -m ebs
#$ -M me@mail.com
# Redirects the standard output to the named file.
#$ -o OUT_FILE
##$ -e ERR_FILE
# merge standard and error outputs
#$ -j y
# choose a parallel environment and run on 60 slots (use $PE_HOSTFILE)
#$ -pe multithread60 16
# Export all my environment variables into job runtime context
#$ -V
###################
