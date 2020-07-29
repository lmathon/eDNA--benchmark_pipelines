###############################################################################
## NOTICE
## script to record memory usage of a bash command
## input : process ID of the command, number of iteration (time)
## output : table with the memory usage by iteration
## example of usage :
## (here i want to know memory usage of 
## the background command "sleep 42" for 5 iterations )
## `sleep 42 &`
## `bash record_memory_usage.sh 22444 5 sleep_42_memory_usage.csv`
## or alternatively ($! is taking the PID of the last ran BACKGROUND process)
## `bash record_memory_usage.sh $! 5 sleep_42_memory_usage.csv`
## usage :
## bash record_memory_usage.sh PID NUMBER_ITER TABLE_MEM
## PID : process ID of the command (input)
## NUMBER_ITER : number of iterations to record (input)
## TABLE_MEM : file name to write the table of memory usage (output)

###############################################################################
## INPUT
## process ID of the command
PID=$1
## number of iteration to record (time)
NUMBER_ITER=$2
## a temporary file to store the output of `top`
temp_file=$(mktemp)
###############################################################################
## OUTPUT
TABLE_MEM=$3
###############################################################################
## MAIN
## catch top information relative to PID process
top -b -n $NUMBER_ITER | grep $PID > $temp_file
## get VIRT RES and %MEM colon
awk '{ print $5"," $6","$7}' $temp_file > $TABLE_MEM
## remove temporary file
rm $temp_file

