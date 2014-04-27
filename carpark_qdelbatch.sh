
# list_of_jobs=$1
# if [ $# -ne 1 ]
#  then
#    echo "Need one input argument to specify the file which lists the currently running jobs."
#    exit
# fi

qstat > temp_qsout

size=$(cat temp_qsout |wc -l)
numjobs=$(echo "$size -2" |bc -l ) 

awk '{print $1 }' temp_qsout | tail -$numjobs >temp_joblist.txt
rm temp_qsout

i=1
while [ $i -le $numjobs ]
do

jobID=$( head  -$i temp_joblist.txt |tail -1)
qdel $jobID

i=$(echo "$i +1" |bc -l)

done

rm temp_joblist.txt

