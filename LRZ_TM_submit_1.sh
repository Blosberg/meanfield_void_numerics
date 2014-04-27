#!/bin/bash
#$ -cwd 
#$ -l h_vmem=100M,s_rt=1:0:0
#$ -m n
#$ -t 1-90
#$ -N maple_mesh
#
# The lines above are instructions to SGE.
# 
# Just copy this sample script to a directory where you have write
# permission (e. g., a sub-directory of your home-directory).
# Then edit it to suit your needs.
#
# Start batch jobs via:
# qsub ./myjob.sh
# --------------------------------------------------------------------

# create a local directory for this job only
if [ ! -e  job$JOB_ID.$SGE_TASK_ID/ ]; then
     mkdir job$JOB_ID.$SGE_TASK_ID/
else
# case it already exists (highly unusual)
echo "Clean up local directory on $HOSTNAME"\
| mail -s "SGE-error" $USER
exit 1
fi

OLDDIR=`pwd`

#cp $PROGRAM $OPT_TMP/job$JOB_ID.$SGE_TASK_ID/$PROG

cp  $OLDDIR/TM_scan_${SGE_TASK_ID}.mpl job$JOB_ID.$SGE_TASK_ID/
#OPTIONS=$(cat $OPTIONFILE | head -n $SGE_TASK_ID | tail -n 1)
cd job$JOB_ID.$SGE_TASK_ID/

maple  TM_scan_${SGE_TASK_ID}.mpl | grep -v "bytes used"  > TM_scan_1.stdout 2> TM_scan_1.stderr

# copy all output files back to your home directory
# and clean up
#cd /$OPT_TMP

#mv $OPT_TMP/job$JOB_ID.$SGE_TASK_ID/ $OLDDIR/
#rm -r $OPT_TMP/job$JOB_ID.$SGE_TASK_ID/

