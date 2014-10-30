#!/bin/bash
#$ -cwd
#$ -m eas
#$ -l h_vmem=500M,s_rt=36:59:0
#$ -t 1-30
#
# The lines above are instructions to SGE.
# 
# Just copy this sample script to a directory where you have write
# permission (e. g., a sub-directory of your home-directory).
# Then edit it to suit your needs.
#
# Start batch jobs via:
# qsub ./myjob.sh

NGtype="LNG"


#----------------------------------------------------------------------
muN=$(head -1 array_muN_a50.txt | tail -1)
eps=$(head -${SGE_TASK_ID} array_eps_a50.txt | tail -1 )
k=50

# eps=28.0 --- this is for the k-array convergence.


#----- LOG WHEN THE JOB STARTED-------
D_before=$(date)
echo " executing job " ${JOB_ID}.${SGE_TASK_ID} " with paramaters : " $NGtype $muN $a $eps  " at " $D_before >> job_${JOB_ID}.${SGE_TASK_ID}.log


#-----------GO TO THE WORK DIRECTORY AND EXECUTE ----------------
cd $WORKDIR
./void_numerics.x $NGtype $muN $a $eps


# copy all output files back to your home directory
# and clean up
# mkdir $OLDDIR/job_${JOB_ID}
# cp -R $WORKDIR/ $OLDDIR/job_${JOB_ID} && rm -r $WORKDIR

