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

WORKDIR=./job_smallp_phasespace_v2_${NGtype}/job.${SGE_TASK_ID}

# create a local directory for this job only
if [ ! -e $WORKDIR ]; then
    mkdir -p $WORKDIR
else
    # this should never happen
    echo "Clean up /data/$USER/ directory on $HOSTNAME"\
        | mail -s "SGE-error" $USER
    exit 1
fi

# create a local copy of the program and start the job
OLDDIR=`pwd`

#---------COPY NECESSARY FILES OVER TO THE WORK DIRECTORY-------
cp void_numerics.x  $WORKDIR/
cp void_numerics.in $WORKDIR/
# cp mu_v_irho_folder/mu_v_*           $WORKDIR/


#-----------GO TO THE WORK DIRECTORY AND EXECUTE ----------------
cd $WORKDIR
./void_numerics.x ${SGE_TASK_ID} $NGtype


# copy all output files back to your home directory
# and clean up
# mkdir $OLDDIR/job_${JOB_ID}
# cp -R $WORKDIR/ $OLDDIR/job_${JOB_ID} && rm -r $WORKDIR

