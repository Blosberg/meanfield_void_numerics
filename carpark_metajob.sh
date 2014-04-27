#--- carpark_metajob.sh -a script that sets up the run for the phasespace scan and then submits it.
# ---last updated on  Sun Apr 27 14:57:27 CEST 2014  by  ga79moz  at location  TUM , lagrange2

#  changes from  Sun Apr 27 14:57:27 CEST 2014 : set this file to be the main control -also determines input for void_numerics.in before copying it over to the WORKDIR

#  changes from  Sun Apr 27 14:21:11 CEST 2014 : fixed the sub_tail file with the correct modulo operator (-1)

#  changes from  Sun Apr 27 13:21:33 CEST 2014 : added the #$ -S /bin/sh command at the beginning to allow me to set variables within the qsub script.

#---------------------------------------------------------------------------------------------------------------------------

NGtype=LNG

Llim=200
t0="0.0000000001"
tf="80"
a="10"  

should_plot_Vdist_v_t="0"
should_plot_rhos="1"
should_check_neg="1"

should_import_IC="0"
should_export_IC="1"

should_peakterminate="1"
t_export="0.1"

BZcond="boltzmann_on_add"
parity_check="88855888"
opath="./"


job_sub_script=carpark_jobsub_TUM.sh
WORKDIR=./void_EQ_output_on_space/job_phasespace_peakterminated_a-${a}_output_${NGtype}/

#----------------------------------------------------------------------------------------

eps_array_file="/home/t30/ger/ga79moz/grad_research_phd/project/Nucl/void_numerics/array_eps_a"${a}".txt"

if [ -f ${eps_array_file} ];
then
   echo "reading epsilon array from file  ${eps_array_file} exists."
else
   echo "ERROR: File ${eps_array_file} does not exist. Exiting"
   exit 
fi

muN_array_file="/home/t30/ger/ga79moz/grad_research_phd/project/Nucl/void_numerics/array_muN_a"${a}".txt"

if [ -f ${muN_array_file} ];
then
   echo "reading epsilon array from file  ${muN_array_file} exists."
else
   echo "ERROR: File ${muN_array_file} does not exist. Exiting"
   exit
fi


num_eps=$(cat $eps_array_file |wc -l)
num_mus=$(cat $muN_array_file |wc -l)

num_tasks=$(( ${num_eps}*${num_mus} ))

echo "num_eps=" ${num_eps} ", num_mus=" ${num_mus}
echo "preparing to submit num_tasks=" ${num_tasks}

#-----------------SET UP THE C++ INPUT FILE ---------------------------
echo $Llim    >  void_numerics.in
echo $t0 $tf  >> void_numerics.in
echo $a       >> void_numerics.in

echo ""       >>  void_numerics.in

echo $should_plot_Vdist_v_t  $should_plot_rhos $should_check_neg >> void_numerics.in
echo $should_import_IC   $should_export_IC                       >> void_numerics.in

echo $should_peakterminate   $t_export >> void_numerics.in


echo ""                                >>  void_numerics.in

echo  $BZcond                          >>  void_numerics.in
echo  $parity_check                    >>  void_numerics.in
echo  $opath                           >>  void_numerics.in

#---------------------------------------------------------------------

# create a local directory for this job only
if [ ! -e $WORKDIR ]; then
    mkdir -p $WORKDIR
else
    # this should never happen
#   echo "Clean up /data/$USER/ directory on $HOSTNAME"
#        | mail -s "SGE-error" $USER
    exit 1
fi

# create a local copy of the program and start the job
OLDDIR=`pwd`

#---------COPY NECESSARY FILES OVER TO THE WORK DIRECTORY-------
cp void_numerics.x  $WORKDIR/
cp void_numerics.in $WORKDIR/

#-----------------  NOW BUILD THE SCRIPT FILE  -----------------

echo "#!/bin/bash"                     >   ${WORKDIR}${job_sub_script}
echo "#$ -S /bin/sh"                   >>  ${WORKDIR}${job_sub_script}

echo "#$ -cwd"                         >>  ${WORKDIR}${job_sub_script}
echo "#$ -m eas"                       >>  ${WORKDIR}${job_sub_script}
echo "#$ -l h_vmem=500M,s_rt=12:59:0"  >>  ${WORKDIR}${job_sub_script}
echo "#$ -t 1-"${num_tasks}            >>  ${WORKDIR}${job_sub_script}


echo ""                                 >>  ${WORKDIR}${job_sub_script}
echo "#-----------------------------"   >>  ${WORKDIR}${job_sub_script}
echo ""                                 >>  ${WORKDIR}${job_sub_script}
echo ""                                 >>  ${WORKDIR}${job_sub_script}

echo "NGtype=\""${NGtype}"\""           >>  ${WORKDIR}${job_sub_script}
echo "eps_array_file="${eps_array_file} >>  ${WORKDIR}${job_sub_script}
echo "muN_array_file="${muN_array_file} >>  ${WORKDIR}${job_sub_script}
echo "num_tasks="${num_tasks}           >>  ${WORKDIR}${job_sub_script}

echo  "num_eps="${num_eps}              >>  ${WORKDIR}${job_sub_script}
echo  "num_mus="${num_mus}              >>  ${WORKDIR}${job_sub_script}

cat carpark_jobsub_tail.sh              >>  ${WORKDIR}${job_sub_script}

#-----------GO TO THE WORK DIRECTORY AND EXECUTE ----------------
cd $WORKDIR

qsub ${job_sub_script}

#--------------------------------------------

