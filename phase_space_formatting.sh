nj=10
#-----the number of independent jobs we consider

if [ $# -ne 1 ]
  then
    echo "Need one input argument to specify NGtype"
    exit
fi

#---------------------------------------------------------

nt=130
Llim=3000

Fixedref=302
endline=$( echo "$Llim+$Fixedref"|bc -l)	#--- last line from the source that makes 
						#--- it to the final comparison.

# the fixed particle is at line 302, so we start counting from 303
# the system is 3000 long so we count up to line 3303


#-----the number of columns (time points) that we consider
ind_filename="Fixed-occ-hist-ti_x_N.txt"
NGtype=$1
output_filename="delta2_vs_t"

echo "executing bash script for NGtype= " $NGtype
#-------------------------------------------------


#---this is the final output file.

TM_eqreffile_irho155="TM_"${NGtype}"_eq_irho155.txt"
TM_eqreffile_irho165="TM_"${NGtype}"_eq_irho165.txt"
TM_eqreffile_irho180="TM_"${NGtype}"_eq_irho180.txt"

MC_tifile_irho155="Fixed_average_"${NGtype}"_it_irho155.txt"
MC_tifile_irho165="Fixed_average_"${NGtype}"_it_irho165.txt"
MC_tifile_irho180="Fixed_average_"${NGtype}"_it_irho180.txt"

awk '{print $2}' ${TM_eqreffile_irho155} > temp_eq_155
awk '{print $2}' ${TM_eqreffile_irho165} > temp_eq_165
awk '{print $2}' ${TM_eqreffile_irho180} > temp_eq_180


#----------------------------------------------


if [ -d "temp_deltati_155" ]
then
    rm    temp_deltati_155  
fi

if [ -d "temp_deltati_165" ]
then
    rm    temp_deltati_165  
fi

if [ -d "temp_deltati_180" ]
then
    rm    temp_deltati_180  
fi

#----------------------------------------------

rho_eq155=$(echo "1.0/155.0" | bc -l)
rho_eq165=$(echo "1.0/165.0" | bc -l)
rho_eq180=$(echo "1.0/180.0" | bc -l)

#---loop over number of time points.

t=1
while [ $t -le $nt ]
do

  
  awk '{print $'$t'}' ${MC_tifile_irho155} | head -${endline} |tail -${Llim} > temp_ti_155
  awk '{print $'$t'}' ${MC_tifile_irho165} | head -${endline} |tail -${Llim} > temp_ti_165
  awk '{print $'$t'}' ${MC_tifile_irho180} | head -${endline} |tail -${Llim} > temp_ti_180

  #---loop over number of jobs being considered

  pr -m -t -s\  temp_ti_155 temp_eq_155 > temp_dif_155
  pr -m -t -s\  temp_ti_165 temp_eq_165 > temp_dif_165
  pr -m -t -s\  temp_ti_180 temp_eq_180 > temp_dif_180


  #--- consolidates all of the different jobs' data 
  #--- at this time point into a single file.

  awk '{s+=(sqrt(($2-$1)^2))/(('$Llim'+1)*'$rho_eq155')}END{print s}'  temp_dif_155 >> temp_deltati_155
  awk '{s+=(sqrt(($2-$1)^2))/(('$Llim'+1)*'$rho_eq165')}END{print s}'  temp_dif_165 >> temp_deltati_165
  awk '{s+=(sqrt(($2-$1)^2))/(('$Llim'+1)*'$rho_eq180')}END{print s}'  temp_dif_180 >> temp_deltati_180

  echo "finished time point " $t " of " $nt

  t=$(echo "$t +1" |bc -l)

done #-----------------------------------------------------------------------------------

awk '{print $1}'  "filling_averaged_irho155_"${NGtype}".txt" > timepoints

#--------------------FINISHED GETTING deltas into file------------------


pr -m -t -s\  timepoints temp_deltati_155 > patternconverge-irho155_t_d2.txt
pr -m -t -s\  timepoints temp_deltati_165 > patternconverge-irho165_t_d2.txt
pr -m -t -s\  timepoints temp_deltati_180 > patternconverge-irho180_t_d2.txt

# rm timepoints
rm temp_*

