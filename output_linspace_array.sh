dx=0.1
numx=50
x0=0.1

xi=0

outputfile="array_muN_a2.txt"

while [ $xi -le $numx ]
do

  x=$( echo "$x0 + $xi*$dx" |bc -l )

  echo $x >> $outputfile

  xi=$(echo "$xi +1" |bc -l)
done #-----------------------------------------------------------------------------------


