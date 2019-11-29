#!/bin/bash
#
# See if required file is present
#
check1=`ls | grep band.yaml | wc -l`
if [ $check1 -eq 0 ]
then
   echo "Cannot find phonopy file band.yaml - bye!!"
   exit 1
fi
echo "Found band.yaml. Extracting distances along path in BZ..."
#
# Get the number of q points
#
nqpt=`head -1 band.yaml | awk '{print $2}'`
#
# Get the number of frequencies
#
str1=`grep -n "q-position" band.yaml | head -2 | tail -1 | awk '{print $1}'`
len1=${#str1}
len1=`expr $len1 - 2`
nline1=${str1:0:$len1}
nfreq=`head -$nline1 band.yaml | grep "frequency" | wc -l | awk '{print $1}'`
#
# Write header line to output
#
echo "# distance along path vs fequencies, number of points:" $nqpt  > disp.dat
#
# Extract the distances along the path
#
grep distance band.yaml | awk '{print $2}' > temp_dist
#
# Set an initial old distance to be very large
#
distold=10000000
echo "...now extracting frequencies..."
#
# Loop over q-points
#
for (( i=1; i<=$nqpt; i++ ))
do
#
# Extract appropriate distance and check that it is not redundant (i.e. alrealdy 
# printed on the previous line)
#
   dist=`tail -n +$i temp_dist | head -1 | awk '{print $1}'`
   diffdist=`echo "$dist - $distold" | bc -l`
   len1=${#diffdist}
   if [ $len1 -eq 1 ] && [ $diffdist -eq 0 ]
   then
      continue
   fi 
   echo $dist > temp1
#
# Now extract frequencies for each q-point. Initially write them to a temp file and 
# convert so that they appear on the same line. Then append this line to output
#
   for (( j=1; j<=$nfreq; j++ ))
   do
      num1=$(($j + ($i - 1) * $nfreq))
      freqtmp=`grep "frequency" band.yaml | tail -n +$num1 | head -1 | awk '{print $2}'`
      echo "    " $freqtmp >> temp1
   done
   num1=`wc -l temp1 | awk '{print $1}'`
   while [ $num1 -gt 1 ]
   do
      sed '$!N;s/\n//' temp1 > temp2
      mv temp2 temp1
      num1=`wc -l temp1 | awk '{print $1}'`
   done
   cat temp1 >> disp.dat
   rm temp1
   distold=$dist
done
echo "...done. File format is distance vs frequencies for each q point."
rm temp*   
#xmgrace -nxy disp.dat &
