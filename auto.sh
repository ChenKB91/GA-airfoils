#!/bin/bash
# /Applications/MATLAB_R2015b.app/bin/matlab -nodesktop -nosplash -r init

# --- TIME TABLE ---
# 1 hr = 3600 seconds
# 2 hr = 7200 seconds
# 3 hr = 10800 seconds
# 4 hr = 14400 seconds
# 5 hr = 18000 seconds
# 6 hr = 21600 seconds
# 7 hr = 25200 seconds
# 8 hr = 28800 seconds
# 9 hr = 32400 seconds
# read -p "Input Seconds to run... -> " tmp
# end=$((SECONDS+$tmp))
read -p "Input Generations to run... ->" gen
# while [ $SECONDS -lt $end ]; do
for i in $(seq 1 $gen);
do
    # Do what you want.
    # python move_file.py 0
    cp -f data/ib*.inp input/
    make cleandata
    bin/ib
    cp -f output/responds/gen_force.dat data/force.txt
    /Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nojvm -nosplash -r "gen_foil; exit;"
    
    # echo "gen_foil($i)"
done
echo "finishing time $SECONDS seconds"

