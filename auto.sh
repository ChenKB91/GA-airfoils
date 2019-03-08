#!/bin/bash
# /Applications/MATLAB_R2015b.app/bin/matlab -nodesktop -nosplash -r init

end=$((SECONDS+21600))
a=0
while [ $SECONDS -lt $end ]; do
    # Do what you want.
    # python move_file.py 0
    cp -f data/ib*.inp input/
    bin/ib
    cp -f output/responds/gen_force.dat data/force.txt
    /Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nosplash -r gen_foil
    
done