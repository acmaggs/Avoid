#!/bin/bash
set -e -x

(cd ../CC; make avoid)

cd ..
mkdir Avoid.$$
cd Avoid.$$
../CC/avoid 512 1024   # system size, k_iterations
ln -s ../Matlab/nu1.m
ln -s ../Matlab/nu2.m
ln -s ../Matlab/airy0.m
ln -s ../Matlab/numone.m
/Applications/MATLAB_R2023b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu1
/Applications/MATLAB_R2023b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu2
