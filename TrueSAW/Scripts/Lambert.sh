#!/bin/bash
set -e -x

(cd ../Paper2/lambert; make)

cd ../
mkdir Lambert.$$
cd Lambert.$$
../Paper2/lambert/lambert 
ln -s ../Matlab/nu1.m
ln -s ../Matlab/nu2.m
ln -s ../Matlab/airy0.m
ln -s ../Matlab/numone.m
/Applications/MATLAB_R2024b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu1
/Applications/MATLAB_R2024b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu2
