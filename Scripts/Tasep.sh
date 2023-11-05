#!/bin/bash
set -e -x

(cd ../CC; make lattice)

cd ..
mkdir Lattice.$$
cd Lattice.$$
../CC/lattice 2048 1024  # system size, k_iterations, note sweep is divided by 2 in code to avoid wrap-around in small systems
ln -s ../Matlab/nu1.m
ln -s ../Matlab/nu2.m
ln -s ../Matlab/airy0.m
ln -s ../Matlab/numone.m
ln -s ../Matlab/tasep.m
/Applications/MATLAB_R2023b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu1
/Applications/MATLAB_R2023b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu2
