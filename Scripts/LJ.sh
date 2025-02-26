#!/bin/bash
set -e -x

(cd ../Paper2/lj; make)

cd ../
mkdir LJ.$$
cd LJ.$$
ln -s ../Paper2/lj/bond_dist.m
ln -s ../Paper2/lj/w1.m
#create the set of bond-length from boltzmann distribution
/Applications/MATLAB_R2024b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch bond_dist

../Paper2/lj/simul
ln -s ../Matlab/nu1.m
ln -s ../Matlab/nu2.m
ln -s ../Matlab/airy0.m
ln -s ../Matlab/numone.m
/Applications/MATLAB_R2024b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu1
/Applications/MATLAB_R2024b.app/bin/matlab  -nodisplay -nodesktop -nosplash -batch nu2
