# Non-reversible Monte Carlo: an example of “true” self-repelling motion

Available from:
[https://hal.science/hal-04265583](https://hal.science/hal-04265583) or
[https://arxiv.org/abs/2310.19494](https://arxiv.org/abs/2310.19494)

## Abstract
We link the large-scale dynamics of non-reversible Monte Carlo algorithms as well as a lifted TASEP to an exactly soluble model of self-repelling motion. We present arguments for the connection between the problems and perform simulations, where we show that the empirical distribution functions generated from Monte Carlo are well described by the analytic solution of self-repelling motion.

### Code
`CC/avoid.cc`  code for harmonic chain.

`CC/lattice.cc`  code for TASEP chain.

The TASEP code can be compiled with, or without graphical display, by changing the variable `Param::GRAPH`.

`Matlab/nu1.m` for the analysis of data, comparing to $\rho_1$ of eq.(4)

`Matlab/nu2.m` for the analysis of data, comparing to $\rho_2$ of eq.(5)

`Scripts/Harmonic.sh` and `Scripts/Tasep.sh` 
run the two c++ programs, and compares
the results with the analytic curves computed using Matlab.

The scripts are for MacOS but can trivially be changed to run on Linux.

Airy functions are from  Meg Noah (2023). Negative Zeros of the Airy Functions (https://www.mathworks.com/matlabcentral/fileexchange/88848-negative-zeros-of-the-airy-functions).


# Event-chain Monte Carlo and the true self-avoiding walk
## Abstract
We study the large-scale dynamics of event chain Monte Carlo algorithms in one dimension, and their relation to the true self-avoiding walk. In particular, we study the influence of stress, and different forms of interaction on the equilibration and sampling properties of algorithms with global balance, but no local balance. 

[https://arxiv.org/abs/2410.08694](https://arxiv.org/abs/2410.08694)
Code for second paper in directory `Paper2/lj` and  `Paper2/lambert`


# Test runs

The `Scripts` directory contains four simple bash scripts that run harmonic, tasep, lambert, lennard-jones systems for just a few minutes and plot the distribution functions using Matlab.
