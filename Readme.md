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
