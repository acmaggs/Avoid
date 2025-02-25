function y = ww(x,P,beta)  % distribution function for LJ bonds
x2 = 1 ./(x.^2);
x6 = x2 .^ 3;
x12 = x6.^2;
V= x12 -x6;
y= exp( -beta*(V +P*x));



