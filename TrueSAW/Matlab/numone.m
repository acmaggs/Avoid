%Refs to https://doi.org/10.1016/j.spa.2012.11.011
function s=numone(x, roots, pk)
% eq(15), eq(43), eq(47) Using N roots of the Airy function
     x=abs(x);
     N=1:20;
     roots = -airy0(1,N)/2^(1/3); % 2^{1/3} since the equation in the paper
                                  % eq(41) has a factor 2, compared to standard
                                  % definitions of Airy.
     up=6^(1/3)* gamma(2/3)/gamma(1/3);
     pk= up^2/2 ./ roots.^4;
     s= sum( pk .* roots .* mittag( x' * roots) , 2)/2.;
end

function m=mittag(x)
% eq (53)
    m = (2^(1/3)/sqrt(3*pi)) .* x .* exp(-4*x.^3/27) .* kummerU(1/6,4/3, 4 *x.^3/27);
end
