clear
clc
N=1.e8;
system('rm delta.h5')

P=0.03;
beta=8;
limit=20/P;

addpath('/Users/tony2/github/chebfun'), savepath 


dist=chebfun( @(x) w1(x,P,beta), [.1,limit]);
normalization=sum(dist);
m1=chebfun( @(x) x.*w1(x,P,beta)/normalization, [.1,limit]);
mean_bond_length=sum ( m1 )

dist=dist/normalization;

cumul=cumsum(dist);
cdfinv = inv(cumul, 'splitting','on');
title('distribution')

figure(10)
plot(cdfinv);
title('inverse cumulative distribution');

tic
r= squeeze(rand ([1 N ]));
sample=single (cdfinv( r ));

toc
mean_samples =mean( sample )
delta=mean_samples-mean_bond_length % should be small

h5create("delta.h5","/delta",[  N ], 'Datatype','single')
h5create("delta.h5","/P",[ 1 ])
h5create("delta.h5","/beta",[ 1 ])
h5write('delta.h5','/delta',sample)
h5write('delta.h5','/P',P)
h5write('delta.h5','/beta',beta)
