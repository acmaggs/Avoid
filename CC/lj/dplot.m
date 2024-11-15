clear
figure(1);
load displ.dat; 
histogram(displ);
mean(displ)

figure(2);
m=bootstrp(300,@mean,displ);
[fi,xi] = ksdensity(m);
plot(xi,fi)
