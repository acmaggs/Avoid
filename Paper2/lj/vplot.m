clear
figure(1);
load virial.dat; 
virial=virial / 1.58376;
histogram(virial);
mean(virial)

figure(2);
m=bootstrp(300,@mean,virial);
[fi,xi] = ksdensity(m);
plot(xi,fi)
