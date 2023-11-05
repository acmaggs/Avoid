clear 
load 32/blist.dat
x(1)=32
y(1)=skewness(blist)
z(1)=var(blist)

load 64/blist.dat
x(2)=64
y(2)=skewness(blist)
z(2)=var(blist)

load 128/blist.dat
x(3)=128
y(3)=skewness(blist)
z(3)=var(blist)

load 256/blist.dat
x(4)=256
y(4)=skewness(blist)
z(4)=var(blist)

load 512/blist.dat
x(5)=512
y(5)=skewness(blist)
z(5)=var(blist)

load 1024/blist.dat
x(6)=1024
y(6)=skewness(blist)
z(6)=var(blist)

load 2048/blist.dat
x(7)=2048
y(7)=skewness(blist)
z(7)=var(blist)

load 4k/blist.dat
x(8)=4*1024
y(8)=skewness(blist)
z(8)=var(blist)

load 8k/blist.dat
x(9)=8*1024
y(9)=skewness(blist)
z(9)=var(blist)

load 16k/blist.dat
x(10)=16*1024
y(10)=skewness(blist)
z(10)=var(blist)

load 32k/blist.dat
x(11)=32*1024
y(11)=skewness(blist)
z(11)=var(blist)

load 64k/blist.dat
x(12)=64*1024
y(12)=skewness(blist)
z(12)=var(blist)

if(0)
load 65k/blist.dat
x(13)=64*1024
y(13)=skewness(blist)
z(13)=var(blist)
end

load 128k/blist.dat
x(13)=1024*128
y(13)=skewness(blist)
z(13)=var(blist)

load 256k/blist.dat
x(14)=256*1024
y(14)=skewness(blist)
z(14)=var(blist)

load 512k/blist.dat
x(15)=512*1024
y(15)=skewness(blist)
z(15)=var(blist)
%%%%%%%%%%%%%%
size=16
set(groot,'defaultAxesfontsize',12)
y2=abs(y); % skew
x2=x;
z2=sqrt(z); %  var

if(0)
%x2(13)=x(12)
y2(13)=y(12)
z2(13)=sqrt(z(12))
end


figure(70)
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');

nexttile
loglog(x2,z2,'+', x2, x2.^(2/3)/3 ,'LineWidth', 2)
xlabel('t', 'FontSize', size)
ylabel('standard deviation', 'FontSize', size)

nexttile
loglog(x2,y2,'+',x2, 3*x2.^(-.333),'LineWidth', 2)
xlabel('t', 'FontSize', size)
ylabel('skewness', 'FontSize', size)
print('-dpdf', 'moments.pdf')
