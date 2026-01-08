clear
size1=8
size2=14
size3=12
set(groot,'defaultAxesfontsize',size1)

x3=load('512/blist.dat') ;
x4=load('128k_tanning/blist.dat') ;

xx = 0: 0.05: 10;
xx = xx+1.e-6; %avoiding problem evaluating at 0
f = numone(xx);

figure(70)
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
t=256;

[values, edges]=histcounts(x3, 'BinMethod','integers','BinWidth', 1);
l=1:length(values);
p0= (edges(l)+edges(l+1))/2.;
p0= p0/t^(2/3);
integral= sum(values)*(edges(2)-edges(1))/t^(2/3);
values=values/integral;
bb=plot( p0 , values, 'LineWidth',1.5);

moment= 2*sum( values .* (p0(2)-p0(1)) );
hold on
f1=f*moment;
xx1  = xx/moment;
plot(xx1,f1,'r-.', 'LineWidth', 1.5)
plot(-xx1,f1,'r-.', 'LineWidth', 1.5)

xlabel('$x/t^{2/3}$', 'fontsize',size2,'interpreter','latex')
ylabel('$t^{2/3}\rho_1(x)$', 'fontsize',size2 ,'interpreter','latex' )
text(1,0.5,'t=256', 'interpreter','latex', 'fontsize', size3) 
axis([-2 2 0 0.6])
yticks([0.0 0.2 0.4 0.6])
pbaspect([1 .66 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

t=32*1024*8;
[values, edges]=histcounts(x4, 'BinMethod','integers', 'BinWidth', 60);
v0=values(1:10)
l=1:length(values);
p0= (edges(l)+edges(l+1))/2.;
e1=edges(1)
e2=edges(2)
de=e2-e1;
p0= p0/t^(2/3);
integral= sum(values)*(edges(2)-edges(1))/t^(2/3);
values=values/integral;
bb=plot( p0 , values, 'LineWidth',1.5);
moment= 2*sum( values .* (p0(2)-p0(1)) );
hold on
f2=f*moment;
xx  = xx/moment;
plot(xx,f2, 'r-.', LineWidth=1.5)
plot(-xx,f2, 'r-.', LineWidth=1.52)

axis([-2 2 0 0.8])
yticks([0.0 0.2 0.4 0.6])
xlabel('$x/t^{2/3}$', 'fontsize',size2,'interpreter','latex')
ylabel('$t^{2/3}\rho_1(x)$', 'fontsize',size2 ,'interpreter','latex' )
text(1,0.5,'t=262,144', 'interpreter','latex','fontsize', size3)
pbaspect([1 .66 1])

print('-dpdf', 'tasep.pdf')
