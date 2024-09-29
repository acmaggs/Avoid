clear
expon=1./4.;
load resume
size=resume(:,1);
ratio=resume(:,3) ./ resume(:,2);
ratio = ratio .^ 1.5;

ratio0=ratio;
ratio0=ratio0(4:end);
size0=size(4:end);
modelfun = @(b,x)( b(1) + b(2)*(x.^(-b(3))) );

opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';

beta0 = [.5;-0.4;.25];
x=size0;
y=ratio0;
beta = nlinfit(x,y,modelfun,beta0,opts)                                                                                          
b=beta;
z=( b(1) + b(2)*(x.^(-b(3))) ); 
figure(123)
l=plot (1. ./ size.^expon ,ratio, 'o', 1. ./x.^expon, z, '-.','linewidth',2);
l(1).MarkerSize=8;
hold on
xlabel('$1/t^{1/4}$','interpreter','latex','fontsize',16)
ylabel('$(\sigma_{\textrm{hot}}/\sigma_{\textrm{cold}})^3$','interpreter','latex','fontsize',16)
%axis([0 0.17 0.4 0.62])
axis([0 0.27 0.4 0.62])
text(.05, .55,'lifted TASEP','fontsize', 16)
text(.05, .45,'Gaussian','fontsize', 16)


load gaussian
size=gaussian(:,1);
ratio=gaussian(:,3) ./ gaussian(:,2);
ratio = ratio .^ 1.5;


ratio0=ratio;
ratio0=ratio0(4:end);
size0=size(4:end);
modelfun = @(b,x)( b(1) + b(2)*(x.^(-b(3))) );
set(gca,'ColorOrderIndex',1)
beta0 = [.5;+0.4;.25];
x=size0;
y=ratio0;
beta = nlinfit(x,y,modelfun,beta0,opts)
b=beta;
z=( b(1) + b(2)*(x.^(-b(3))) );
l=plot (1. ./ size.^expon ,ratio, ' o', 1. ./x.^expon, z, '-.','linewidth',2);
l(1).MarkerSize=8;
hold off
exportgraphics(gca,'sigmaratio.pdf','BackgroundColor','none')
