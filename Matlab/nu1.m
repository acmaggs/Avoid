function y=nu1()  % symmetric displacement curve
    set(groot,'defaultAxesfontsize',16)
    load blist.dat
    width = sqrt(var(blist))     
    binwidth = 2; %must be even
    binwidth2 = 0.9;
    lw=3;
    [values, edges] = histcounts(abs(blist), 'BinWidth',binwidth);

    v_integral = sum (values) - values(1)/2;
    x = 0:length(values)-1;
    x = x+.5;
    x = x/v_integral/2; % rescale x-axis to unit normalization in [-infty, infty]
    v_moment = sum( x .* values) *( x(2)-x(1) );
   
    xx = 0: 0.05: 10;
    xx = xx+1.e-6; %avoiding problem evaluating at 0
    f = numone(xx); % theoretical curve
    integral_nu = (sum(f)-f(1)/2) * (xx(2)-xx(1)) % should be 1/2
    nu_moment = sum (xx*f)*(xx(2)-xx(1));
    
    v_moment = v_moment/nu_moment; %scale data to fix first moment
    values = values * v_moment;
    x = x / v_moment;
    
    figure(6)
    bar(x, values, 'EdgeColor' , 'none', 'BarWidth',binwidth2) % rescaled histogram
    hold on
    plot(xx,f, 'linewidth',lw)
    hold off
    axis([0 5 0 .25])

    xlabel('$|x|$') %, 'fontsize',18,'interpreter','latex')
    ylabel('$\nu_1(x)$', 'fontsize',24 ,'interpreter','latex' )
    print('-dpdf', 'f4.pdf')

    
    %%%%%%%%%%%%%%%%%
    biggest=max(blist);
    smallest=min(blist);
    if (abs(biggest) > abs(smallest))
        blist(end+1) = -biggest;
    else
        blist(end+1) = -smallest;
    end
    [values, edges]=histcounts(blist, 'BinMethod','integers', 'BinWidth',binwidth);

    x = -length(values)/2: length(values)/2-1;
    values = 2*v_moment*values;
    x = x/v_integral/2/v_moment;
    
    figure(7)    
    bar(x, values, 'EdgeColor' , 'none', 'BarWidth',binwidth2);
    hold on
    plot(xx,f, 'r', 'linewidth', lw);
    plot(-xx,f, 'r', 'linewidth', lw);
    hold off

    xlabel('$x$', 'fontsize',24 ,'interpreter','latex')
    ylabel('$\nu_1(x)$',  'fontsize',24 ,'interpreter','latex')
    axis([-4 4 0 .3])
    exportgraphics(gca,'nu1.pdf','BackgroundColor','none')
end
