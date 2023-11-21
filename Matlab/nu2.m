function nu2()
    binwidth = 1
    binwidth2 = 0.9;
    lw=3;
    set(groot,'defaultAxesfontsize',16)


    load list.dat
    width = sqrt(var(list))
    list = list-1;
    [values, edges] = histcounts(list, 'BinWidth',binwidth);
    mn_entry = mean(list)
  
    step=0.01;
    h = 0: step : 2;
    u = exp(-8*h.^3/9).*kummerU(1/6,2/3, 8* h.^3/9);
    u = u * (2 * 6^(1/3) * sqrt(pi)) /( (gamma(1/3) )^2 );
    integral_U = (sum(u)-u(1)/2)*step %check this should be unity
    exact_moment = sum(h .* u) * step;
    
    vnorm = sum(values) - values(1)/2. ; %  Trapezoidal rule
    x = 0:length(values)-1;
    x = x/vnorm; % scale x-axis so integral is now unity

    %rescale curve to impose "exact_moment"
    v_moment = sum( x .* values ) * (x(2)- x(1))/exact_moment;
    values = values * v_moment;
    x = x / v_moment;
    
    figure(3)
    bar(x, values, 'EdgeColor', 'none', 'BarWidth',binwidth2) 
    hold on
    plot(h, u, 'linewidth', lw) % analytic result
    hold off

    xlabel('$h$', 'fontsize', 24,'interpreter','latex')
    ylabel('$\nu_2(h)$', 'fontsize', 24,'interpreter', 'latex')
    axis([0 2 0 1.4])
    exportgraphics(gca,'nu2.pdf','BackgroundColor','none')
end
