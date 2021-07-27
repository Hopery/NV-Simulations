function RF(rho)
    global H
    
    time = 0:1:200;
    
    population_13c = zeros(1,length(time));
    
    for kk = 1:length(time)
        rho_0 = time_evolution(H,time(kk), rho);
        rho_meas = PartialTrace(rho_0,1,[2,2]);
        population_13c(kk) = rho_meas(1,1);
    end
    
    % Plots
    figure('name','RF');
    plot(time,population_13c, 'LineWidth',1);
    xlabel('Time')
    ylabel('spin state population')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end

