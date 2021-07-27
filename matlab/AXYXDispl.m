% Display the AXYX sequence
% Input: 
%   * `x0` - Components of the fourier transform.
%   * `t` - Time for the time evolution. Usually 1/(2*Omega)
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Returns: time where the first lower point in the graphic happens.
% Expected: Wave functions between 1 and 0, similar o ha rabi but much
% slower.

function SEF = AXYXDispl(x0,t,rho, times)
 global H_x % A_zz A_zx gamma_13c B 
 
    % Array with results
    nv_population = zeros(1,times);
    % Fourier times
    tau_f = AXYSearch(x0);
    
    % This should be closer to 1/ sqrt((gamma_13c*B-A_zz)^2 + 2*A_zx^2),
    % where the denominator is equal to omega. The value was taken from the
    % graphic of SpinEchoFourier from its lower point.
    tau = 1.885;
    
    % To get time of first lower point in the graphic.
    time = 0;
    check = true;
    
    % Similar structure as SpinEchoN and SpinEchoFourier
    for kk = 1:times
        rho_0 = time_evolution(H_x, 0.5*t, rho);
        % XYXYYXYX
        for ii = 1:kk
            rho_0 = XYX(t,rho_0,tau_f,tau);
        end
        % -------
        rho_10= time_evolution(H_x,3*0.5*t,rho_0);
        rho_meas = PartialTrace(rho_10,2,[2,2]);
        nv_population(kk) = rho_meas(1,1);
        
        % To get time of first lower point in the graphic.
        if kk ~= 1 && check && nv_population(kk) >  nv_population(kk-1)
            time = kk;
            check = false;
        end 
    end
    
    SEF = time;
    
    % Plots
    figure('name','AXYX')
    plot(1:times,nv_population, 'LineWidth',1);
    xlabel('Ntimes')
    ylabel('spin state population')
    legend('m_s=0')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end

