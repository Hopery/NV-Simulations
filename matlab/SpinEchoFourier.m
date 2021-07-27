% A SpinEcho experiment but  we got the times by fourier componentsans the
% secuency XYXYYXYX is componsed by sequency each pulse (for more info see
% `pi_x` and `pi_y` functions). The difference between XYXYYXYX simple
% sequences and composed XYXYYXYX sequences is that the noise is reduced in
% the second one.
% Input: 
%   * `x0` - Components of the fourier transform.
%   * `t` - Time for the time evolution. Usually 1/(2*Omega)
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Expected: Graphic similar to an ODMR. This graphic shows us the same 
% graphic as in SpinEcho but with less noise (interactions)

function SpinEchoFourier(x0,t,rho,times)
    global H_x % A_zz A_zx gamma_13c B 
    
    % Fixed thime for time-evolution and graphic, times are bigger than
    % before bc it need to be twice the time
    tau = 1.6:0.005:2.1;
    % Array with results
    nv_population = zeros(1,length(tau));
    % Fourier times
    tau_f = AXYSearch(x0);
    
    % Same as in SpinEcho but using Fourier timesandfor the XYXYYXYX we use
    % other sequences for x and y
    parfor kk = 1:length(tau)
        rho_0 = time_evolution(H_x, 0.5*t, rho);
        % XYXYYXYX
        for ii = 1:times
            rho_2 = pi_x(t,rho_0,tau_f,tau(kk));
            rho_3 = pi_y(t,rho_2,tau_f,tau(kk));
            rho_4 = pi_x(t,rho_3,tau_f,tau(kk));
            rho_5 = pi_y(t,rho_4,tau_f,tau(kk));
            rho_6 = pi_y(t,rho_5,tau_f,tau(kk));
            rho_7 = pi_x(t,rho_6,tau_f,tau(kk));
            rho_8 = pi_y(t,rho_7,tau_f,tau(kk));
            rho_9 = pi_x(t,rho_8,tau_f,tau(kk));
            rho_0 = rho_9;
        end
        % -------
        rho_10= time_evolution(H_x,3*0.5*t,rho_9);
        rho_meas = PartialTrace(rho_10,2,[2,2]);
        nv_population(kk) = rho_meas(1,1);
    end
    
    % Plots
    figure('name','SpinEchoFourier')
    plot(tau,nv_population, 'LineWidth',1);
    xlabel('tau [\mus]')
    ylabel('spin state population')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end

