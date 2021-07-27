% Display the AXYY sequence
% Input: 
%   * `x0` - Components of the fourier transform.
%   * `t` - Time for the time evolution. Usually 1/(2*Omega)
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Expected: Wave functions between 1 and 0, similar to the rabi but much
% slower and equal to the AXYXDispl graphic

function AXYYDispl(x0,t,rho,times)
    global H_x
    
    % This should be closer to 1/ sqrt((gamma_13c*B-A_zz)^2 + 2*A_zx^2),
    % where the denominator is equal to omega. The value was taken from the
    % graphic of SpinEchoFourier from itslower point. 
    tau = 1.885;
    
    % Fourier times
    tau_f = AXYSearch(x0);
    
    % Get taus for AXYY function
    tau_1 = (tau_f(1)*tau)-t/2;
    tau_2 = (tau_f(2)- tau_f(1))*tau-t;
    tau_3 = (tau_f(3)- tau_f(2))*tau-t;
    
    % Array with results
    nv_population = zeros(1,times);
    
    % Same structure as AXYXDispl
    for kk = 1:times
        rho_00 = time_evolution(H_x, 0.5*t, rho);
        for ii = 1:kk
            rho_00 = AXYY(t,rho_00,tau_f,tau,[tau_1 tau_2 tau_3]);
        end
        rho_2 = time_evolution(H_x,3*0.5*t,rho_00);
        rho_meas = PartialTrace(rho_2,2,[2,2]);
        nv_population(kk) =  rho_meas(1,1);
    end
    
    % Plots
    figure('name','AXYY')
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

