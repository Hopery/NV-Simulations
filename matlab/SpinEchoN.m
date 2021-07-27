% Same as SpinEcho but  with a fixed tau and being repeited the experiments
% the desired times.
% Input: 
%   * `t` - Time for the time evolution. Usually 1/(2*Omega)
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Expected: Wave graphic going down in diagonal. Which means we are looking
% as with zoom a very big wave 

function SpinEchoN(t, times, rho)
    global H_y H_free A_zz A_zx  gamma_13c B
    
    % Tau time
    omega = sqrt((gamma_13c*B-A_zz)^2 + A_zx^2);
    tau = 1/(2*omega);
    nv_population = zeros(1,times);
    
    % Same as i SpinEcho but with fixed tau
    for kk = 1:times
        rho_1 = time_evolution(H_y, 0.5*t, rho); % y/2
        
        rho_2 = time_evolution(H_free,tau/2,rho_1); % tau
        
        rho_3 = SpinEchoAux(kk,t,tau,rho_2);  % xyxyyxyx

        rho_4 = time_evolution(H_y,3*0.5*t,rho_3); % y/2
        rho_meas = PartialTrace(rho_4,2,[2,2]);
        nv_population(kk) = rho_meas(1,1);
    end
    
    % Plot
    figure('name','SpinEchoN')
    plot(1:times,nv_population, 'LineWidth',1);
    xlabel('tau [\mus]')
    ylabel('spin state population')
    legend('m_s=0')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end
