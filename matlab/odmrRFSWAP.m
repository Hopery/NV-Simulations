% ODMR of function RFandSWAP.
% Input: 
%   * `x0` - Components of the fourier transform.
%   * `t` - Time for the time evolution. Usually 1/(2*Omega)
%   * `timesSW` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Expected: Wave functions with one big peack that reachs zero

function odmrRFSWAP(x0,rho,timesSW,t)
global H_13c_w B gamma_13c
    
    % Time where the first lower point of the RFandSWAP wave happens
    tau = 51;
    % Nucleus omega that will be use in the Hamiltonian
    omega = gamma_13c * B-0.1:5e-3:gamma_13c * B+0.1;
    % Array with results
    population = zeros(1,length(omega));
    
    % Same as in RFandSWAP function
    for kk = 1:length(omega)
        % Sequence
        rho_0 = SWAP(x0,rho,timesSW,t);
        rho_13c = PartialTrace(rho_0,1,[2,2]);
        rho_1 = kron([1 0;0 0], rho_13c);
        % RF
        rho_1 = time_evolution(H_13c_w(omega(kk)),tau, rho_1);
        rho_2 = SWAP(x0,rho_1,timesSW,t);
        rho_meas = PartialTrace(rho_2,2,[2,2]);
        population(kk) = rho_meas(1,1);
    end
    
    % Plot
    figure('name','RFandSWAPODMR')
    plot(omega,population, 'LineWidth',1);
    xlabel('Frequency')
    ylabel('spin state population')
    legend('m_s=0')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end

