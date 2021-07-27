% Does an SpinEcho experiment and displays it. Spin echo is used to
% decouple the electro sping from the nuecleus spin and the enviroment 
% (XY8 secuence), also with the dephasing it increase the decay time so we
% have more time to apply gates.
% Input: 
%   * `t` - Time for the time evolution. Usually 1/(2*Omega)
%   * `times` - Times the sequence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Expected: Wave graphic similar to an odmr but with several peaks. That 
% represents we have to much interactions and noise.

function SpinEcho(t,times, rho)
    global H_y H_free
    
    % Fixed thime for time-evolution and graphic
    tau = 0.7:1.e-3:1.2;
    % Array with results
    nv_population = zeros(1,length(tau));
    
    for kk = 1:length(tau)
        % First, do a time evolution woth the Hamiltonian tha has the y
        % drive. 
        % Apply the fixed time 0.5*t(piPulse - Laser for t time). This
        % to get  1/Omega/4
        rho_1 = time_evolution(H_y, 0.5*t, rho); % y/2
        % After it apply a time evolution to the free Hamiltonian 
        % (no drives) with the corresponding time of the loop divided by 2.
        rho_2 = time_evolution(H_free,tau(kk)/2,rho_1); % tau
        % Apply XYXYYXYX the desired  number of times with the 
        % corresponding time of the loop. (See  `SpinEchoAux` function for 
        % sequence explanation)
        rho_3 = SpinEchoAux(times,t,tau(kk),rho_2);  % xyxyyxyx
        % Finnaly apply time evolution to the Hamiltonian with drive y
        % and apply a time 3*0.5*t, since we want to recover our initial
        % position
        rho_4 = time_evolution(H_y,3*0.5*t,rho_3); % y/2
        % To measure it do partial trace of the nucleus to keep the
        % electron
        rho_meas = PartialTrace(rho_4,2,[2,2]);
        % Like a proyection in zero, we only keep the first element of the
        % matrix
        nv_population(kk) = rho_meas(1,1);
        
    end
    
    % Plot
    figure('name','SpinEcho')
    plot(tau,nv_population, 'LineWidth',1);
    xlabel('tau [\mus]')
    ylabel('spin state population')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end



