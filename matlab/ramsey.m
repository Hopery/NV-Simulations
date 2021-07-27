% Does a Ramsey experiment and displays it.
% Input: 
%   * `t` - Time for the time evolution. Usually 1/(40)/4, which means we
%   use a quater of the Omega, since a full period is 1/Omega we need
%   1/Omega/4, where Omega = 40
%   * `rho` - Density matrix `rho` (electron)
% Expected: Wave graphic between 1 and 0. If we compare it with the Rabi
% Rabi experiment must befaster than it but similar in the other ways.
% In the simulation we have a perfect wave, in reality it will decay
% as time keep going on so the wave will be smaller and smaller. 
% This represent the time we haveto perform gates until the decay rate 
% is too big (no visible in simulator since there is no noise)

function ramsey(t,rho)
    global H_x H_free
    
    % Fixed thime for time-evolution and graphic
    x = 0:1:1000;
    % Array with results
    nv_population = zeros(1,length(x));
    
    for kk = 1:length(x)
        % First, do a time evolution woth the Hamiltonian tha has the X
        % drive. This is because we want to rotate in the x axis
        % Apply the fixed time (piPulse - Laser for t time)
        rho_1 = time_evolution(H_x, t, rho); 
        % After it apply a time evolution to the free Hamiltonian 
        % (no drives) with the corresponding time of the loop.
        % This is because we dont want a graphic describing a circle, with
        % this we let the electron continue and we'll have a wave at the
        % end
        rho_2 = time_evolution(H_free,x(kk),rho_1); 
        % Finnaly apply time evolution to the Hamiltonian with drive x
        % and apply a time 3*t, since we want to go back to our current
        % position
        rho_3 = time_evolution(H_x,3*t,rho_2); 
        % To measure it do partial trace of the nucleus to keep the
        % electron
        rho_meas = PartialTrace(rho_3,2,[2,2]);
        % Like a proyection in zero, we only keep the first element of the
        % matrix
        nv_population(kk) = rho_meas(1,1);
        
    end
    
    % Plot
    figure('name','Ramsey')
    plot(x,nv_population, 'LineWidth',1);
    xlabel('Time')
    ylabel('spin state population')
    legend('m_s=0')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end


