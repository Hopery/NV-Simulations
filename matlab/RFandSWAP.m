% Applies a Swap then put electron's density matrix to [1 0;0 0], do an RF
% applies again a Swap ad measure.
%         SWAP - Ptr & kron -   RF      -  SWAP
% |0,?> - |?,up> - |0,up> - |0,up-down> - |0-(-1),up>
% Input: 
%   * `x0` - Components of the fourier transform.
%   * `t` - Time for the time evolution. Usually 1/(2*Omega)
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Returns: Time in the first lower point of the graphic
% Expected: Wave functions between 1 and 0, similar to the AXYYDispl
% function but slower.

function rfsw = RFandSWAP(x0,rho,times,t)
global H_13c

    time = 0:1:200;
    % Array with results
    population = zeros(1,length(time));

    % To get time of first lower point in the graphic.
    timeR = 0;
    check = true;
    
    for kk = 1:length(time)
        % Sequence
        rho_0 = SWAP(x0,rho,times,t);
        rho_13c = PartialTrace(rho_0,1,[2,2]);
        rho_1 = kron([1 0;0 0], rho_13c);
        % RF
        rho_1 = time_evolution(H_13c,time(kk), rho_1);
        rho_2 = SWAP(x0,rho_1,times,t);
        rho_meas = PartialTrace(rho_2,2,[2,2]);
        population(kk) = rho_meas(1,1);
        
       % To get time of first lower point in the graphic.
        if kk ~= 1 && check && population(kk) >  population(kk-1)
            timeR = time(kk);
            check = false;
        end 
    end
    
    rfsw = timeR;
    
    % Plots
    figure('name','RFandSWAP')
    plot(time,population, 'LineWidth',1);
    xlabel('Time')
    ylabel('spin state population')
    legend('m_s=0')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end

