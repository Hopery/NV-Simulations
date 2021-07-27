% Does a Rabi experiment and displays it.
% Input: Density matrix `rho` (electron)
% Expected: Wave graphic between 1 and 0. Representing how the spin is
% changing between those two values.

function rabi(rho)
    global H
    
    % Fixed time for time-evolution and graphic
    time = 0.001:0.0001:0.05;
    % Array with results
    measurements = zeros(1,length(time));
    
    for i = 1:length(time)
        % Do a time evolution with the full hamiltonian and measure it
        meas = time_evolution(H,time(i),rho);
        % To measure it do partial trace of the nucleus to keep the
        % electron
        meas = PartialTrace(meas,2,[2,2]);
        % Like a proyection in zero, we only keep the first element of the
        % matrix
        measurements(i) = meas(1,1);
    end
       
    % Plot
    figure('name','Rabi')
    plot(time, measurements, 'LineWidth',1);
    xlabel('Time')
    ylabel('Measurements')
    box on;
    ax=gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.ZAxis.FontSize = 15;
end

