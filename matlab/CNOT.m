% Does a CNOT, being the nucleous the tarjet and the electron the controlled
% Input: 
%   * `x0` - Components of the fourier transform.
%   * `t` - Time for the time evolution. Usually 1/(2*Omega), which means ...
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Returns: Measurement of the given density matrix with the CNOT applied.

function cnt = CNOT(x0,rho,times,t)
    global H_x H_y
    
    % Fourier times
    tau_f = AXYSearch(x0);
    
    % This should be closer to 1/ sqrt((gamma_13c*B-A_zz)^2 + 2*A_zx^2),
    % where the denominator is equal to omega. The value was taken from the
    % graphic of SpinEchoFourier from itslower point.
    tau = 1.885;

    % Do a first time evolutionwith the Hamiltonian with the X drive and
    % multiplying the time by  0.5.
    rho_0 = time_evolution(H_x, 0.5*t, rho);
	% Apply XYXYYXYX A sequence
    for ii = 1:times
    	rho_2 = pi_x(t,rho_0,tau_f,tau);
    	rho_3 = pi_y(t,rho_2,tau_f,tau);
    	rho_4 = pi_x(t,rho_3,tau_f,tau);
     	rho_5 = pi_y(t,rho_4,tau_f,tau);
    	rho_6 = pi_y(t,rho_5,tau_f,tau);
    	rho_7 = pi_x(t,rho_6,tau_f,tau);
    	rho_8 = pi_y(t,rho_7,tau_f,tau);
        rho_9 = pi_x(t,rho_8,tau_f,tau);
        rho_0 = rho_9;
    end
    % -------
    % Last time evolution with the Hamiltonian with the Y drive in order to 
    % do a reverse, if we apply an X drive we won't have our desired
    % result. Multiplying the time by 0.5 in this way we "flip" the spin in
    % y axis.
	rho_10= time_evolution(H_y,0.5*t,rho_9);
    % To measure it do partial trace of the nucleus to keep the
    % electron
	rho_meas = PartialTrace(rho_10,2,[2,2]);
    % Like a proyection in zero, we only keep the first element of the
    % matrix
	cnt = rho_meas(1,1);
	disp(rho_meas)
end

