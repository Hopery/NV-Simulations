% SWAP gate, do a swap between the nucleus and the electron
% Input: 
%   * `x0` - Components of the fourier transform.
%   * `t` - Time for the time evolution. Usually 1/(2*Omega), which means ...
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
% Returns: Measurement of the given density matrix with the SWAP applied.

function sw = SWAP(x0,rho,times,t)
global H_x H_y
    
    % Fourier times
    tau_f = AXYSearch(x0);
    
    % This should be closer to 1/ sqrt((gamma_13c*B-A_zz)^2 + 2*A_zx^2),
    % where thedenominator isequal to omega. The value was talen from the
    % graphic of SpinEchoFourier from itslower point. This is like this bc
    % ...
    tau = 1.885;

    % Get taus for AXYY function
    tau_1 = (tau_f(1)*tau)-t/2;
    tau_2 = (tau_f(2)- tau_f(1))*tau-t;
    tau_3 = (tau_f(3)- tau_f(2))*tau-t;
    
    % Sequence
    % pi/2 y
    rho_0 = time_evolution(H_y, 0.5*t, rho);
    % AXY pi/2 x
    for ii = 1:times
        rho_0 = XYX(t,rho_0,tau_f,tau);
    end
    % pi/2 y
    rho_1 = time_evolution(H_y, 0.5*t, rho_0 );
    % pi/2 x
    rho_2 = time_evolution(H_x, 0.5*t, rho_1 );
    % AXY pi/2 y
    for ii = 1:times
        rho_2 = AXYY(t,rho_2,tau_f,tau,[tau_1,tau_2,tau_3]);
    end
    % pi/2 x
    rho_3 = time_evolution(H_x, 0.5*t, rho_2 );
    
    sw = rho_3;
    rho_meas = PartialTrace(rho_3,2,[2,2]);
    % disp(rho_meas)
    rho_meas = rho_meas(1,1);
end

