% XYXYYXY where X and Y are composed by same sequence
% Input: 
%   * `t` - Time for the time evolution. Usually 1/(2*Omega), which means ...
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho_0` - Density matrix `rho` (electron)
%   * `tau_f` - Tau times given by the fourier components
%   * `tau` - Time given for the loop of last function
% Returns: density matrix after sequence.

function xyy = XYY(t,rho_0,tau_f,tau)
    rho_2 = pi_x(t,rho_0,tau_f,tau);
    rho_3 = pi_y(t,rho_2,tau_f,tau);
    rho_4 = pi_x(t,rho_3,tau_f,tau);
    rho_5 = pi_y(t,rho_4,tau_f,tau);
    rho_6 = pi_y(t,rho_5,tau_f,tau);
    rho_7 = pi_x(t,rho_6,tau_f,tau);
    rho_8 = pi_y(t,rho_7,tau_f,tau);
    xyy = rho_8;
end

