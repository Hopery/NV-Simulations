% Creates a PiPulse Y
% Input: 
%   * `t` - Time for the time evolution. Usually 1/(2*Omega), which means ...
%   * `tau_f` - Tau times given by the fourier components
%   * `rho` - Density matrix `rho` (electron)
%   * `tau` - Time given for the loop of last function
% Returns: density matrix after sequence. For more info see 
% `Robust dynamical decoupling sequences for individual-nuclear-spin 
% addressing` paper.

function py = pi_y(t,rho,tau_f,tau)
global H_phase H_free
    % Taus
    tau_1 = (tau_f(1)*tau)-t/2;
    tau_2 = (tau_f(2)- tau_f(1))*tau-t;
    tau_3 = (tau_f(3)- tau_f(2))*tau-t;

    % Sequence
    rho_0 = time_evolution(H_free,tau_1,rho); 
    rho_1 = time_evolution(H_phase(pi/6 + pi/2), t, rho_0);      
    rho_2 = time_evolution(H_free,tau_2,rho_1); 
    rho_3 = time_evolution(H_phase(pi/2), t, rho_2);      
    rho_4 = time_evolution(H_free,tau_3,rho_3); 
    rho_5 = time_evolution(H_phase(pi/2 + pi/2), t, rho_4); 
    rho_6 = time_evolution(H_free,tau_3,rho_5); 
    rho_7 = time_evolution(H_phase(pi/2), t, rho_6); 
    rho_8 = time_evolution(H_free,tau_2,rho_7); 
    rho_9 = time_evolution(H_phase(pi/6 + pi/2), t, rho_8); 
    rho_10 = time_evolution(H_free,tau_1,rho_9); 
    
    py= rho_10;
end

