% Auxiliar function for `SpinEcho`. Here the sequence XYXYYXYX is applied
% the desired number of times.
% Input: 
%   * `t` - Time for the time evolution.
%   * `times` - Times the secuence XYXYYXYX will be aplied
%   * `rho` - Density matrix `rho` (electron)
%   * `tau` - Time given for the loop of last function
% Returns: density matrix after sequence. This sequence is usefull for
% decopling electro from enviroment and nucleus (less noise)

function a = SpinEchoAux(times,t,tau,rho)
global H_x H_y H_free 
    for i = 1:times
            rho_1 = time_evolution(H_x, t,rho);% x 1
            rho_2 = time_evolution(H_free,tau,rho_1); % tau
            rho_3 = time_evolution(H_y, t,rho_2); % y 2
            rho_4 = time_evolution(H_free,tau,rho_3); % tau
            rho_5 = time_evolution(H_x, t,rho_4); % x 3
            rho_6 = time_evolution(H_free,tau,rho_5); % tau
            rho_7 = time_evolution(H_y, t,rho_6); % y 4
            rho_8 = time_evolution(H_free,tau,rho_7); % tau
            rho_9 = time_evolution(H_y, t,rho_8); % y 5
            rho_10 = time_evolution(H_free,tau,rho_9); % tau
            rho_11 = time_evolution(H_x, t,rho_10); % x 6
            rho_12 = time_evolution(H_free,tau,rho_11); % tau
            rho_13 = time_evolution(H_y, t,rho_12); % y 7
            rho_14 = time_evolution(H_free,tau,rho_13); % tau
            rho_15 = time_evolution(H_x, t,rho_14); % x 8
            rho_16 = time_evolution(H_free,0.5*tau,rho_15); % tau/2
            rho = rho_16;
    end
    a = rho;
end