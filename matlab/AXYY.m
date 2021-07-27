function axyy = AXYY(t,rho,tau_f,tau,taus)
global H_free H_phase
    % First pulse block
    rho_0 = time_evolution(H_free,taus(3),rho);
    rho_1 = time_evolution(H_phase(0),t,rho_0);
    rho_2 = time_evolution(H_free,taus(2),rho_1);
    rho_3 = time_evolution(H_phase(pi/6),t,rho_2);
    rho_4 = time_evolution(H_free,taus(1),rho_3);

     % XYXYYXY
     rho_5 = XYY(t,rho_4,tau_f,tau);
        
     % Third pulse block
     rho_12 = time_evolution(H_free,taus(1),rho_5);
     rho_13 = time_evolution(H_phase(pi/6),t,rho_12);
     rho_14 = time_evolution(H_free,taus(2),rho_13);
     rho_15 = time_evolution(H_phase(0),t,rho_14);
     rho_16 = time_evolution(H_free,taus(3),rho_15);
     rho_17 = time_evolution(H_phase(pi/2),t,rho_16);
     axyy = rho_17;
end

