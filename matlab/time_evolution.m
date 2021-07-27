function U = time_evolution(H,t,rho)
% calculate the time evolution operator
U = expm(-1i * 2*pi*H * t);

U = U * rho * U';

end

