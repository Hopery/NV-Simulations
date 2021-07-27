close all;
clear;

% ----------------------------------------------------------------------
% Operators

global S_x S_y S_z I_x I_y I_z 

% Electron operators
S_x = 1/2 * [0 1;1 0];
S_x = kron(S_x, eye(2));
S_y = 1/(2*1i) * [0 1;-1 0];
S_y = kron(S_y, eye(2));
S_z = [0 0;0 -1];
S_z = kron(S_z, eye(2));

% Nucleous operators
I_x = 1/2 * [0 1;1 0];
I_x = kron(eye(2),I_x);
I_y = 1/(2*1i) * [0 1;-1 0];
I_y = kron(eye(2),I_y);
I_z = 1/2 * [1 0; 0 -1];
I_z = kron(eye(2),I_z);

% Density matrix
rho_nv = [1 0; 0 0];
rho_13c = [0.5 0; 0 0.5];
rho = kron(rho_nv, rho_13c);

% ----------------------------------------------------------------------
% Parameters

global gamma_13c B D

% Energy between spin state 0 and +-1 (zfs in ground state)(MHz)
D = 2870;
% Magnetic field (G)
B = 500;
% gyromagnetic ratio of the NV center’s electron spin(MHz/G)
gamma_nv = 2.8;
gamma_13c = 1.07e-3;
%10e-3
Omega = 40;
Omega_13c = 10e-3;
% omega electron (MHz)
w = D - gamma_nv * B;
% omega nucleous (MHz)
w_l = gamma_13c*B;

global A_zz A_zx 

% Hyperfine
A_zx = 50e-3;
A_zy = 0;
A_zz = 10e-3;

% ----------------------------------------------------------------------
% Hamiltonians

% To get a Hamiltonian with the energy levels in the diagonal (zfs = zero
% field splitting)
H_zfs = D * S_z^2;

% Interactions between electron and nucleous. Zeeman= when +1, -1 and 0 are
% splittig, if not -1 and +1 are in the same position.
H_zeeman_nv = gamma_nv * B * S_z;
H_zeeman_13c = gamma_13c * B * I_z;

% To rotate in X or Y axis, Omega bc is the "line of the wave"
H_drive_x = Omega * S_x;
H_drive_y = Omega * S_y;
H_drive_x_13c = Omega_13c * I_x;
H_drive_y_13c = Omega_13c * I_y;

% Hyperfine Hamiltonian. Describes hyperfine interaction between the
% NV centre’s electron spin and the nitrogen’s nuclear spin
H_hyperfine = A_zx * S_z * I_x + A_zy * S_z * I_y + A_zz * S_z * I_z;

% w means the different between state 1 and -1 when a magnetic flied is
% aplied, we multiply it by same HAmiltonian as H_zfs to represent the
% different spin states.
H_drive_freq = - w * S_z^2;
H_drive_freq_13c = - w_l * I_z;
% -----------------------------------------------------------------------
% Hamiltonian of the system

global H H_free H_x H_y H_phase H_13c H_13c_w

% Hamiltonian for the nv nucleus
H_13c = H_zeeman_13c + H_drive_freq_13c + H_drive_x_13c + H_hyperfine;
% Hamiltonian for the nv nucleus with frequency (omega) as parameter
H_13c_w = @(omega) H_zeeman_13c - omega*I_z + H_drive_x_13c + H_hyperfine;

% Full electron hamiltonian with both drives
H = H_zfs + H_zeeman_nv + H_zeeman_13c + H_drive_x + H_drive_y + ...
    H_hyperfine + H_drive_freq;

% Electron hamiltonian with no drives
H_free = H_zfs + H_zeeman_nv + H_zeeman_13c + H_hyperfine ...
    + H_drive_freq;

% Electron hamiltonian with one drive
H_x = H_zfs + H_zeeman_nv + H_zeeman_13c + H_drive_x ...
    + H_hyperfine + H_drive_freq;

H_y = H_zfs + H_zeeman_nv + H_zeeman_13c + H_drive_y ...
    + H_hyperfine + H_drive_freq;

% Electron hamiltonian with phase as parameter (so can be x or y with some
% phase)
H_phase = @(phase) H_zfs + H_zeeman_nv + H_zeeman_13c + ...
    kron([0 (Omega*exp(1i*phase))/2;(Omega*exp(-1i*phase))/2 0], eye(2)) ...
    + H_hyperfine + H_drive_freq;


%------------------------------------------------------------------------
% Fourier components -> taus

global f1e f2e f3e f4e

f1e = 0.2;
f2e = 0;
f3e = 0;
f4e = 0;

theta1=0.1*pi;
theta2=0.3*pi;
theta3=0.6*pi;
theta4=0.9*pi; 
x0=[theta1, theta2,theta3,theta4];

% ------------------------------------------------------------------------
% rho CNOT

rho_nv = [1 0; 0 0];
rho_13c = 0.5*[1 -1; -1 1];
rho_CNOT = kron(rho_nv, rho_13c);

% -----------------------------------------------------------------------
%% Run
% rabi(rho)
% ramsey(1/(40)/4,rho)
SpinEcho(1/(2*Omega), 10, rho)
% SpinEchoN(1/(2*Omega), 100, rho)
% SpinEchoFourier(x0,1/(2*Omega),rho,16);
% times = AXYXDispl(x0,1/(2*Omega),rho,100);
% disp(rho_CNOT)
% times = 27;
% CNOT(x0,rho_CNOT,times/2,1/(2*Omega))
% AXYYDispl(x0,1/(2*Omega),rho,100);
%disp(PartialTrace(rho,2,[2,2]));
% disp(rho)
%SWAP(x0,rho,times/2,1/(2*Omega));
% RFandSWAP(x0,rho,times/2,1/(2*Omega))
% odmrRFSWAP(x0,rho,times/2,1/(2*Omega))