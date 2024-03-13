%% Discrete-time Implementantion of ETI 1 Transistor Fuzz
%  Shijie Yang - 08/03/2024
clear 
close all
clc

%% Input Audio Signal
[Vin, Fs] = audioread('input_fuzz.wav');
Vin = Vin(:,1);

%% Sampling Period
Ts = 1/Fs;

%% Number of Samples
Nsamp = length(Vin);

%% Simulated time
tstop = Nsamp*Ts;

%% User Parameters


%% Circuit Parameters


%% Adaptation Conditions

% r-type port
Z3 = Ts / (2*11.2e-9);
Z4 = Ts / (2*11.2e-9);
Z5 = 20e3;
Z6 = 100e3;
Z7 = 10e3;

%% Solve R-type Junction
Q = [ 1,  0,  0;
      0,  1,  0;
     -1,  0,  1;
      0,  1, -1;
     -1,  1,  0;
      0,  1,  0;
      0,  0,  1;]';

syms Z11 Z12 Z21 Z22

Z_t = [Z11, Z12; Z21, Z22];
Z_b = diag([Z3,Z4,Z5,Z6,Z7]);

Z = sym(zeros(7,7));
Z(1:2,1:2) = Z_t;
Z(3:7,3:7) = Z_b;

S = 2*Q' /(Q/Z*Q') * Q/Z - eye(7);

[Z11_sol,Z12_sol,Z21_sol,Z22_sol] = solve([S(1,1)==0,S(1,2)==0,S(2,1)==0,S(2,2)==0], [Z11,Z12,Z21,Z22]);
S = subs(S, [Z11,Z12,Z21,Z22], [Z11_sol,Z12_sol,Z21_sol,Z22_sol]);

Z11 = double(Z11_sol);
Z12 = double(Z12_sol);
Z21 = double(Z21_sol);
Z22 = double(Z22_sol);
S = double(S)

%% Solve Series and Parallel Junction


%% Initialize Output Signals
n_ports = 7;
Vout = zeros(size(Vin));
a = zeros(n_ports, 1);
b = zeros(n_ports, 1);

a_Cm = 0;
a_Ch = 0;
A0 = 100e3;

%% Simulation Loop
ii=0;
while (ii < Nsamp)
    ii = ii+1;
    
    %% Scattering

    % C and L
    b_Cm = a_Cm;
    b_Ch = a_Ch;

    % Forward
    a(3) = b_Cm;
    a(4) = b_Ch;
    a(5) = 0;
    a(6) = 0;
    a(7) = Vin(ii);

    % Local
    b(1) = S(1,:) * a;
    b(2) = S(2,:) * a;
    a(1) = (Z22+A0*Z12)/(Z22-A0*Z12)*b(1) + -2*Z12/(Z22-A0*Z12)*b(2);
    a(2) = 2*A0*Z22/(Z22-A0*Z12)*b(1) + (-Z22-A0*Z12)/(Z22-A0*Z12)*b(2);

    % Backward
    b(3) = S(3,:) * a;
    b(4) = S(4,:) * a;
    b(5) = S(5,:) * a;
    b(6) = S(6,:) * a;
    b(7) = S(7,:) * a;

    a_Cm = b(3);
    a_Ch = b(4);

    %% Read Output
    Vout(ii) = (a(6)+b(6)) / 2;

end

%% LTSpice Files
[Vout_LTSpice,Fs_LTspice] = audioread('output_abp.wav');
Vout_LTSpice = Vout_LTSpice(:,1);
N_LTspice = length(Vout_LTSpice);
time_LTSpice = (0:(N_LTspice-1)) / Fs_LTspice;

figure
% set(gcf, 'Color', 'w');
plot(time_LTSpice,Vout_LTSpice,'r','Linewidth',2); hold on; 
plot(double(Ts)*[0:Nsamp-1],Vout,'b--','Linewidth',1); grid on;
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [Volt]','Fontsize',16,'interpreter','latex');
xlim([0,tstop]);
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
set(gca,'FontSize',15);
