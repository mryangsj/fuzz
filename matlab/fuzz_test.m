%% Discrete-time Implementantion of T-Bridge
%  Shijie Yang - 07/03/2024
clear 
close all
clc

%% Input Audio Signal
[Vin,Fs] = audioread('input_fuzz_static.wav');
Vin = Vin(:,1);

%% Number of Samples
Nsamp = length(Vin);

%% Simulated time
tstop = Nsamp/Fs;

%% User Parameters
alpha_lin = 0.5;
alpha_vol = 0.5;

%% Circuit Parameters
Rlin = 10e3;
Rvol = 50e3;
R1 = 23e3;
R2 = 12e3;
R3 = 220e3;
R4 = (1 - alpha_lin) * Rlin;
R5 = 22e3;
R6 = alpha_lin * Rlin;
Rin = 22e3;
C71 = 4.7e-6;
R72 = 50e3;
R81 = 22e3;
C81 = 4.7e-6;
C91 = 4.7e-6;
Rout = 100e3;
R92 = alpha_vol * Rvol;
R93 = (1 - alpha_vol) * Rvol;
C94 = 4.7e-6;
R10 = 1e-2;
Vcc = 9;

%% Adaptation Conditions
% r-type junction
Z1 = R1;
Z2 = R2;
Z3 = R3;
Z4 = R4;
Z5 = R5;
Z6 = R6;
Z10 = R10;

% port 7
Z_R71 = Rin;
Z_C71 = 1 / (2*Fs*C71);
Z_S71 = Z_R71 + Z_C71;
Z_R72 = R72;
Z_P72 = (Z_S71 * Z_R72) / (Z_S71 + Z_R72);
Z7 = Z_P72;

% port 8
Z_R81 = R81;
Z_C81 = 1 / (2*Fs*C81);
Z_P81 = (Z_R81 * Z_C81) / (Z_R81 + Z_C81);
Z8 = Z_P81;

% port 9
Z_R91 = Rout;
Z_C91 = 1 / (2*Fs*C91);
Z_S91 = Z_R91 + Z_C91;
Z_R92 = R92;
Z_P92 = (Z_R92 * Z_S91) / (Z_R92 + Z_S91);
Z_R93 = R93;
Z_S93 = Z_R93 + Z_P92;
Z_C94 = 1 / (2*Fs*C94);
Z_S94 = Z_C94 + Z_S93;
Z9 = Z_S94;

%% Solve R-type Junction
B = [ 1,  0,  0,  0,  0;
      0,  1,  0,  0,  0;
      0,  0,  1,  0,  0;
      0,  0,  0,  1,  0;
      0,  0,  0,  0,  1;
      0,  0, -1,  1,  0;
     -1,  1,  1,  0,  0;
      1,  0,  0,  0,  0;
      0, -1,  0,  0,  1;
      0,  0,  0, -1, -1;]';

syms Z11 Z12 Z21 Z22
Z_t = [Z11,Z12;Z21,Z22];
Z_b = diag([Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10]);
Z = sym(zeros(10,10));
Z(1:2,1:2) = Z_t;
Z(3:10,3:10) = Z_b;
S = eye(10) - 2*Z*B'/(B*Z*B')*B;

[Z11_sol,Z12_sol,Z21_zol,Z22_sol] = solve([S(1,1)==0,S(1,2)==0,S(2,1)==0,S(2,2)==0], [Z11,Z12,Z21,Z22]);
Z_t = double([Z11_sol,Z12_sol;Z21_zol,Z22_sol]);
Z = subs(Z, [Z11,Z12,Z21,Z22], [Z11_sol,Z12_sol,Z21_zol,Z22_sol]);
S = subs(S, [Z11,Z12,Z21,Z22], [Z11_sol,Z12_sol,Z21_zol,Z22_sol]);
Z = double(Z);
S = double(S);

S_R = (diag([R1,R2])/Z_t - eye(2)) / (diag([R1,R2])/Z_t + eye(2));

%% Solve Series and Parallel Junction
S_S71 = scattering_series([Z_S71, Z_C71, Z_R71]);
S_P72 = scattering_parallel([Z_P72, Z_R72, Z_S71]);
S_P81 = scattering_parallel([Z_P81, Z_C81, Z_R81]);
S_S91 = scattering_series([Z_S91, Z_R91, Z_C91]);
S_P92 = scattering_parallel([Z_P92, Z_R92, Z_S91]);
S_S93 = scattering_series([Z_S93, Z_R93, Z_P92]);
S_S94 = scattering_series([Z_S94, Z_C94, Z_S93]);

%% Initialize Output Signals
n_ports = 10;
Vout = zeros(size(Vin));
a = zeros(n_ports, 1);
b = zeros(n_ports, 1);
b_S71 = zeros(3, 1);
b_P72 = zeros(3, 1);
b_P81 = zeros(3, 1);
b_S91 = zeros(3, 1);
b_P92 = zeros(3, 1);
b_S93 = zeros(3, 1);
b_S94 = zeros(3, 1);

a_C71 = 0;
a_C81 = 0;
a_C91 = 0;
a_C94 = 0;

%% Initialization of Waves
unstable_flag = 1;
epsilon = 1e-16;
j = 0;
while (unstable_flag)
    j = j+1;

   %% Scattering

    % C and L
    b_C71 = a_C71;
    b_C81 = a_C81;
    b_C91 = a_C91;
    b_C94 = a_C94;

    % Forward
    b_S71(1) = S_S71(1,2) * b_C71;
    b_P72(1) = S_P72(1,3) * b_S71(1);
    a(7) = b_P72(1);

    b_P81(1) = S_P81(1,2) * b_C81;
    a(8) = b_P81(1);

    b_S91(1) = S_S91(1,3) * b_C91;
    b_P92(1) = S_P92(1,3) * b_S91(1);
    b_S93(1) = S_S93(1,3) * b_P92(1);
    b_S94(1) = S_S94(1,2) * b_C94 + S_S94(1,3) * b_S93(1);
    a(9) = b_S94(1);
    
    a(10) = Vcc;

    % Local
    b(1) = S(1,7)*a(7) + S(1,8)*a(8) + S(1,9)*a(9) + S(1,10)*a(10);
    b(2) = S(2,7)*a(7) + S(2,8)*a(8) + S(2,9)*a(9) + S(2,10)*a(10);
    a(1) = S_R(1,1)*b(1) + S_R(1,2)*b(2);
    a(2) = S_R(2,1)*b(1) + S_R(2,2)*b(2);

    % Backward
    b(7) = S(7,1)*a(1) + S(7,2)*a(2) + S(7,7)*a(7) + S(7,8)*a(8) + S(7,9)*a(9) + S(7,10)*a(10);
    b_P72(3) = S_P72(3,1) * b(7) + S_P72(3,3) * b_S71(1);
    b_S71(2) = S_S71(2,1) * b_P72(3) + S_S71(2,2) * b_C71;

    b(8) = S(8,1)*a(1) + S(8,2)*a(2) + S(8,7)*a(7) + S(8,8)*a(8) + S(8,9)*a(9) + S(8,10)*a(10);
    b_P81(2) = S_P81(2,1) * b(8) + S_P81(2,2) * b_C81;
    b_P81(3) = S_P81(3,1) * b(8) + S_P81(3,2) * b_C81;

    
    b(9) = S(9,1)*a(1) + S(9,2)*a(2) + S(9,7)*a(7) + S(9,8)*a(8) + S(9,9)*a(9) + S(9,10)*a(10);
    b_S94(2) = S_S94(2,1) * b(9) + S_S94(2,2) * b_C94 + S_S94(2,3) * b_S93(1);
    b_S94(3) = S_S94(3,1) * b(9) + S_S94(3,2) * b_C94 + S_S94(3,3) * b_S93(1);
    b_S93(3) = S_S93(3,1) * b_S94(3) + S_S93(3,3) * b_P92(1);
    b_P92(3) = S_P92(3,1) * b_S93(3) + S_P92(3,3) * b_S91(1);
    b_S91(3) = S_S91(3,1) * b_P92(3) + S_S91(3,3) * b_C91;


    if(abs(b_C71-b_S71(2))>epsilon || ...
       abs(b_C81-b_P81(2))>epsilon || ...
       abs(b_C91-b_S91(3))>epsilon || ...
       abs(b_C94-b_S94(2))>epsilon)
        a_C71 = b_S71(2);
        a_C81 = b_P81(2);
        a_C91 = b_S91(3);
        a_C94 = b_S94(2);
    else
        unstable_flag = 0;
    end
end

%% Simulation Loop
ii=0;
while (ii < Nsamp)
    ii = ii+1;
    
    %% Scattering

    % C and L
    b_C71 = a_C71;
    b_C81 = a_C81;
    b_C91 = a_C91;
    b_C94 = a_C94;

    % Forward
    b_S71(1) = S_S71(1,2) * b_C71 + S_S71(1,3) * Vin(ii);
    b_P72(1) = S_P72(1,3) * b_S71(1);
    a(7) = b_P72(1);

    b_P81(1) = S_P81(1,2) * b_C81;
    a(8) = b_P81(1);

    b_S91(1) = S_S91(1,3) * b_C91;
    b_P92(1) = S_P92(1,3) * b_S91(1);
    b_S93(1) = S_S93(1,3) * b_P92(1);
    b_S94(1) = S_S94(1,2) * b_C94 + S_S94(1,3) * b_S93(1);
    a(9) = b_S94(1);
    
    a(10) = 9;

    % Local
    b(1) = S(1,:)*a;
    b(2) = S(2,:)*a;
    a(1) = S_R(1,1)*b(1) + S_R(1,2)*b(2);
    a(2) = S_R(2,1)*b(1) + S_R(2,2)*b(2);

    % Backward
    b(7) = S(7,:)*a;
    b_P72(3) = S_P72(3,1) * b(7) + S_P72(3,3) * b_S71(1);
    b_S71(2) = S_S71(2,1) * b_P72(3) + S_S71(2,2) * b_C71 + S_S71(2,3) * Vin(ii);
    a_C71 = b_S71(2);

    b(8) = S(8,:)*a;
    b_P81(2) = S_P81(2,1) * b(8) + S_P81(2,2) * b_C81;
    b_P81(3) = S_P81(3,1) * b(8) + S_P81(3,2) * b_C81;
    a_C81 = b_P81(2);
    
    b(9) = S(9,:)*a;
    b_S94(2) = S_S94(2,1) * b(9) + S_S94(2,2) * b_C94 + S_S94(2,3) * b_S93(1);
    b_S94(3) = S_S94(3,1) * b(9) + S_S94(3,2) * b_C94 + S_S94(3,3) * b_S93(1);
    b_S93(3) = S_S93(3,1) * b_S94(3) + S_S93(3,3) * b_P92(1);
    b_P92(3) = S_P92(3,1) * b_S93(3) + S_P92(3,3) * b_S91(1);
    b_S91(2) = S_S91(2,1) * b_P92(3) + S_S91(2,3) * b_C91;
    b_S91(3) = S_S91(3,1) * b_P92(3) + S_S91(3,3) * b_C91;
    a_C91 = b_S91(3);
    a_C94 = b_S94(2);

    %% Read Output
    Vout(ii) = b_S91(2) / 2;

end

%% LTSpice Files
[Vout_LTSpice,Fs_LTspice] = audioread('output_fuzz_static.wav');
Vout_LTSpice = Vout_LTSpice(:,1);
N_LTspice = length(Vout_LTSpice);
time_LTSpice = (0:(N_LTspice-1)) / Fs_LTspice;

figure
set(gcf, 'Color', 'w');
plot(time_LTSpice,Vout_LTSpice,'r','Linewidth',2); hold on;
plot((0:Nsamp-1)/Fs,Vout,'b--','Linewidth',1); grid on;
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [Volt]','Fontsize',16,'interpreter','latex');
xlim([0,tstop]);
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
set(gca,'FontSize',15);
