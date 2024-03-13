%% Discrete-time Implementantion of the ETI Fuzz
%  Shijie Yang - 04/03/2024
clear 
close all
clc

%% Input Audio Signal
[Vin,Fs] = audioread('ManoucheGuitarCompr.wav');

%% Sampling Period
Ts = sym(1/Fs);

%% Number of Samples
Nsamp = length(Vin);

%% Simulated time
tstop = Nsamp*Ts;

%% User Parameters
alpha_lin = sym(0.5);
alpha_vol = sym(0.8);

%% Circuit Parameters 
R3 = sym(220e3);
R4 = sym((1-alpha_lin) * 10e3 + 1e-3);
R5 = sym(22e3);
R6 = sym(alpha_lin * 10e3 + 1e-3);
Rin = sym(22e3);
C71 = sym(4.7e-6);
R72 = sym(50e3);
R81 = sym(22e3);
C81 = sym(4.7e-6);
C91 = sym(4.7e-6);
Rout = sym(100e3);
R92 = sym(alpha_vol * 50e3 + 1e-3);
R93 = sym((1 - alpha_vol) * 50e3 + 1e-3);
C94 = sym(4.7e-6);
R10 = sym(1e-3);
Vcc = sym(9);

%% Adaptation Conditions

% r-type port
Z3 = R3;
Z4 = R4;
Z5 = R5;
Z6 = R6;
Z10 = R10;

% port 7
Z_R71 = Rin;
Z_C71 = Ts / (2*C71);
Z_S71 = Z_R71 + Z_C71;
Z_R72 = R72;
Z_P72 = (Z_S71 * Z_R72) / (Z_S71 + Z_R72);
Z7 = Z_P72;

% port 8
Z_R81 = R81;
Z_C81 = Ts / (2*C81);
Z_P81 = (Z_R81 * Z_C81) / (Z_R81 + Z_C81);
Z8 = Z_P81;

% port 9
Z_R91 = Rout;
Z_C91 = Ts / (2*C91);
Z_S91 = Z_R91 + Z_C91;
Z_R92 = R92;
Z_P92 = (Z_R92 * Z_S91) / (Z_R92 + Z_S91);
Z_R93 = R93;
Z_S93 = Z_R93 + Z_P92;
Z_C94 = Ts / (2*C94);
Z_S94 = Z_C94 + Z_S93;
Z9 = Z_S94;

%% Solve R-type Junction
B_R = [ 1,  0,  0,  0,  0;
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
Z_t = [Z11, Z12; Z21, Z22];
Z_b = diag([Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10]);
Z_R = sym(zeros(10, 10));
Z_R(1:2, 1:2) = Z_t;
Z_R(3:10, 3:10) = Z_b;
S_R = eye(10) - 2*Z_R*B_R'/(B_R*Z_R*B_R')*B_R;
[Z11, Z12, Z21, Z22] = solve([S_R(1,1)==0,S_R(1,2)==0,S_R(2,1)==0,S_R(2,2)==0], [Z11,Z12,Z21,Z22]);
Z_R(1:2, 1:2) = [Z11,Z12;Z21,Z22];
S_R = eye(10) - 2*Z_R*B_R'/(B_R*Z_R*B_R')*B_R;
S_R = double(S_R);

%% Solve Series and Parallel Junction
S_S71 = double(scattering_series([Z_S71, Z_C71, Z_R71]));
S_P72 = double(scattering_parallel([Z_P72, Z_R72, Z_S71]));
S_P81 = double(scattering_parallel([Z_P81, Z_C81, Z_R81]));
S_S91 = double(scattering_series([Z_S91, Z_R91, Z_C91]));
S_P92 = double(scattering_parallel([Z_P92, Z_R92, Z_S91]));
S_S93 = double(scattering_series([Z_S93, Z_R93, Z_P92]));
S_S94 = double(scattering_series([Z_S94, Z_C94, Z_S93]));

%% Initialize Output Signals
n_ports = 10;
Vout = zeros(size(Vin));
a_R = zeros(n_ports, 1);
b_R = zeros(n_ports, 1);
a_S71 = zeros(3, 1);
b_S71 = zeros(3, 1);
a_P72 = zeros(3, 1);
b_P72 = zeros(3, 1);
a_P81 = zeros(3, 1);
b_P81 = zeros(3, 1);
a_S91 = zeros(3, 1);
b_S91 = zeros(3, 1);
a_P92 = zeros(3, 1);
b_P92 = zeros(3, 1);
a_S93 = zeros(3, 1);
b_S93 = zeros(3, 1);
a_S94 = zeros(3, 1);
b_S94 = zeros(3, 1);

a_C71 = 0;
a_C81 = 0;
a_C91 = 0;
a_C94 = 0;
a_R(10) = 9;

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
    b_S71(1) = S_S71(1,:) * [b_P72(3); b_C71; Vin(ii)];
    b_P72(1) = S_P72(1,:) * [b_R(7); 0; b_S71(1)];
    a_R(7) = b_P72(1);

    b_P81(1) = S_P81(1,:) * [b_R(8); b_C81; 0];
    a_R(8) = b_P81(1);

    b_S91(1) = S_S91(1,:) * [b_P92(3); 0; b_C91];
    b_P92(1) = S_P92(1,:) * [b_S93(3); 0; b_S91(1)];
    b_S93(1) = S_S93(1,:) * [b_S94(3); 0; b_P92(1)];
    b_S94(1) = S_S94(1,:) * [b_R(9); b_C94; b_S93(1)];
    a_R(9) = b_S94(1);

    % Local
    b_R(1) = S_R(1,:) * a_R;
    b_R(2) = S_R(2,:) * a_R;
    a_R(1) = 0;
    a_R(2) = 0;

    % Backward
    b_R(7) = S_R(7,:) * a_R;
    b_R(8) = S_R(8,:) * a_R;
    b_R(9) = S_R(9,:) * a_R;

    b_P72(3) = S_P72(3,:) * [b_R(7); 0; b_S71(1)];
    b_S71(2) = S_S71(2,:) * [b_P72(3); b_C71; Vin(ii)];
    a_C71 = b_S71(2);

    b_P81(2) = S_P81(2,:) * [b_R(8); b_C81; 0];
    a_C81 = b_P81(2);

    b_S94(2) = S_S94(2,:) * [b_R(9); b_C94; b_S93(1)];
    b_S94(3) = S_S94(3,:) * [b_R(9); b_C94; b_S93(1)];
    b_S93(3) = S_S93(3,:) * [b_S94(3); 0; b_P92(1)];
    b_P92(3) = S_P92(3,:) * [b_S93(3); 0; b_S91(1)];
    b_S91(2) = S_S91(2,:) * [b_P92(2); 0; b_C91];
    b_S91(3) = S_S91(3,:) * [b_P92(3); 0; b_C91];
    a_C94 = b_S94(2);
    a_C91 = b_S91(3);

    %% Read Output
    Vout(ii) = b_S91(2) / 2;

end

%% LTSpice Files
load('Vout_LTspice.txt');
timeLTSpice = Vout_LTspice(:,1);
VoutLTSpice = Vout_LTspice(:,2);

figure
set(gcf, 'Color', 'w');
plot(timeLTSpice,VoutLTSpice,'r','Linewidth',2); hold on; 
plot(Ts*[1:Nsamp],Vout,'b--','Linewidth',1); grid on;
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [Volt]','Fontsize',16,'interpreter','latex');
xlim([0,8]);
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
set(gca,'FontSize',15);
