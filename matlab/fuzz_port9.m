%% Discrete-time Implementantion of the ETI Fuzz
%  Shijie Yang - 04/03/2024
clear 
close all
clc

%% Input Audio Signal
[Vin,Fs] = audioread('input.wav');
Vin = Vin(:,1);

%% Sampling Period
Ts = sym(1/Fs);

%% Number of Samples
Nsamp = length(Vin);

%% Simulated time
tstop = Nsamp*Ts;

%% User Parameters
alpha_vol = sym(0.5);

%% Circuit Parameters
C91 = sym(4.7e-6);
Rout = sym(100e3);
R92 = sym(alpha_vol * 50e3);
R93 = sym((1 - alpha_vol) * 50e3);
C94 = sym(4.7e-6);

%% Adaptation Conditions

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

%% Solve Series and Parallel Junction
S_S91 = double(scattering_series([Z_S91, Z_R91, Z_C91]));
S_P92 = double(scattering_parallel([Z_P92, Z_R92, Z_S91]));
S_S93 = double(scattering_series([Z_S93, Z_R93, Z_P92]));
S_S94 = double(scattering_series([Z_S94, Z_C94, Z_S93]));

%% Initialize Output Signals
n_ports = 10;
Vout = zeros(size(Vin));
a_S91 = zeros(3, 1);
b_S91 = zeros(3, 1);
a_P92 = zeros(3, 1);
b_P92 = zeros(3, 1);
a_S93 = zeros(3, 1);
b_S93 = zeros(3, 1);
a_S94 = zeros(3, 1);
b_S94 = zeros(3, 1);

a_C91 = 0;
a_C94 = 0;

%% Simulation Loop
ii=0;
while (ii < Nsamp)
    ii = ii+1;
    
    %% Scattering

    % C and L
    b_C91 = a_C91;
    b_C94 = a_C94;

    % Forward
    b_S91(1) = S_S91(1,3) * b_C91;
    b_P92(1) = S_P92(1,3) * b_S91(1);
    b_S93(1) = S_S93(1,3) * b_P92(1);
    b_S94(1) = S_S94(1,2) * b_C94 + S_S94(1,3) * b_S93(1);
    a_R(9) = b_S94(1);

    % Local
    b = 2*Vin(ii) - a_R(9);

    % Backward
    b_S94(2) = S_S94(2,1) * b + S_S94(2,2) * b_C94 + S_S94(2,3) * b_S93(1);
    b_S94(3) = S_S94(3,1) * b + S_S94(3,2) * b_C94 + S_S94(3,3) * b_S93(1);
    b_S93(3) = S_S93(3,1) * b_S94(3) + S_S93(3,3) * b_P92(1);
    b_P92(3) = S_P92(3,1) * b_S93(3) + S_P92(3,3) * b_S91(1);
    b_S91(2) = S_S91(2,1) * b_P92(3) + S_S91(2,3) * b_C91;
    b_S91(3) = S_S91(3,1) * b_P92(3) + S_S91(3,3) * b_C91;
    a_C91 = b_S91(3);
    a_C94 = b_S94(2);

    %% Read Output
    Vout(ii) = -b_S91(2) / 2;

end

%% LTSpice Files
[Vout_LTSpice,Fs_LTspice] = audioread('output_port9.wav');
Vout_LTSpice = Vout_LTSpice(:,1);
N_LTspice = length(Vout_LTSpice);
time_LTSpice = (0:(N_LTspice-1)) / Fs_LTspice;

figure
set(gcf, 'Color', 'w');
plot(time_LTSpice,Vout_LTSpice,'r','Linewidth',2); hold on; 
plot(double(Ts)*[0:Nsamp-1],Vout,'b--','Linewidth',1); grid on;
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [Volt]','Fontsize',16,'interpreter','latex');
xlim([0,10]);
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
set(gca,'FontSize',15);
