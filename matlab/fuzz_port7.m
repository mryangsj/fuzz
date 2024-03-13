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


%% Circuit Parameters
Rin = sym(22e3);
C71 = sym(4.7e-6);
R72 = sym(50e3);

%% Adaptation Conditions

% port 7
Z_R71 = Rin;
Z_C71 = Ts / (2*C71);
Z_S71 = Z_R71 + Z_C71;
Z_R72 = R72;
Z_P72 = (Z_S71 * Z_R72) / (Z_S71 + Z_R72);
Z7 = Z_P72;

%% Solve Series and Parallel Junction
S_S71 = double(scattering_series([Z_S71, Z_C71, Z_R71]));
S_P72 = double(scattering_parallel([Z_P72, Z_R72, Z_S71]));

%% Initialize Output Signals
n_ports = 10;
Vout = zeros(size(Vin));
a_S71 = zeros(3, 1);
b_S71 = zeros(3, 1);
a_P72 = zeros(3, 1);
b_P72 = zeros(3, 1);

a_C71 = 0;
Z7 = double(Z7);

%% Simulation Loop
ii=0;
while (ii < Nsamp)
    ii = ii+1;
    
    %% Scattering

    % C and L
    b_C71 = a_C71;

    % Forward
    b_S71(1) = S_S71(1,2) * b_C71 + S_S71(1,3) * Vin(ii);
    b_P72(1) = S_P72(1,3) * b_S71(1);
    a_R(7) = b_P72(1);

    % Local
    b = (100e3 - Z7) / (100e3 + Z7) * a_R(7);

    % Backward
    b_P72(3) = S_P72(3,1) * b + S_P72(3,3) * b_S71(1);
    b_S71(2) = S_S71(2,1) * b_P72(3) + S_S71(2,2) * b_C71 + S_S71(2,3) * Vin(ii);
    a_C71 = b_S71(2);

    %% Read Output
    Vout(ii) = -(a_R(7) + b) / 2;

end

%% LTSpice Files
[Vout_LTSpice,Fs_LTspice] = audioread('output_port7.wav');
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
