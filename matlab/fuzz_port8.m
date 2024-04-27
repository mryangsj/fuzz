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
R81 = sym(22e3);
C81 = sym(4.7e-6);

%% Adaptation Conditions

% port 8
Z_R81 = R81;
Z_C81 = Ts / (2*C81);
Z_P81 = (Z_R81 * Z_C81) / (Z_R81 + Z_C81);
Z8 = Z_P81;

%% Solve Series and Parallel Junction
S_P81 = double(scattering_parallel([Z_P81, Z_C81, Z_R81]));

%% Initialize Output Signals
n_ports = 10;
Vout = zeros(size(Vin));
a_P81 = zeros(3, 1);
b_P81 = zeros(3, 1);

a_C81 = 0;

%% Simulation Loop
ii=0;
while (ii < Nsamp)
    ii = ii+1;
    
    %% Scattering

    % C and L
    b_C81 = a_C81;

    % Forward
    b_P81(1) = S_P81(1,2) * b_C81;
    a_R(8) = b_P81(1);

    % Local
    b = 2*Vin(ii) - a_R(8);

    % Backward
    b_P81(2) = S_P81(2,1) * b + S_P81(2,2) * b_C81;
    a_C81 = b_P81(2);
    b_P81(3) = S_P81(3,1) * b + S_P81(3,2) * b_C81;

    %% Read Output
    Vout(ii) = b_P81(3) / 2;

end

%% LTSpice Files
[Vout_LTSpice,Fs_LTspice] = audioread('output_port8.wav');
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
