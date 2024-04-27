%% Discrete-time Implementantion of T-Twin
%  Shijie Yang - 07/03/2024
clear 
close all
clc

%% Input Audio Signal
[Vin, Fs] = audioread('input_t.wav');
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
Z1 = Ts / (2*(13e-6));
Z2 = Ts / (2*(31e-6));
Z3 = 5e3;
Z5 = 10;
Z6 = 100e3;


%% Solve R-type Junction
B = [ 1,  0,  0;
      0,  1,  0;
      0,  0,  1;
     -1,  0, -1;
     -1,  1,  0;
      0,  1,  1;]';

syms Z4
Z = diag([Z1,Z2,Z3,Z4,Z5,Z6]);
S = eye(6) - 2*Z*B'/(B*Z*B')*B;
S = subs(S, Z4, solve(S(4,4)==0, Z4));
S = double(S)

%% Solve Series and Parallel Junction


%% Initialize Output Signals
n_ports = 6;
Vout = zeros(size(Vin));
a = zeros(n_ports, 1);
b = zeros(n_ports, 1);

a_C1 = 0;
a_C2 = 0;

%% Simulation Loop
ii=0;
while (ii < Nsamp)
    ii = ii+1;
    
    %% Scattering

    % C and L
    b_C1 = a_C1;
    b_C2 = a_C2;

    % Forward
    a(1) = b_C1;
    a(2) = b_C2;
    a(3) = 0;
    a(5) = 0;
    a(6) = 0;
    b(4) = S(4,:) * a;

    % Local
    a(4) = 2*Vin(ii) - b(4);

    % Backward
    b(1) = S(1,:) * a;
    b(2) = S(2,:) * a;
    b(3) = S(3,:) * a;
    b(5) = S(5,:) * a;
    b(6) = S(6,:) * a;

    a_C1 = b(1);
    a_C2 = b(2);

    %% Read Output
    Vout(ii) = (a(6)+b(6)) / 2;

end

%% LTSpice Files
[Vout_LTSpice,Fs_LTspice] = audioread('output_t.wav');
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
