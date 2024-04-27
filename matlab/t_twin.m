%% Discrete-time Implementantion of T-Bridge
%  Shijie Yang - 07/03/2024
clear 
close all
clc

%% Input Audio Signal
[Vin, Fs] = audioread('input_tt.wav');
Vin = Vin(:,1);

%% Number of Samples
Nsamp = length(Vin);

%% Simulated time
tstop = Nsamp/Fs;

%% User Parameters


%% Circuit Parameters


%% Adaptation Conditions

% r-type port
Z1 = 1 / (2*Fs*31e-6);
Z2 = 1 / (2*Fs*71e-6);
Z3 = 32e3;
Z4 = 23e3;
Z5 = 6e-3;
Z7 = 1e3;
Z8 = 100e3;


%% Solve R-type Junction
B = [ 1,  0,  0,  0;
      0,  1,  0,  0;
      0,  0,  1,  0;
      0,  0,  0,  1;
      0,  0,  1,  1;
     -1,  0,  0, -1;
      1,  1,  0,  0;
      0, -1, -1,  0;]';

syms Z6
Z = diag([Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8]);
S = eye(8) - 2*Z*B'/(B*Z*B')*B;
S = subs(S, Z6, solve(S(6,6)==0, Z6));
S = double(S)

%% Solve Series and Parallel Junction


%% Initialize Output Signals
n_ports = 8;
Vout = zeros(size(Vin));
a = zeros(n_ports, 1);
b = zeros(n_ports, 1);

a_C1 = 0;
a_C2 = 0;
a_C3 = 0;
b_C1 = 0;
b_C2 = 0;
b_C3 = 0;

%% Initialization of Waves
unstable_flag = 1;
epsilon = 1e-16;
j = 0;
while (unstable_flag)
    j = j+1;

    %% Scattering

    % C and L
    b_C2 = a_C2;
    b_C3 = a_C3;

    % Forward
    a(1) = b_C2;
    a(2) = b_C3;
    a(5) = b_C1;
    a(7) = 9.0;

    % Local
    b(6) = S(6,1) * a(1) + S(6,2) * a(2) + S(6,5) * a(5) + S(6,7) * a(7);
    a(6) = -b(6);

    % Backward
    b(1) = S(1,1) * a(1) + S(1,2) * a(2) + S(1,5) * a(5) + S(1,7) * a(7);
    b(2) = S(2,1) * a(1) + S(2,2) * a(2) + S(2,5) * a(5) + S(2,7) * a(7);
    b(5) = S(5,1) * a(1) + S(5,2) * a(2) + S(5,5) * a(5) + S(5,7) * a(7);

    if(abs(b_C1-b(5))>epsilon || abs(b_C2-b(1))>epsilon || abs(b_C3-b(2))>epsilon)
        a_C1 = b(5);
        a_C2 = b(1);
        a_C3 = b(2);
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
    b_C1 = a_C1;
    b_C2 = a_C2;
    b_C3 = a_C3;

    % Forward
    a(1) = b_C2;
    a(2) = b_C3;
    a(3) = 0;
    a(4) = 0;
    a(5) = b_C1;
    a(7) = 0;
    a(8) = 0;
    b(6) = S(6,:) * a;

    % Local
    a(6) = 2*Vin(ii) - b(6);

    % Backward
    b(1) = S(1,:) * a;
    b(2) = S(2,:) * a;
    b(3) = S(3,:) * a;
    b(4) = S(3,:) * a;
    b(5) = S(5,:) * a;
    b(7) = S(7,:) * a;
    b(8) = S(8,:) * a;

    a_C1 = b(5);
    a_C2 = b(1);
    a_C3 = b(2);

    %% Read Output
    Vout(ii) = (a(8)+b(8)) / 2;

end

%% LTSpice Files
[Vout_LTSpice,Fs_LTspice] = audioread('output_tt.wav');
Vout_LTSpice = Vout_LTSpice(:,1);
N_LTspice = length(Vout_LTSpice);
time_LTSpice = (0:(N_LTspice-1)) / Fs_LTspice;

figure
% set(gcf, 'Color', 'w');
plot(time_LTSpice,Vout_LTSpice,'r','Linewidth',2); hold on; 
plot((0:Nsamp-1)/Fs,Vout,'b--','Linewidth',1); grid on;
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [Volt]','Fontsize',16,'interpreter','latex');
xlim([0,tstop]);
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
set(gca,'FontSize',15);
