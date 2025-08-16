clear
close all
clc
%% Parameter für CHIRP-Signal
fs = 1000;         % Abtastfrequenz [Hz]
T = 10;            % Gesamtdauer [s]
t = (0:1/fs:T)';   % Zeitvektor

% Chirp-Signal
f0 = 1;            % Startfrequenz [Hz]
f1 = 10;          % Endfrequenz [Hz]
chirp_signal = chirp(t, f0, T, f1);

% Im Workspace speichern
assignin('base', 't_chirp', t);
assignin('base', 'chirp_signal', chirp_signal);

%% Parameter für Modell
% Parameter, die nicht bestimmt werden müssen
m = 0.5;
M = 100;
F_R = 1;
x0 = [0 ; 0];

%Parameter, die bestimmt werden sollen
K = 3;
D = 0.5;

% Abtastzeit
Ts = 0.001;

%% === Simulink-Modell starten ===
model_name = 'nichtlineares_ZRM_sim';
simOut = sim(model_name, 'ReturnWorkspaceOutputs', 'on');

%% === Daten extrahieren ===
x = simOut.get('x');     % x(t) [Nx1] – Position
I = simOut.get('i');     % Eingang u(t)
t = simOut.get('t');     % Zeitvektor
dx_sim = simOut.get('dx_sim');     % Zeitvektor
ddx_sim = simOut.get('ddx_sim');     % Zeitvektor
%% Verschiedene Ableitungsverfahren

% Aus der Simulation
dx1 = squeeze(dx_sim);
ddx1 = squeeze(ddx_sim);

% Gradienten Verfahren
dx2 = gradient(x,Ts);
ddx2 = gradient(dx2,Ts);

%TVRegdiff Verfahren
dx3 = TVRegDiff(x,Ts,10);
dx3(1,:)=[];
dx3 = dx3*1000;

ddx3 = TVRegDiff(dx1,Ts,10);
ddx3(1,:)=[];
ddx3 = ddx3*1000;

%% RBF-Netze zur jeder Differntiation
[FR_RBF1,~,K_hat1,D_hat1] = Function_RBF(dx1,ddx1,F_R,M,m);
[FR_RBF2,~,K_hat2,D_hat2]= Function_RBF(dx2,ddx2,F_R,M,m);
[FR_RBF3,~,K_hat3,D_hat3]= Function_RBF(dx3,ddx3,F_R,M,m);

u1 = dx1;
u2 = dx2;
u3 = dx3;
y = F_R*atan(M*u1);

figure
plot(u1,FR_RBF1+D_hat1*u1)
hold on
plot(u2,FR_RBF2+D_hat2*u2)
hold on
plot(u3,FR_RBF3+D_hat3*u3)
hold on
plot(u1,y+D*u1)
hold on 
