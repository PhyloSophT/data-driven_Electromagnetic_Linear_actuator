clear
close all
clc
ds_factor = 1;
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
dx = squeeze(dx_sim);
ddx_sim = simOut.get('ddx_sim');     % Zeitvektor
ddx = squeeze(ddx_sim);

%Berechnung der Geschwindikeit sowie der Beschleunigung
%dx1 = TVRegDiff(x,Ts,10);
%dx1(1,:)=[];
%dx1 = dx1*1000;
%dx=dx1;

tv_dx = TVRegDiff(x,1,1,[],[],1e-6,1,[],[]);
tv_dx(1,:)=[];
tv_dx = tv_dx/Ts;

tv_ddx = TVRegDiff(tv_dx,1,1,[],[],1e-6,1,[],[]);
tv_ddx(1,:)=[];
tv_ddx = tv_ddx/Ts;
size(t)
size(tv_dx)
size(tv_ddx)
%dx  = gradient(x,Ts);     % Geschwinigkeit
%ddx = gradient(dx,Ts);   % Beschleunigung

%ddx1 = TVRegDiff(dx1,Ts,10);
%ddx1(1,:)=[];
%ddx1 = ddx1*1000;
%ddx=ddx1;




%% Radial Basis Function RBF-Approximation einer Nichtlinearität nonlinearity alpha*x^3
u = dx;
y = F_R*atan(M*u);
%plot(u,y)
%grid on
p = 31;
u_max = max(u);
u_min = min(u);
n_y = length(y);
sigma_n = 0.7;
del_e = (u_max-u_min)/(p-1);
sigma = sigma_n*del_e;

C_matrix = zeros(n_y,p);
u_m_vec = zeros(p,1);
for i = 1:p
u_m = u_min+(i-1)*(u_max-u_min)/(p-1);
C_i = (u-u_m).^2;
C_matrix(:,i) = C_i;
u_m_vec(i) = u_m;
end
A_RBF = exp(-C_matrix/(2*sigma^2));
A_RBF_sum = sum(A_RBF,2);
A_RBF_norm = A_RBF./A_RBF_sum;
Psi_RBF = A_RBF_norm;
%figure(2)
%plot(u,A_RBF_norm)
%hold on

% Extrahiere Signale
x1 = x(:);        % Position
v = dx(:,1);      % Geschwindigkeit
a = ddx(:,1);      % Beschleunigung (Zielgröße)

%% === Least Squares Ansatz ===
% a = - (D/m)*v  + (K/m)*u -(F_R/m)*arctan(M*v)
% => y = Phi * theta, mit:
%    y = a
%    Phi = [-v, u,RBF_netze]
%    theta = [D/m; K/m, ... ]

Phi = [-v,I , -Psi_RBF]; %%%%falls bei der identifikation mehrere therme von x, xpunkt oder xpunktpunkt müssen diese zusammengefasst werden und können nicht mehr einzeln betrachtet werden
Y = a;
% Least-Squares-Schätzung
theta_hat = Phi \ Y;

% Parameter rekonstruieren
theta_hat = theta_hat*m;

D_hat = theta_hat(1,1);
K_hat = theta_hat(2,1);
theta_FR = theta_hat(3:end,1);
FR_RBF = Psi_RBF*theta_FR;
figure
plot(u,FR_RBF+D_hat*u, 'DisplayName', 'Rekonstruktion')
hold on
plot(u,y+D*u, 'DisplayName', 'Reibkennlinie')
xlabel('Zeit [s]');
legend; grid on;
title('Reibkennlinie und Rekonstruktion');

set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontSize', 16);
% Plot: Position
figure;
plot(t, I, '--g', 'DisplayName', 'u(t) (Simulink)');
hold on;
plot(t, x(1:ds_factor:end), '--y', 'DisplayName', 'x(t) (Simulink)');
hold on;
plot(t, dx, '--b', 'DisplayName', 'v(t) (Simulink)');
hold on;
plot(t, ddx, '--r', 'DisplayName', 'a(t) (Simulink)');
hold on;
xlabel('Zeit [s]'); 
legend; grid on;
title('Linearaktor');
