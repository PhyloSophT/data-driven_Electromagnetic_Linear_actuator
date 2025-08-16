
clear
close all
clc
ds_factor = 1;
m = 0.5;
D = 0.2;
K = 2;
C = 0;

Ts = 0.001;

%% Zustandsraumdarstellung
A = [0 1; -(C/m) -(D/m)];
b = [0 ; K/m];
c = [1 0];
d = 0;
x0 = [0 ; 0];

% Parameter für CHIRP-Signal
fs = 1000;         % Abtastfrequenz [Hz]
T = 10;            % Gesamtdauer [s]
t = (0:1/fs:T)';   % Zeitvektor

% Chirp-Signal
f0 = 1;            % Startfrequenz [Hz]
f1 = 100;          % Endfrequenz [Hz]
chirp_signal = chirp(t, f0, T, f1);

% Im Workspace speichern
assignin('base', 't_chirp', t);
assignin('base', 'chirp_signal', chirp_signal);

%% === Simulink-Modell starten ===
model_name = 'Projekt_Test';
simOut = sim(model_name, 'ReturnWorkspaceOutputs', 'on');

%% === Daten extrahieren ===
x = simOut.get('x');     % x(t) [Nx1] – Position
dx = simOut.get('dx');   % dx(t) [Nx2] – Spalte 1: Geschwindigkeit, Spalte 2: Beschleunigung
ddx = simOut.get('ddx');
u = simOut.get('i');     % Eingang u(t)
t = simOut.get('t');     % Zeitvektor

grad_dx = gradient(x,Ts);
grad_ddx = gradient(grad_dx,Ts);

%tvregdiff bessere funktion zum ableiten

tv_dx = TVRegDiff(x,1,1,[],[],1e-6,1,[],[]);
tv_dx(1,:)=[];
tv_dx = tv_dx/Ts;

tv_ddx = TVRegDiff(tv_dx,1,1,[],[],1e-6,1,[],[]);
tv_ddx(1,:)=[];
tv_ddx = tv_ddx/Ts;
size(t)
size(tv_dx)
size(tv_ddx)

% Extrahiere Signale
x1 = x(:);        % Position
v = tv_dx(:,1);      % Geschwindigkeit
a = tv_ddx(:,1);      % Beschleunigung (Zielgröße)

%% === Least Squares Ansatz ===
% a = -(C/m)*x1 - (D/m)*v + (K/m)*u
% => y = Phi * theta, mit:
%    y = a
%    Phi = [-x1, -v, u]
%    theta = [C/m; D/m; K/m]

Phi = [-x1, -v, u];
y = a;

% Least-Squares-Schätzung
theta_hat = Phi \ y;

% Parameter rekonstruieren
C_hat = theta_hat(1) * m;
D_hat = theta_hat(2) * m;
K_hat = theta_hat(3) * m;

set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontSize', 16);
% Plot: Position
figure;
plot(t, x(1:ds_factor:end), '--b', 'DisplayName', 'Original x(t)');
hold on;

xlabel('Zeit [s]');
ylabel('Position x(t)');
legend; grid on;
title('Positionsverlauf des Linearaktors');

% Plot: Geschwindigkeit
figure;
plot(t, dx, '--b', 'DisplayName', 'Original v(t) (Simulink)');
hold on;
plot(t, grad_dx, '-.r', 'DisplayName', 'v(t) aus Gradient');
hold on;
plot(t, tv_dx, '-.g', 'DisplayName', 'v(t) aus TVRegDiff');
hold on;

xlabel('Zeit [s]'); 
ylabel('Geschwindigkeit v(t)');
legend; grid on;
title('Vergleich der Ableitungsverfahren');

% Plot: Beschleunigung
figure;
plot(t, ddx, '--b', 'DisplayName', 'Original a(t) (Simulink)');
hold on;
plot(t, grad_ddx, '-.r', 'DisplayName', 'a(t) aus Gradient');
hold on;
plot(t, tv_ddx, '-.g', 'DisplayName', 'a(t) aus TVRegDiff');
hold on;


xlabel('Zeit [s]'); 
ylabel('Beschleunigung a(t)');
legend; grid on;
title('Vergleich der Ableitungsverfahren');
