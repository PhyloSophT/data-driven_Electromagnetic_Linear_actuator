clear
close all
clc
%% Parameter für CHIRP-Signal
fs = 1000;         % Abtastfrequenz [Hz]
T = 10;            % Gesamtdauer [s]
t = (0:1/fs:T)';   % Zeitvektor

% Verstärkte Amplitude oder variable Verstärkung
f0 = 1; f1 = 10;
amplitude = 1;  % Erhöhte Anregung
chirp_signal = amplitude * chirp(t, f0, T, f1);

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
x_punkt =  simOut.get('dx_sim');
x_punktt =  simOut.get('ddx_sim');
x_punkt = squeeze(x_punkt);
x_punktt = squeeze(x_punktt);
I = simOut.get('i');     % Eingang u(t)
t = simOut.get('t');     % Zeitvektor
%% Verschiedene Ableitungsverfahren

dx = TVRegDiff(x,1,1,[],[],1e-6,1,[],[]);
dx(1,:)=[];
dx = dx/Ts;

ddx = TVRegDiff(dx,1,1,[],[],1e-6,1,[],[]);
ddx(1,:)=[];
ddx = ddx/Ts;

u = dx;
y = F_R*atan(M*u)+D*u;

% Downsampling-Faktor
ds_factor = 1;

% Downsampling der Signale
u = u(1:ds_factor:end);
y = y(1:ds_factor:end);
t = t(1:ds_factor:end);
dx = dx(1:ds_factor:end);
ddx = ddx(1:ds_factor:end);
I = I(1:ds_factor:end);

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

v = dx(:,1);      % Geschwindigkeit
a = ddx(:,1);

Phi = [I , -Psi_RBF]; %%%%falls bei der identifikation mehrere therme von x, xpunkt oder xpunktpunkt müssen diese zusammengefasst werden und können nicht mehr einzeln betrachtet werden
Y = a;
% Least-Squares-Schätzung
theta_hat = Phi \ Y;

% Parameter rekonstruieren
theta_hat = theta_hat*m;

K_hat = theta_hat(1,1);
theta_FR = theta_hat(2:end,1);
FR_RBF = Psi_RBF*theta_FR;

%[FR_RBF1,~,K_hat1,D_hat1]= Function_RBF(dx,ddx,F_R,M,m);

figure(1)
plot(u,FR_RBF) %Implementation von Dämpfung ins RBF-Netz
hold on
%plot(u,FR_RBF1+D_hat1*u) %Ohne Implementation von Dämpfung ins RBF-Netz
hold on
plot(u,y) %Zielverlauf
hold on

%% 7. Interpolationsfunktion der identifizierten Reibung
FR_fun = @(v_in) interp1(u, FR_RBF, v_in, 'spline', 'extrap');

%% 8. Systemrekonstruktion mit Runge-Kutta
n = length(t);
% Systemrekonstruktion mit Runge-Kutta (inkl. Geschwindigkeit und Beschleunigung)
x_rec = zeros(n,1);
v_rec = zeros(n,1);
a_rec = zeros(n,1); % Beschleunigung
x_rec(1) = 0;
v_rec(1) = 0;

for i = 1:n-1
    u_i = K_hat * I(i);
    
    % Runge-Kutta-Schritte (wie vorher)
    f1 = v_rec(i);
    g1 = (1/m) * (u_i - FR_fun(v_rec(i)));

    f2 = v_rec(i) + 0.5 * Ts * g1;
    g2 = (1/m) * (u_i - FR_fun(v_rec(i) + 0.5 * Ts * g1));

    f3 = v_rec(i) + 0.5 * Ts * g2;
    g3 = (1/m) * (u_i - FR_fun(v_rec(i) + 0.5 * Ts * g2));

    f4 = v_rec(i) + Ts * g3;
    g4 = (1/m) * (u_i - FR_fun(v_rec(i) + Ts * g3));

    x_rec(i+1) = x_rec(i) + (Ts/6)*(f1 + 2*f2 + 2*f3 + f4);
    v_rec(i+1) = v_rec(i) + (Ts/6)*(g1 + 2*g2 + 2*g3 + g4);
    
    % Beschleunigung speichern (aus letzter g4-Schätzung oder Mittelwert)
    a_rec(i) = (1/m) * (u_i - FR_fun(v_rec(i)));
end
a_rec(end) = a_rec(end-1);  % Letzter Wert auffüllen

%% 9. Plot
% Plot: Position
figure;
plot(t, x(1:ds_factor:end), '--b', 'DisplayName', 'Original x(t)');
hold on;
plot(t, x_rec, 'r', 'DisplayName', 'Rekonstruiert x(t)');
xlabel('Zeit [s]'); ylabel('Position x(t)');
legend; grid on;
title('Rekonstruktion der Position');

% Plot: Geschwindigkeit
figure;
plot(t, x_punkt, '--b', 'DisplayName', 'Original v(t) (Simulink)');
hold on;
plot(t, dx, '-.g', 'DisplayName', 'v(t) aus TVRegDiff');
plot(t, v_rec, 'r', 'DisplayName', 'Rekonstruiert v(t)');
xlabel('Zeit [s]'); ylabel('Geschwindigkeit v(t)');
legend; grid on;
title('Rekonstruktion der Geschwindigkeit');

% Plot: Beschleunigung
figure;
plot(t, x_punktt, '--b', 'DisplayName', 'Original a(t) (Simulink)');
hold on;
plot(t, ddx, '-.g', 'DisplayName', 'a(t) aus TVRegDiff');
plot(t, a_rec, 'r', 'DisplayName', 'Rekonstruiert a(t)');
xlabel('Zeit [s]'); ylabel('Beschleunigung a(t)');
legend; grid on;
title('Rekonstruktion der Beschleunigung');


%figure(4);
%plot(t,ddx);
%hold on
%plot(t,dx);
%hold off