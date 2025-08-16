function [FR_RBF, Psi_RBF,K_hat,D_hat] = Function_RBF(dx,ddx,F_R,M,m)
%% === Simulink-Modell starten ===
model_name = 'nichtlineares_ZRM_sim';
simOut = sim(model_name, 'ReturnWorkspaceOutputs', 'on');
I = simOut.get('i');     % Eingang u(t)
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

v = dx(:,1);      % Geschwindigkeit
a = ddx(:,1);

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

