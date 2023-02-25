close all
clear all
clc

M = 5; % Mass in kg
M_var = 7.5; % Variable mass in kg
Kc = 500; % Spring constant in N/m
Kb_c = 10; % Motor back-emf constant
Kf = 10; % Motor force constant
R = 0.1; % Motor armature resistance in Ohms

%% Part 1:

% Continuous-time State Matrices:
A = [0 0 1; ((R*Kc)/(Kb_c*Kf)) ((-R*Kc)/(Kb_c*Kf)) 0; (-Kc/M) (Kc/M) 0];
B = [0; (1/Kb_c); 0];
%C = [1 0 0];
C = [1 0 0; 0 1 0];
%D = 0;
D = [0; 0];
sysc = ss(A, B, C, D);

% Discretizing State-Space Equation:
T = 0.01;
sysd = c2d(sysc, T);
G = sysd.A;
H = sysd.B;

% Mass-Varied System:
M_sys = M_var; % Varied mass of the system (can be changed to see how the controller behaves)
A_sys = [0 0 1; ((R*Kc)/(Kb_c*Kf)) ((-R*Kc)/(Kb_c*Kf)) 0; (-Kc/M_sys) (Kc/M_sys) 0];
B_sys = [0; (1/Kb_c); 0];
%C_sys = [1 0 0];
C_sys = [1 0 0; 0 1 0];
%D_sys = 0;
D_sys = [0; 0];
P_s = ss(A_sys, B_sys, C_sys, D_sys);
P_z = c2d(P_s, T);

%% Part 2 - State-Feedback Controller with Reduced-Order Deadbeat Observer:

% Partitioning the matrices for the reduced-order observer:
n_o = 2;
Gaa = G((1:n_o), 1:n_o);
Gab = G((1:n_o), ((n_o+1):end));
Gba = G(((n_o+1):end), (1:n_o));
Gbb = G(((n_o+1):end), ((n_o+1):end));
Ha = H(1:n_o);
Hb = H((n_o+1):end);

%%% State-Feedback Controller:
pd_k = [(-25 + 25i), (-25 - 25i), -25];
pd_k = exp(pd_k * T);

K = place(G, H, pd_k); % Calculating state-feedback gains

% Calculating reference input gain:
G_cl = G - (H*K); 
gain = C * inv(eye(3) - G_cl) * H;
gain = gain(1);

Ka = K(1:n_o);
Kb = K((n_o+1):end);

%%% Reduced-Order Deadbeat Observer:

% Forming the deadbeat poles:
if(n_o == 1)    
    pd_o = [0, 0];
elseif(n_o == 2)
    pd_o = 0;
end

L = place(Gbb', Gab', pd_o);
L = L';

% Calculating more system matrices:
Gr = Gbb - (L*Gab);
Hu = Hb - (L*Ha);
Hy = (Gr*L) + Gba - (L*Gaa);

%%% Feedback System:
F2 = ss(Gr, Hy, Kb, (Ka+(Kb*L)), T);
%F2 = tf(F2);
%F2 = F2(1);

%%% Feedforward System:
ff1 = ss(Gr, Hu, Kb, 0, T);
F1 = feedback(1, ff1);
%F1 = tf(F1);

%% Part 3:

% Controlled System:
G_cltf = feedback((F1*sysd), F2);
G_cltf_m = feedback((F1*P_z), F2);

t = [0:T:0.6]; % Time vector for linear simulation

% Calculating the reference input:
yd = 1;
r = yd./gain; % Reference input r
u = ones(1, length(t)).*r; % Reference input vector

% Actual system:
[y, tOut, x_states] = lsim(G_cltf, u, t);
y = y(:, 1);

U_R_z = F1/(1 + (F2*F1*sysd));
u_input = ones(1, length(t)) .* r;
[u_ctrl, ~,  ~] = lsim(U_R_z, u_input, t);

% Varied-mass system:
[y_m, ~, x_states_m] = lsim(G_cltf_m, u, t);
y_m = y_m(:, 1);

U_R_z_m = F1/(1 + (F2*F1*P_z));
[u_ctrl_m, ~,  ~] = lsim(U_R_z_m, u_input, t);

% Plotting Results:
lgnd_lbl = sprintf('M = %0.1f kg', M_var);

figure
subplot(2, 1, 1)
stairs(tOut, y, 'b', 'LineWidth', 1)
hold on
stairs(tOut, y_m, 'r', 'LineWidth', 1)
xlabel('Time [s]')
ylabel('Controller Output [m]')
title('Combined State Feedback/ROO Control Response vs. Time')
legend('M = 5 kg', lgnd_lbl)
set(gca, 'FontSize', 15)
grid on

subplot(2, 1, 2)
stairs(tOut, u_ctrl, 'b', 'LineWidth', 1)
hold on
stairs(tOut, u_ctrl_m, 'r', 'LineWidth', 1)
title('Combined State Feedback/ROO Control Input vs. Time')
xlabel('Time [s]')
ylabel('u')
legend('M = 5 kg', lgnd_lbl)
set(gca, 'FontSize', 15)
grid on

%% Re-defining Plant Systems:

% Actual Plant:
M_sys = 5; % Mass of the system (can be changed to see how the controller behaves)

A_sys = [0 0 1; ((R*Kc)/(Kb_c*Kf)) ((-R*Kc)/(Kb_c*Kf)) 0; (-Kc/M_sys) (Kc/M_sys) 0];
B_sys = [0; (1/Kb_c); 0];
C_sys = [1 0 0];
%C_sys = [1 0 0; 0 1 0];
D_sys = 0;
%D_sys = [0; 0];
P_s = ss(A_sys, B_sys, C_sys, D_sys);
P_s = tf(P_s);
P_z = c2d(P_s, T); % Discretized plant transfer function

% Mass-Varied Plant:
M_sys2 = M_var; % Mass of the system (can be changed to see how the controller behaves)

A_sys = [0 0 1; ((R*Kc)/(Kb_c*Kf)) ((-R*Kc)/(Kb_c*Kf)) 0; (-Kc/M_sys2) (Kc/M_sys2) 0];
B_sys = [0; (1/Kb_c); 0];
C_sys = [1 0 0];
%C_sys = [1 0 0; 0 1 0];
D_sys = 0;
%D_sys = [0; 0];
P_s = ss(A_sys, B_sys, C_sys, D_sys);
P_s = tf(P_s);
P_z_m = c2d(P_s, T); % Discretized plant transfer function

%% Part 4 - Polynomial Controller:

P_tf = tf(sysd); % Converting plant system to transfer function
P_tf = P_tf(1); % Taking only the first output (which is y)

[A_z, B_z] = tfdata(P_tf, 'v'); % Extracting numerator and denominator coefficients from plant transfer function

% Creating phi(z):
syms z
phi_z = expand(prod(z - pd_k));
phi_coef = double(coeffs(phi_z));
phi_coef = phi_coef(end:-1:1);

lamd_z = [1 0 0]; % Lambda(z) for a deadbeat observer

% Controller coefficients:
[numc, denc] = fbcp(A_z, B_z, phi_coef); % Designing the polynomial controller

% Creating the controller feedforward and feedback blocks:
F1_z = tf(lamd_z, denc, T); % Feedforward
F2_z = tf(numc, lamd_z, T); % Feedback

%%% Performing Linear Simulations:

% Actual system:
G_cltf_poly = feedback((F1_z*P_z), F2_z);
G_cltf_poly = minreal(G_cltf_poly);
[y_poly, tOut, ~] = lsim(G_cltf_poly, u, t);

U_R_z = F1_z/(1 + (F1_z*F2_z*P_z));
u_input = ones(1, length(t)) .* r;
[u_poly, ~,  ~] = lsim(U_R_z, u_input, t);

% Varied-mass system:
G_cltf_poly_m = feedback((F1_z*P_z_m), F2_z);
G_cltf_poly_m = minreal(G_cltf_poly_m);
[y_poly_m, ~, ~] = lsim(G_cltf_poly_m, u, t);

U_R_z_m = F1_z/(1 + (F1_z*F2_z*P_z_m));
[u_poly_m, ~,  ~] = lsim(U_R_z_m, u_input, t);

% Plotting Results:
figure
subplot(2, 1, 1)
stairs(tOut, y_poly, 'b', 'LineWidth', 1)
hold on
stairs(tOut, y_poly_m, 'r', 'LineWidth', 1)
xlabel('Time [s]')
ylabel('Controlled System Output [m]')
title('Polynomial Controlled-System Output vs. Time')
legend('M = 5 kg', lgnd_lbl)
set(gca, 'FontSize', 15)
grid on

subplot(2, 1, 2)
stairs(tOut, u_poly, 'b', 'LineWidth', 1)
hold on
stairs(tOut, u_poly_m, 'r', 'LineWidth', 1)
xlabel('Time [s]')
ylabel('u')
title('Polynomial Controlled-System Input vs. Time')
legend('M = 5 kg', lgnd_lbl)
set(gca, 'FontSize', 15)
grid on

%% Part 5 - Addition of ZPET Controller:

G_cltf_poly = tf(G_cltf_poly);
[B_z, D_z] = tfdata(G_cltf_poly, 'v'); % Extracting numerator and denominator coefficients from polynomial-controlled transfer function
D_z_tf = tf(D_z, [1], T); % Creating transfer function using denominator coefficients

[~, ~, gain2] = zpkdata(G_cltf_poly, 'v'); % Extracting the gain of the transfer function

G_zeros = roots(B_z); % Obtaining all zeros of the closed loop system
B_z_plus = G_zeros(abs(G_zeros) <= 1); % Obtaining the good zeros
B_z_minus = G_zeros(abs(G_zeros) > 1); % Obtaining the bad zeros

%%% Creating Q(z):

% Creating B+(z)'s z-polynomial:
syms z
B_z_plus = expand(prod(z - B_z_plus));
B_z_plus = double(coeffs(B_z_plus));
B_z_plus = B_z_plus(end:-1:1);
B_z_plus = tf(B_z_plus, [1], T); % B+(z)

% Creating B-(z)'s z-polynomial, B-(z^-1)'s z-polynomial, and B-(1):
B_z_minus = expand(prod(z - B_z_minus)); 
B_z1_minus = double(coeffs(B_z_minus)); 
B_z_minus = B_z1_minus(end:-1:1);
B_z_minus = tf(B_z_minus, [1], T); % B-(z)
B_z1_minus = tf(B_z1_minus, [1 0], T); % B-(z^-1)
B_minus1 = evalfr(B_z_minus, 1); % B-(-1)

Q_z = (B_z1_minus * D_z_tf)/((B_minus1.^2) * B_z_plus * gain2); % Formula for Q(z) from slides
Q_z = tf([1], [1 0 0 0 0], T) * Q_z; % Adding poles at z = 0 until the transfer function's numerator degree equals its denominator degree

% Formulating time and desired trajectory vectors:
t_sim1 = [0:T:0.5];
t_sim2 = [(0.5+T):T:1];
t_sim3 = [(1+T):T:1.5];
t_sim4 = [(1.5+T):T:2];

y1_t = 4.*(t_sim1.^2).*(3 - (4.*t_sim1));
y2_t = ones(1, length(t_sim2));
y3_t = y2_t - y1_t(1:(end-1));
y4_t = zeros(1, length(t_sim4));

t_tr1 = [t_sim1 t_sim2 t_sim3 t_sim4]; % Time vector for 1 cycle
y_tr_t1 = [y1_t, y2_t, y3_t, y4_t]; % Desired trajectory vector for 1 cycle

n_des = 5; % Number of desired cycles
for i = 0:(n_des-1)
    if(i == 0)
        t_tr = t_tr1;
        y_tr = y_tr_t1;
    else
        t_tr_cur = t_tr1 + (i.*2);
        t_tr = [t_tr, t_tr_cur(2:end)];
        y_tr = [y_tr, y_tr_t1(2:end)];
    end
end

% Performing linear simulations:
[r_t, ~, ~] = lsim(Q_z, y_tr, t_tr);
[y_sim_zept, ~, ~] = lsim(G_cltf_poly, r_t, t_tr);
[y_sim_zept_m, ~, ~] = lsim(G_cltf_poly_m, r_t, t_tr);

%% Part 6 - Addition of Stable Repetitive Controller 

% Learning parameters for the controller:
Kl = 0.8;
N = 200;

% Creating the (z^N - 1) term:
N_poly = [1, zeros(1, (N-1)), -1];
N_poly = tf([1], N_poly, T); % (z^N - 1)

C_z = Kl * Q_z * N_poly; % Stable repetitive controller

% Systems' OLTFs:
G_oltf_st = C_z * G_cltf_poly; % Actual system
G_oltf_st_m = C_z * G_cltf_poly_m; % Varied-mass system

% Systems' CLTFs:
G_cltf_st = feedback(G_oltf_st, 1); % Actual system
G_cltf_st_m = feedback(G_oltf_st_m, 1); % Varied-mass system

% Performing linear simulations:
[y_sim_st, ~, ~] = lsim(G_cltf_st, y_tr, t_tr);
[y_sim_st_m, ~, ~] = lsim(G_cltf_st_m, y_tr, t_tr);

%% Part 7 - Plotting Results from 5 and 6

% Plotting Desired Trajectory:
figure
stairs(t_tr, y_tr, 'm', 'LineWidth', 1);
xlabel('Time [s]')
ylabel('y_d')
title('Desired Trajectory')
set(gca, 'FontSize', 15)
grid on

% Plotting ZEPT Controller Results:
figure
stairs(t_tr, y_sim_zept, 'b', 'LineWidth', 1)
%stairs(t_tr, y_sim_zept_m, 'c', 'LineWidth', 1)
xlabel('Time [s]')
ylabel('Output')
title('ZPET Controller Output vs. Time') 
legend(lgnd_lbl)
set(gca, 'FontSize', 15)
grid on

% Plotting Stable Repetitive Tracking Results w/ Desired Trajectory:
figure
stairs(t_tr, y_tr, 'r', 'LineWidth', 1);
hold on
stairs(t_tr, y_sim_st, 'b', 'LineWidth', 1)
%stairs(t_tr, y_sim_st, 'c', 'LineWidth', 1)
xlabel('Time [s]')
ylabel('Output')
title('Stable Repetitive Controller Output vs. Time') 
legend('Desired Trajectory', lgnd_lbl)
set(gca, 'FontSize', 15)
grid on