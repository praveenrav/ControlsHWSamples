%% Problem 3

close all
clear
clc

% Defining State Matrices:
G = [0.5 1; 0.5 0.7];
H = [0.2; 0.1];
C = [1 0];

yd = 1; % Desired output

%%% Part c:
p_des = [(0.4 + 0.3i), (0.4 - 0.3i)]; % Desired poles
K = place(G, H, p_des); % Calculating state feedback controller gains
G_cl = (G - (H * K));
gain = C * inv(eye(2) - G_cl) * H; % Calculating the reference input gain
r = 1./gain; % Calculating the reference input r for tracking control

%%% Part d:
G_mod = [G [0; 0]; C 1]; % Modified G matrix
H_mod = [H; 0]; % Modified H matrix
p_des = [(0.4 + 0.3i), (0.4 - 0.3i), 0.5]; % Desired poles
K_int = place(G_mod, H_mod, p_des); % Calculating state feedback controller gains

%%% Part e:
d = 1; % Disturbance
F = [0.1; -0.1]; % Disturbance matrix

% Recalculating the state matrices:
G_st = G_mod - (H_mod * K_int);
H_st = [(((d./yd) * F) + (H./gain)); -1];
C_st = [1 0 0];
D_st = 0;

% Linear Simulation:
sys_int = ss(G_st, H_st, C_st, D_st, 1); % Creating a discrete state-space model with T = 1 second
t_mod = [0:20]; % Time steps
u_mod = ones(1, length(t_mod)); % Desired constant output y_d (y_d = 1 in this case)

[output, tOut, states] = lsim(sys_int, u_mod, t_mod); % Performing discretized linear simulation
states = states'; 
x1_k = states(1, :); % x1(k)
x2_k = states(2, :); % x2(k)
u_k = (-K_int * states) + r; % Calculating u(k)

% Plotting Results:
figure
subplot(3, 1, 1)
stairs(tOut, x1_k, 'b', 'LineWidth', 0.75) % Plotting x1(k) with ZOH
xlabel('k')
ylabel('x1(k)')
set(gca, 'FontSize', 15)
grid on

subplot(3, 1, 2)
stairs(tOut, x2_k, 'b', 'LineWidth', 0.75) % Plotting x2(k) with ZOH
xlabel('k')
ylabel('x2(k)')
set(gca, 'FontSize', 15)
grid on

subplot(3, 1, 3)
stairs(tOut, u_k, 'r', 'LineWidth', 0.75) % Plotting u(k) with ZOH
xlabel('k')
ylabel('u(k)')
set(gca, 'FontSize', 15)
grid on
