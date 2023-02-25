%% Question 1

close all
clear
clc

syms w g1 g2 g3 s;

% Setting up the system matrices:
A = [0 1 0; 0 0 1; 0 (-w^2) 0];
C = [1 0 0];

% Setting up observer:
L = [g1; g2; g3];
A_LC = A - (L*C);

% Calculating the characteristic equation of A_LC:
A_LC_eig = A_LC - (eye(3) .* s);
det(A_LC_eig);

%% Question 2

close all
clear
clc

syms k1 k2 k3 lam w s; % Necessary symbols for calculations

K = 0.5 .* (w.^2); % Setting the expression for K

% Setting up system matrices:
A = [0 0 -K; 0 0 K; 1 -1 0];
B = [0; 1; 0];
C = [1 0 0];

% Calculating gains for full-state feedback:
k = [k1 k2 k3];

ABk = A - (B*k);
ABk_lam = ABk - (eye(3) * lam);
det(ABk_lam); % Characteristic Equation of A - (B * k)

k = [w, 3*w, (-3 .* (w.^2))]; % Gains for full-state feedback calculated
ABk = A - (B*k);

G_s = C * inv((s*eye(3)) - ABk) * B; % Finding the transfer function Y(s)/R(s)

% Designing the full-state observer:
L = [k1; k2; k3];
A_LC = A-(L*C);
A_LC_eig = (A_LC - (eye(3) * lam));
det(A_LC_eig); % Characteristic Equation of A - (L * C)

% Designing the reduced-order observer:
A11 = 0;
A12 = [0 -K];
A21 = [0; 1];
A22 = [0 K; -1 0];
B1 = 0;
B2 = [1; 0];

L_r = [399; (-40/w)]; % Calculated gains for the reduced-order observer

% System matrices for the reduced-order observer:
A_r = A22 - (L_r * A12);
B_u = B2 - (L_r * B1);
B_y = (A_r * L_r) + A21 - (L_r * A11);
K1 = k(1);
K2 = k(2:3);

% Calculating F1(s) and F2(s):
F1_s = simplify(inv(1 + (K2 * inv(s * eye(2) - A_r) * B_u)));
F2_s = simplify((K2*inv(s * eye(2) - A_r)*B_y) + K1 + (K2 * L_r));

% Calculating Gc1_s and Gc2_s:
Gc1_s = C * inv(s * eye(3) - A + (B * k)) * B;

P = C * inv(s * eye(3) - A) * B;
Gc2_s = simplify(inv(1 + (P * F1_s * F2_s)) * P * F1_s);

%% Question 3

close all
clear
clc

syms s;

% Setting up system matrices:
A = [0 1 0 0; 30 0 0 20; 0 0 0 1; -2 0 0 -10];
B = [0; -25; 0; 10];
C = [40 0 0 0; 0 0 10 0; 0 0 0 -5];
D = zeros(3, 1);

%%% Combined full-state feedback and full-order observer:

% Designing the controller:
Pd_U = [-1; -2; -3; -4];
K = place(A, B, Pd_U);

% Designing the observer:
Pd_O = 10 .* Pd_U;
L = place(A', C', Pd_O)';

% Setting up the system matrices of the combined full-state feedback and
% full-order observer system:
A_f = [(A - B*K), (-B*K); zeros(4, 4), (A - L*C)];
B_f = [B; zeros(4, 1)];
C_f = [C, C];
D_f = zeros(3, 1);

csys1 = ss(A_f, B_f, C_f, D_f); % Creating a continuous system using the above system matrices

% Performing linear simulation:
t = [0:0.1:10];
u = ones(length(t),1);
x0 = [1 1 1 1 -1 -1 -1 -1];
[y1, t1, x1] = lsim(csys1,u,t,x0);
S1 = lsiminfo(y1,t1);

% Plotting results:
figure
subplot(3, 1, 1)
plot(t1, y1(:, 1), 'b')
xlabel('Time (s)')
ylabel('Y_1(t)')

subplot(3, 1, 2)
plot(t1, y1(:, 2), 'b')
xlabel('Time (s)')
ylabel('Y_2(t)')

subplot(3, 1, 3)
plot(t1, y1(:, 3), 'b')
xlabel('Time (s)')
ylabel('Y_3(t)')

sgtitle('Step Response w/ Full-Order Observer')

figure
subplot(4, 1, 1)
plot(t1, x1(:, 5), 'r')
xlabel('Time (s)')
ylabel('X_1(t) Error')

subplot(4, 1, 2)
plot(t1, x1(:, 6), 'r')
xlabel('Time (s)')
ylabel('X_2(t) Error')

subplot(4, 1, 3)
plot(t1, x1(:, 7), 'r')
xlabel('Time (s)')
ylabel('X_3(t) Error')

subplot(4, 1, 4)
plot(t1, x1(:, 8), 'r')
xlabel('Time (s)')
ylabel('X_4(t) Error')

sgtitle('Observer State Error for Full-Order Observer')

%%% Combined state feedback and reduced-order observer:

% Performing transformation on system matrices:
R = [0 1 0 0];
T = [C; R];
A_t = T * A * inv(T);
B_t = T * B;
C_t = C * inv(T);

K = place(A_t, B_t, Pd_U); % Recalculating the gains of the state-feedback controller since the state matrices are transformed

% Dividing into measurable and unmeasurable states:
A11 = A_t(1:3, 1:3);
A12 = A_t(1:3, 4);
A21 = A_t(4, 1:3);
A22 = A_t(4, 4);
B1 = B_t(1:3);
B2 = B_t(4);

% Calculating the reduced-observer gain:
L_r = place(A22', A12', [-10]);
L_r = L_r';

% Splitting up the state-feedback controller gains:
K1 = K(1:3);
K2 = K(4);

% Calculating the matrices of the combined state feedback and reduced-order
% observer system:
A_r = A22 - (L_r * A12);
A_red = [(A_t - B_t*K) (-B_t*K2); zeros(1, 4) A_r];
B_red = [B_t; 0];
C_red = [C_t, zeros(3, 1)];
D_red = zeros(3, 1);

csys2 = ss(A_red, B_red, C_red, D_red); % Creating a continuous system using the above system matrices

% Performing Linear Simulation:
t2 = [0:0.1:10];
u2 = ones(length(t2),1);
x0_2 = [1 1 1 1 -1];
[y2, t2, x2] = lsim(csys2, u2, t2, x0_2);
S2 = lsiminfo(y2, t2);

% Plotting results:
figure
subplot(4, 1, 1)
plot(t2, y2(:, 1), 'b')
xlabel('Time (s)')
ylabel('Y_1(t)')
title('Step Response w/ Reduced-Order Observer')

subplot(4, 1, 2)
plot(t2, y2(:, 2), 'b')
xlabel('Time (s)')
ylabel('Y_2(t)')

subplot(4, 1, 3)
plot(t2, y2(:, 3), 'b')
xlabel('Time (s)')
ylabel('Y_3(t)')

subplot(4, 1, 4)
plot(t2, x2(:, 5), 'r')
xlabel('Time (s)')
ylabel('X_1(t) Error')
title('Observer State Error for Reduced-Order Observer')

% Determining which output the system can be observed using rank test:
A = [0 1 0 0; 30 0 0 20; 0 0 0 1; -2 0 0 -10]; % State matrix

% First output:
C1 = [40 0 0 0; 0 0 0 0; 0 0 0 0];
Q1 = [C1; C1*A; C1*A*A; C1*A*A*A];
rQ1 = rank(Q1);

% Second output:
C2 = [0 0 0 0; 0 0 10 0; 0 0 0 0];
Q2 = [C2; C2*A; C2*A*A; C2*A*A*A];
rQ2 = rank(Q2);

% Third output:
C3 = [0 0 0 0; 0 0 0 0; 0 0 0 -5];
Q3 = [C3; C3*A; C3*A*A; C3*A*A*A];
rQ3 = rank(Q3);