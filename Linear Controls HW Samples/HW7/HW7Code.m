%% Question 1

close all
clear
clc

syms w s

K = (w.^2)./2;

% Importing State Matrices:
A = [0 0 -K; 0 0 K; 1 -1 0];
B = [0; 1; 0];
C = [1 0 0];

G_s = C * inv(s * eye(3) - A) * B; % Obtains the plant's open-loop transfer function

% Solving for D(s)
lam_s = s^2 + ((20*w*s) + (200 * (w.^2))); % Desired characteristic equation for the observer
phi_s = s^3 + ((3 * w) * (s^2)) + ((4 * (w.^2))*s) + (2 * (w.^3)); % Desired characteristic equation for the state feedback controller
D_s = expand(lam_s * phi_s);

% Solving the Matrix Equation A_mat * x = B_mat for this system of equations:
A_mat = [1 0 0 0 0 0; 0 1 0 0 0 0; (w.^2) 0 1 0 0 0; 0 (w.^2) 0 (0.5*w.^2) 0 0; 0 0 (w.^2) 0 (0.5*w.^2) 0; 0 0 0 0 0 (0.5*w.^2)];
B_mat = [1; (23 * w); (264 * (w.^2)); (682 * (w.^3)); (840 * (w.^4)); (400 * (w.^5))];
x_mat = inv(A_mat) * B_mat; % Solving for the alpha and beta parameters

% Extracting the alpha and beta coefficients:
alph = x_mat(1:3);
beta = x_mat(4:6);

%% Question 2

close all
clear
clc

% Importing state matrices and performance indices:
B = 1;
Q = 1;
R = [1:10:10000]; % Vector of possible R values
syms P;

% A = 0:

A = 0; % State matrix for this case
K1 = []; % Vector of K gains
p = []; % Vector of P matrices (scalars, in this case)

% Solving the P matrix and gains K for each value of r:
for i = 1:length(R)
    % Solving for the P matrix using ARE:
    ARE = (Q - ((P.^2)/R(i))) == 0;
    ARE_soln = double(solve(ARE, P));
    p_cur = ARE_soln(ARE_soln > 0);
    p = [p; p_cur];
    
    % Solving for the gains:
    K1_cur = inv(R(i)) * B' * p_cur;
    K1 = [K1; K1_cur];
end

% Plotting results:
figure
plot(R, p, 'b')
hold on
plot(R, K1, 'r')
xlabel('r')
ylabel('Resulting P Matrix and K Gains')
title('A = 0')
legend('P', 'K')
set(gca, 'FontSize', 18)

% A > 0:

A = 5; % State matrix for this case
K2 = []; % Vector of K gains
p = []; % Vector of P matrices (scalars, in this case)

% Solving the P matrix and gains K for each value of r:
for i = 1:length(R)
    % Solving for the P matrix using ARE:
    ARE = ((2 * A * P) + Q - ((P.^2)/R(i))) == 0;
    ARE_soln = double(solve(ARE, P));
    p_cur = ARE_soln(ARE_soln > 0);
    p = [p; p_cur];
    
    % Solving for the gains:
    K2_cur = inv(R(i)) * B * p_cur;
    K2 = [K2; K2_cur];
end

% Plotting results:
figure
plot(R, p, 'b')
hold on
plot(R, K2, 'r')
xlabel('r')
ylabel('Resulting P Matrix and K Gains')
title('A > 0')
legend('P', 'K')
set(gca, 'FontSize', 18)

% A < 0:

A = -5; % State matrix for this case
K3 = []; % Vector of K gains
p = []; % Vector of P matrices (scalars, in this case)

% Solving the P matrix and gains K for each value of r:
for i = 1:length(R)
    % Solving for the P matrix using ARE:
    ARE = ((2 * A * P) + Q - ((P.^2)/R(i))) == 0;
    ARE_soln = double(solve(ARE, P));
    p_cur = ARE_soln(ARE_soln > 0);
    p = [p; p_cur];
    
    % Solving for the gains:
    K3_cur = inv(R(i)) * B * p_cur;
    K3 = [K3; K3_cur];
end

% Plotting results:
figure
plot(R, p, 'b')
hold on
plot(R, K3, 'r')
xlabel('r')
ylabel('Resulting P Matrix and K Gains')
title('A < 0')
legend('P', 'K')
set(gca, 'FontSize', 18)

%% Problem 4

close all
clear
clc

% Importing system transfer functions:
Gs = tf([0 1 2], [1 2 5]);
Gnegs = tf([0 -1 2], [1 -2 5]);

% Plotting Root-Locus of System using "1/r" as the gain:
rlocus(Gs*Gnegs)

%% Problem 5

close all
clear
clc

% Importing State Matrices:
A = [-10 0 -10 0; 0 -0.7 9 0; 0 -1 -0.7 0; 1 0 0 0];
B = [20 2.8; 0 -3.13; 0 0; 0 0];

% Importing Performance Index Matrices:
Q = [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1];
R = [1 0; 0 1];

% Calculating the state-feedback gains using LQR Control
[P, K] = icare(A, B, Q, R);
