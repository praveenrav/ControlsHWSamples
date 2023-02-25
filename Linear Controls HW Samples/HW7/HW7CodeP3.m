close all
clear
clc

% Importing State Matrices:
A = [-2 1; -5 0];
B = [1; 2];
C = [1 0];
D = 0;

% Performance Index Matrices
Q = zeros(2, 2);
P = [1 0; 0 0];

% Call ode45 to solve for Pb
T=[0:0.01:5]'; 
XP0=P(:);

% Can change r-value as necessary:
R = 10;

[T, XPb] = ode45('ricfun', T, XP0, [], A, B, Q, R);
XP = flipud(XPb);

for k = 1:length(T)
    Pk = reshape(XP(k,:)',2,2);
    K(k,:) = inv(R)*B'*Pk;
end

% Call ode45 to simulate the system response
[T,X]=ode45('sysfun', T, [1;0], [], A, B, K, T);

inp = -K.*X; % Solving for the input trajectory

% Plotting results:
figure
plot(T, inp(:, 1), 'b')
xlabel('Time (s)')
ylabel('Input Trajectory')
set(gca, 'FontSize', 18)

figure
plot(T, X(:, 1), 'r')
xlabel('Time (s)')
ylabel('Output Trajectory')
set(gca, 'FontSize', 18)

% Helper Functions:
function dx = ricfun(t,x,flag,A,B,Q,R)
%Ricatti Eq. Function to be called by ode45
Pb=reshape(x,2,2);
dPb=Pb*A+A'*Pb+Q-Pb*B*inv(R)*B'*Pb;
dx=dPb(:);
end

function dx=sysfun(t,x,flag,A,B,K,T)
%System Function to be called by ode45
k=max(find(t>=T)); %determine the index k
Kk=K(k,:);
dx=(A-B*Kk)*x;
end
