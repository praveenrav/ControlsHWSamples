function dx = ricfun(t,x,flag,A,B,Q,R)
%Ricatti Eq. Function to be called by ode45
Pb=reshape(x,2,2);
dPb=Pb*A+A'*Pb+Q-Pb*B*inv(R)*B'*Pb;
dx=dPb(:);
end

