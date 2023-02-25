function dx=sysfun(t,x,flag,A,B,K,T)
%System Function to be called by ode45
k=max(find(t>=T)); %determine the index k
Kk=K(k,:);
dx=(A-B*Kk)*x;
end
