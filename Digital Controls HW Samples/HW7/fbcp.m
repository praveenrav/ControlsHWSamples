function [numc, denc, E, ND] = fbcp(nump, denp, dend)

n = length(denp) - 1;
E = zeros(2*n,2*n);

for jj = 1:n
    E(jj:jj+n,jj) = denp';
    E(jj:jj+n,jj+n) = nump';
end

if(length(dend)~=(2*n))
    y = [dend'; zeros(n-1,1)];
else
    y = dend';
end
Ei = inv(E);
ND = Ei(:, 1:n+1);
x = Ei*y;
denc = x(1:n).';
numc = x((n+1):(2*n)).';
end
