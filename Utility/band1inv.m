% This function computes the inverse of band matrix

function [dz0, dz1] = band1inv(K)
n = size(K,1);
C = chol(K,'lower');
is = 1./diag(C);
L = C*sparse(1:n,1:n,is);
u = full(diag(L,-1));
A = speye(n) + sparse(1:n-1,2:n,-u.^2,n,n);
dz0 = A\(is.^2);    % main diagonal 
dz1 = -u.*dz0(2:n); % 1st diagonal
end