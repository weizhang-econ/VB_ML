    % This function computes the density of y under FSV marginal of f and A
function lden = deny_fsv(X,Y,L,h,Hyper)
[T,n] = size(Y);
npr = size(h,2);
r = npr-n;
bigX = SURform2(X,n);
y = reshape(Y',T*n,1);
k_alp = size(bigX,2);
Omega = sparse(1:T*r,1:T*r,reshape(exp(h(:,n+1:n+r))',T*r,1));
Sig = sparse(1:T*n,1:T*n,reshape(exp(h(:,1:n))',T*n,1));

Sy = kron(speye(T),L)*Omega*kron(speye(T),L') + Sig;
XiSy = bigX'/Sy;
Kalp = sparse(1:k_alp,1:k_alp,1./Hyper.Valp) + XiSy*bigX;
CKalp = chol(Kalp,'lower');
tmpc = CKalp\(Hyper.alp0./Hyper.Valp + XiSy*y);      
lden = -T*n/2*log(2*pi) -.5*ldet(Sy) -.5*sum(log(Hyper.Valp)) -sum(log(diag(CKalp)))...
    -.5*(y'*(Sy\y) +sum(Hyper.alp0.^2./Hyper.Valp) -sum(tmpc.^2));            
end


