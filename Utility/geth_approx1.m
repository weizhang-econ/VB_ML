% approximate h using a linear state space model where the log chi^2 errors
% are modeled as N(-1.27,4.94)
function [h_hat, Kh] = geth_approx1(s2,h0,sigh2)
ystar = log(s2);
T = length(s2);
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HH = H'*H;
tmp = 4.94;
Kh = HH/sigh2 + speye(T)/tmp ;
h_hat = Kh\(h0/sigh2*HH*ones(T,1) + (ystar+1.27)/tmp );
end