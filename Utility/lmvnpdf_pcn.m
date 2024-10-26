function lden = lmvnpdf_pcn(x,mu,K)
n = length(mu);
CK = chol(K,'lower');
e = CK'*(x-mu); 
lden = - n/2*log(2*pi) + sum(log(diag(CK))) - .5*(e'*e);
end