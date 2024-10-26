%% This script returns the log density of the inverse-Wishart
%  distribuiton

function lden = liwpdf(Sig,nu0,S0,n)
cSig = -nu0*n/2*log(2)  -n*(n-1)/4*log(pi)- sum(gammaln((nu0+1-(1:n))/2))...
    + .5*nu0*ldet(S0);

lden = cSig - .5*(n+nu0+1)*ldet(Sig) - .5*trace(Sig\S0);
end