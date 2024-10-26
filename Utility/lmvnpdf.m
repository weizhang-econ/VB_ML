function lden = lmvnpdf(x, mu, Sig)
    n = length(mu);
    CSig = chol(Sig,'lower');
    e = CSig\(x-mu);
    lden = - n/2*log(2*pi) - sum(log(diag(CSig))) - .5*(e'*e);
end