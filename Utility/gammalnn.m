function loggamn = gammalnn(x,n)
loggamn = n*(n-1)/4*log(pi) + sum(gammaln((x+1-(1:n))/2));
end