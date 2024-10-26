%% E(st^(-2)) = E((1+(s1-1)rho^(t-tstar-2))^(-2))
function meanstsq = E_stsq(s2grid,rhogrid,ptheta,stfunc,T,tstar)
tt = tstar+3:T;
tmp = zeros(length(s2grid),length(rhogrid),T-(tstar+3)+1);
for jj = 1:length(s2grid)
    s2j = s2grid(jj);
    for kk = 1:length(rhogrid)
        rhok = rhogrid(kk);
        tmp(jj,kk,:) = stfunc(s2j,rhok,tt).^(-2);
    end
end
tmpptheta = repmat(ptheta,1,1,T-(tstar+3)+1);
meanstsq = squeeze(sum(sum(tmpptheta.*tmp),2));
end