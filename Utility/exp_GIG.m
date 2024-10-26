function [meanx, meaninvx, meanlogx] = exp_GIG(a,b,p)
% tmp = double(besselk(sym(p+1), sym(sqrt(a*b)))/besselk(sym(p), sym(sqrt(a*b))));
% %logmeanx = .5*log(b/a)+double(log(besselk(sym(p+1), sym(sqrt(a*b)))))-double(log(besselk(sym(p), sym(sqrt(a*b)))));
% meanx = sqrt(b/a)*tmp;
% meaninvx = sqrt(a/b)*tmp-2*p/b;

%% If using numerical differetiation
% h=1;
% pgrid = [p-h,p+h];
% tmpy = double(log(besselk(sym(pgrid), sym(sqrt(a*b)))));
% tmpderiv = diff(tmpy)/(2*h);
% meanlogx = log(sqrt(b/a))+tmpderiv;

% If using random draws from gigrnd and compute the mean
x = gigrnd(p,a,b,5000);
meanlogx = mean(log(x));
meanx = mean(x);
meaninvx = mean(1./x);

% %% If using numerical integration
% loglik = @(x) p/2*log(a/b)-log(2)-double(log(besselk(sym(p), sym(sqrt(a*b)))))+(p-1)*log(x)-0.5*(a*x+b./x);
% loglikx = loglik(xgrid);
% % i=1;
% % while exp(max(loglikx)) < 1e-10
% %     xgrid = [10*i:0.001:30*i];
% %     loglikx = loglik(xgrid);
% %     i=i+1;
% % end
% tmplikx = exp(loglikx - max(loglikx));
% meanlogx = trapz(xgrid,tmplikx.*log(xgrid))*exp(max(loglikx));
end