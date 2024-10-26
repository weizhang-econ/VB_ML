%% theta = [s2, rho], row of thetagrid is grids for s2, and column of thetagrid is grids for rho
function logptheta = approx_logptheta_bytime(Y,Z,alp,iSigma,tstar,T,n,s2x,rhox,iKalp,stfunc,alpha,beta,test)


tt = tstar+3:T;
sttmp = stfunc(s2x,rhox,tt);
EM1 = sum(log(sttmp));
tmp = 0;
for tt = tstar+3:T
    xt = Z(tt,:);
    Xt = kron(speye(n),xt);
    yt = Y(tt,:)';
    stsq = sttmp(tt-tstar-2)^2;
    tmp = tmp + ((yt-Xt*alp)'*iSigma*(yt-Xt*alp)+trace(iSigma*Xt*iKalp*Xt'))/stsq;   
end
tt=tstar+2;
xt = Z(tt,:);
Xt = kron(speye(n),xt);
yt = Y(tt,:)';
logptheta =  -(.5*n+1)*log(s2x^2)-n*EM1-((yt-Xt*alp)'*iSigma*(yt-Xt*alp)+trace(iSigma*Xt*(iKalp*Xt')))/(2*s2x^2)-.5*tmp+(alpha-1)*log(rhox) + (beta-1)*log(1-rhox)-betaln(alpha,beta);

% %% if using iKalp
% %CiKalp = chol(iKalp,'lower'); somehow when using chol(iKalp, 'lower'),
% %sometimes Matlab returns the error: matrix must be positive definite.
% %However, if using chol(iKalp)', it will be fine.
% iKalp = (iKalp+iKalp')/2;
% CiKalp = chol(iKalp,'lower');
% tmp1 = Xttmp*CiKalp;
% tmp = resid'*tmpkr*resid + trace(tmpkr*(tmp1*tmp1'));
% 
% %% if using Kalp
% % CKalp = chol(Kalp,'lower');
% % tmp1 = CKalp\Xttmp';
% %tmp = resid'*tmpkr*resid + trace(tmpkr*(tmp1'*tmp1));
% 
% logptheta =  -(.5*n+1)*log(s2x^2)-n*EM1-.5*tmp+(alpha-1)*log(rhox) + (beta-1)*log(1-rhox)-betaln(alpha,beta);
end