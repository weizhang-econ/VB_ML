%% theta = [s2, rho], row of thetagrid is grids for s2, and column of thetagrid is grids for rho
function logptheta = approx_logptheta(resid,CiSigma,tstar,T,n,s2x,rhox,tmp1,stfunc,alpha,beta)
% ytstar3 = reshape(Y(tstar+3:end,:)',n*(T-(tstar+3)+1),1);
% Xstar3 = X((tstar+3-1)*n+1:end,:);
% resid = ytstar3 - Xstar3*alp;
% tt = tstar+3:T;
% sttmp = stfunc(s2x,rhox,tt);
% EM1 = sum(log(sttmp));
% M2 = 1./(sttmp.^2);
% tmpkr = kron(sparse(1:(T-(tstar+3)+1), 1:(T-(tstar+3)+1),M2),iSigma);
% 
% % tmp = 0;EM1 = 0;
% % for tt = tstar+3:T
% %     xt = Z(tt,:);
% %     Xt = kron(speye(n),xt);
% %     yt = Y(tt,:)';
% %     meanrhot = meanrhoT(tt); meanrhot2 = meanrhoT2(tt); varrhot = varrhoT(tt);
% %     sttmp = 1+(s2x-1)*meanrhot;
% %     M1 = log(sttmp)-.5*(s2x-1)^2*varrhot/sttmp^2;
% %     EM1 = EM1 + M1;
% % 
% %     tmpM2 = meanrhot2*s2x^2-(2*meanrhot2-2*meanrhot)*s2x+(meanrhot-1)^2;
% %     M2 = 1/tmpM2;
% %     if test
% %         tmp = tmp + M2*((yt-Xt*alp)'*iSigma*(yt-Xt*alp));
% %     else
% %         tmp = tmp + M2*((yt-Xt*alp)'*iSigma*(yt-Xt*alp)+trace(iSigma*Xt*(Kalp\Xt')));
% %     end
% % end
% tt=tstar+2;
% xt = Z(tt,:);
% Xt = kron(speye(n),xt);
% yt = Y(tt,:)';
% if test
%     tmp =  resid'*tmpkr*resid;
%     logptheta =  -(.5*n+1)*log(s2x^2)-((yt-Xt*alp)'*iSigma*(yt-Xt*alp))/(2*s2x^2)-n*EM1-.5*tmp+(alpha-1)*log(rhox) + (beta-1)*log(1-rhox);
% else
%     tmp = resid'*tmpkr*resid + trace(tmpkr*Xstar3*(Kalp\Xstar3'));
%     logptheta =  -(.5*n+1)*log(s2x^2)-((yt-Xt*alp)'*iSigma*(yt-Xt*alp)+trace(iSigma*Xt*(Kalp\Xt')))/(2*s2x^2)-n*EM1-.5*tmp+(alpha-1)*log(rhox) + (beta-1)*log(1-rhox)-betaln(alpha,beta);
% end





%tmpkr = kron(sparse(1:(T-(tstar+2)+1), 1:(T-(tstar+2)+1),1./(ststar2.^2)),iSigma);

%% if using iKalp
%CiKalp = chol(iKalp,'lower'); somehow when using chol(iKalp, 'lower'),
%sometimes Matlab returns the error: matrix must be positive definite.
%However, if using chol(iKalp)', it will be fine.

%CiSigma = chol(iSigma,'lower');

%tmp = resid'*tmpkr*resid + trace(tmpkr*tmp2);

%% if using Kalp
% CKalp = chol(Kalp,'lower');
% tmp1 = CKalp\Xttmp';
%tmp = resid'*tmpkr*resid + trace(tmpkr*(tmp1'*tmp1));

end