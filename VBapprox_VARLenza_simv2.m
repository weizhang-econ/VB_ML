% This script estimates the BVAR in reduced-form with stationary autoregressive SV using VB
%clear all; %clc;
%function [A_VB,kappa_VB,rho_VB, s0sq_VB, s1sq_VB, s2_VB, mllbsv,lml_lenza,lmlstd_lenza,time_spent,time_spent_ml] = VBapprox_VARLenza_simv2(Y,Y0,ml_dummy,p,tstar)
[T,n] = size(Y);
y = reshape(Y',n*T,1);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p);
for i=1:p
    Z(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
Z = [ones(T,1) Z];
X = SURform2(Z,n);
k_alp = n^2*p + n;        % dimension of A
k = k_alp/n;

% priors
kappa1 = .2^2; kappa3 = 100; kappa4 = .2^2;kappa2 = (.2^2)^2;
[Hyper.alp0,Hyper.Valp,sig2_hat,~] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
Hyper.iValp = 1./Hyper.Valp;
Hyper.c0 = [1,1/.2^2; 1,1/.2^2];
[C_alp,idx_kappa1,idx_kappa2] = get_C(n,p,sig2_hat);


nus0 = .5*(n+1); nus1 = .5*(n+1);
alpha = 3; beta=1.5;
% initialize for storage
stfunc = @(s2x,rhox,tx) 1+(s2x-1).*(rhox.^(tx-tstar-2));
akappa = zeros(2,1);bkappa = zeros(2,1);ckappa = zeros(2,1);
meaninvkappa = (.2^2)^(-1)*ones(2,1);meanlogkappa = log(.2^2)*ones(2,1);meankappa = .2^2*ones(2,1);
%logk = @(nu,zz) 1/2*log(pi./(2*nu))-nu.*log(exp(1)*zz./(2*nu));
varalp = zeros(n,k);
sig2 = zeros(n,T);
A = zeros(k,n);
for ii=1:n
    iValpi = sparse(1:k,1:k,1./Hyper.Valp((ii-1)*k+1:ii*k));
    %     A(ii,:) = (XX + iValpi)\(X'*Y(:,ii));
    KAi = iValpi + Z'*Z;
    A(:,ii) = KAi\(sparse(1:k,1:k,Hyper.Valp((ii-1)*k+1:ii*k))\Hyper.alp0((ii-1)*k+1:ii*k) + Z'*Y(:,ii));
    varalp(ii,:) = diag(inv(KAi))';
    sig2(ii,:) = (Y(:,ii) - Z*A(:,ii)).^2;
end
meansig = sparse(1:n,1:n,(mean(sig2(:,1:tstar-1),2)));
%nusig = median(diag(Psisig))/meansig+1+n;
nusig = n+2;
Psisig = (nusig-n-1)*meansig; 

alp = reshape(A,k_alp,1);
invSig =  sparse(1:n,1:n,(mean(1./sig2(:,1:tstar-1),2)));
varA = reshape(varalp',k_alp,1);
s0 = 10;s1=30;s2=15;rho=0.6;
s0t_1 = s0; s1t_1 = s1; s2t_1 = s2; rhot_1 = rho;
meaninvssq = [ones(tstar-1,1);[s0^(-2),s1^(-2),s2^(-2)]';1./(1+(s2-1)*rho.^(tstar+3:T)).^2'];

B0 = eye(n);
rhogrid = (0.5:0.01:0.9)';s2grid = (50:0.5:100)';
logptheta = zeros(length(s2grid),length(rhogrid));
Like_alp = zeros(n,1);
% VB iteration starts here
TOL =  10; TOL_para = 1e-1;% tolerance level
mllb = -1e10; delta = 1e10; delta_para = 1e10; start_time = tic; iteration_idx = 0;
disp(['Starting VB Estimation...'])
ckappa(1) = Hyper.c0(1,1)-n*p/2; akappa(1) = 2*Hyper.c0(1,2); bkappa(1) =  sum(((alp(idx_kappa1).^2+varA(idx_kappa1))./C_alp(idx_kappa1)));
ckappa(2)= Hyper.c0(2,1)-(n-1)*n*p/2; akappa(2) = 2*Hyper.c0(2,2); bkappa(2) =  sum(((alp(idx_kappa2).^2+varA(idx_kappa2))./C_alp(idx_kappa2)));
Kalp_VB = zeros(k,k,n);
while delta > TOL|| delta_para > TOL_para

    for jj = 1:2
        [meankappa(jj), meaninvkappa(jj), meanlogkappa(jj)] = exp_GIG(akappa(jj),bkappa(jj),ckappa(jj));
    end
    [~,Hyper.iValp,logV_Minn] = update_iValp(p,meaninvkappa(1),meaninvkappa(2),100,Y0,Y,meanlogkappa(1),meanlogkappa(2));

    % update A
    Lambda_i = kron(chol(invSig,'lower'),sparse(1:T,1:T,sqrt(meaninvssq)));
    temp1 = zeros(T,n);

    for ii = 1:n
        A(:,ii) = 0;
        Zi = Lambda_i'*vec(Y - Z*A);
        iValpi = sparse(1:k,1:k,Hyper.iValp((ii-1)*k+1:ii*k));
        alpi0 = Hyper.alp0((ii-1)*k+1:ii*k);
        Wi = kron(B0(:,ii),Z)'*Lambda_i;
        Kalpi = iValpi + Wi*Wi';
        Kalp_VB(:,:,ii)=Kalpi;
        CKalpi = chol(Kalpi,'lower');

        alpi_hat = (CKalpi')\(CKalpi\(iValpi*alpi0 + Wi*Zi));
        A(:,ii) = alpi_hat;

        temp1(:,ii) = diag(Z*(Kalpi\Z'));

        Like_alp(ii) = -.5*(alpi_hat-alpi0)'*iValpi*(alpi_hat-alpi0) -.5*trace(CKalpi'\(CKalpi\iValpi)) ...
            -.5*sum(logV_Minn((ii-1)*k+1:ii*k)) - .5*ldet(Kalpi);

        tmp = Kalpi\speye(size(Kalpi,1));
        varalpi = diag(tmp);
        varalp(ii,:)=varalpi;
    end
    alp = reshape(A,k_alp,1); varA = reshape(varalp',k_alp,1);

    % update kappa
    bkappa(1) =  sum(((alp(idx_kappa1).^2+varA(idx_kappa1))./C_alp(idx_kappa1)));
    bkappa(2) =  sum(((alp(idx_kappa2).^2+varA(idx_kappa2))./C_alp(idx_kappa2)));
    Like_kappa = sum(Hyper.c0(:,1).*log(Hyper.c0(:,2))-gammaln(Hyper.c0(:,1))-ckappa./2.*(log(akappa./bkappa))+log(2)+logBesselK(-ckappa,sqrt(akappa.*bkappa))+...
        (Hyper.c0(:,1)-ckappa).*meanlogkappa-(Hyper.c0(:,2)-.5*akappa).*meankappa+.5*bkappa.*meaninvkappa);

    % update Sig
    nuhat = nusig+T;
    YZA = (Y-Z*A)'*diag(sqrt(meaninvssq));
    tmpYZA = YZA*YZA';

    tmpXKX = diag(sum(repmat(meaninvssq,1,n).*temp1));
    tmp2 = tmpYZA+tmpXKX;
    Psisighat = Psisig + tmp2;
    %Sig = Psisighat/(nuhat-n-1);
    invSig = (nuhat+n-1)*(Psisighat\speye(n));
    CiSigma = chol(invSig,'lower');

    %update s
    %s0
    ytstar = Y(tstar,:)'; Xtstar = kron(speye(n),Z(tstar,:));
    tmp = ytstar - Xtstar * alp;

    tmptr = sum(diag(invSig).*temp1(tstar,:)');
    phis0 = .5*(tmp'*invSig*tmp + tmptr);
    s02mean = phis0/(nus0-1);s02invmean = nus0/phis0;

    %s1
    ytstar1 = Y(tstar+1,:)'; Xtstar1 = kron(speye(n),Z(tstar+1,:));
    tmp = ytstar1 - Xtstar1 * alp;

    tmptr = sum(diag(invSig).*temp1(tstar+1,:)');
    phis1 = .5*(tmp'*invSig*tmp +  tmptr);

    s12mean = phis1/(nus1-1);s12invmean = nus1/phis1;
    Like_s12sq = -nus0*log(phis0) + gammaln(nus0)-nus1*log(phis1)+gammaln(nus1);

    yttmp = reshape(Y(tstar+2:end,:)',n*(T-(tstar+2)+1),1);
    Xttmp = X((tstar+2-1)*n+1:end,:);
    resid = yttmp - Xttmp*alp;

    ngridrho = length(rhogrid);

    %s2 and rho
    temp3 = repmat(diag(invSig)',T-(tstar+2)+1,1);
    ttindex = tstar+3:T;

    for kk = 1:length(s2grid)
        s2k = s2grid(kk);
        for jj = 1:ngridrho
            rhoj = rhogrid(jj);
            sttmp = stfunc(s2k,rhoj,ttindex);ststar2 = [s2k;sttmp'];
            temp2 = temp1(tstar+2:end,:)./repmat(ststar2.^2,1,n);
            EM1 = sum(log(sttmp));

            temp4 = temp3.*temp2;
            temptr = sum(sum(temp4,2));
            CSsqinv = sparse(1:(T-(tstar+2)+1), 1:(T-(tstar+2)+1),1./ststar2);
            C = kron(CSsqinv, CiSigma);
            Ctmp = C'*resid;
            tmpquadra = sum(Ctmp.^2);
            tmp = tmpquadra+temptr;
            logptheta(kk,jj) =  -(.5*n+1)*log(s2k^2)-n*EM1-.5*tmp+(alpha-1)*log(rhoj) + (beta-1)*log(1-rhoj)-betaln(alpha,beta);
        end
    end


    ptheta = exp(logptheta-max(logptheta(:)));
    Ctheta = -log(sum(sum(ptheta)))-max(logptheta(:));
    ptheta = ptheta./sum(sum(ptheta));
    means2 = sum(sum(ptheta,2).*s2grid);invmeans2sq = sum(sum(ptheta,2).*(s2grid.^(-2)));
    meanrho = sum(sum(ptheta).*rhogrid');
    invmeanstsq = E_stsq(s2grid,rhogrid,ptheta,stfunc,T,tstar);


    meaninvssq = [ones(tstar-1,1);[s02invmean,s12invmean,invmeans2sq]';invmeanstsq];


    % update VLB
    % Sig

    %tmpkr = kron(sparse(1:(T-(tstar)+1), 1:(T-(tstar)+1),meaninvssq(tstar:end)),invSig);
    temp2 = temp1(tstar:end,:).*repmat(meaninvssq(tstar:end),1,n);
    temp3 = repmat(diag(invSig)',T-tstar+1,1);
    temp4 = temp3.*temp2;
    temptr = sum(sum(temp4,2));

    yttmpsig = reshape(Y(tstar:end,:)',n*(T-(tstar)+1),1);
    Xttmpsig = X((tstar-1)*n+1:end,:);
    residsig = yttmpsig - Xttmpsig*alp;

    CSsqinv = sparse(1:(T-(tstar)+1), 1:(T-(tstar)+1),sqrt(meaninvssq(tstar:end)));
    C1 = kron(CSsqinv, CiSigma);
    Ctmp = C1'*residsig; tmpquadra = sum(Ctmp.^2);
    tmp = tmpquadra + temptr;
    Like_Sig = -.5*(nusig+T)*ldet(Psisighat)+.5*tmp;

    % constant
    C = -.5*n*T*log(2*pi) +.5*(n^2*p+n) + .5*nusig*ldet(Psisig)-.5*n*nusig*log(2)-gammalnn(nusig,n)-2*log(2)+...
        .5*(nusig+T)*n*log(2)+gammalnn((nusig+T),n);

    new_mllb = C - Ctheta + Like_kappa + sum(Like_alp) + Like_Sig + Like_s12sq ;

    %delta = abs(new_mllb - mllb)
    delta = new_mllb - mllb;
    %store_vlb(iteration_idx+1) = mllb;
    mllb = new_mllb;

    iteration_idx = iteration_idx + 1;
    if ( mod(iteration_idx, 100) == 0 )
        disp(  [ num2str(iteration_idx) ' iterations with delta = ', num2str(delta),'...' ] )
    end
    %disp(['E(rho) = ', num2str(meanrho), '; E(s0) = ' , num2str(sqrt(s02mean)),'; E(s1) = ',num2str(sqrt(s12mean)),'; E(s2) = ',num2str(means2)])
    rho=meanrho; s0=sqrt(s02mean);s1=sqrt(s12mean);s2=means2;
    delta_para = sum(abs([rho-rhot_1;s0-s0t_1;s1-s1t_1;s2-s2t_1]));
    rhot_1 = rho; s0t_1=s0;s1t_1=s1;s2t_1=s2;
end
time_spent = toc(start_time);


alp_VB = alp;A_VB = reshape(alp_VB,k,n)';
kappa_VB = meankappa;
rho_VB = meanrho; s0sq_VB = s02mean; s1sq_VB = s12mean; s2_VB = means2;
mllbsv = mllb;
disp(['VB estimation takes ',num2str(time_spent), ' seconds.'])

disp(['E(rho) = ', num2str(meanrho), '; E(s0) = ' , num2str(sqrt(s02mean)),'; E(s1) = ',num2str(sqrt(s12mean)),'; E(s2) = ',num2str(means2)])


lml_lenza=nan;lmlstd_lenza=nan;time_spent_ml=nan;
if ml_dummy
    M=10000;
    disp(['Starting ML Estimation...'])
    time_ml = tic;
    [lml_lenza,lmlstd_lenza] = ml_var_lenza_redu(X,Y,Y0,y,Z,p,rhogrid,s2grid,M,Hyper,alp_VB,Kalp_VB,...
        ptheta,bkappa,nuhat,Psisighat,nusig,Psisig,nus0,phis0,nus1,phis1,alpha,beta,tstar,stfunc);

    time_spent_ml = toc(time_ml);
    disp(['log-ML estimation takes ',num2str(time_spent_ml), ' seconds.'])
end
%end

% cd E:\VB\0519\VB\Xuewen
% VBapprox_VARSVredu
% disp(['SVO: ' num2str(mllbsvo) '; SV: ' num2str(mllb)])