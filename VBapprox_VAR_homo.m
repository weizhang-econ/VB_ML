% This script estimates the BVAR in reduced-form with no SV using VB
%clear all; %clc;
%function [A_VB,kappa_VB,rho_VB, s0sq_VB, s1sq_VB, s2_VB, mllbsv,lml_lenza,lmlstd_lenza,time_spent,time_spent_ml] = VBapprox_VAR_homo(Y,Y0,ml_dummy,p,tstar)
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
meansig = sparse(1:n,1:n,(mean(sig2,2)));
%nusig = median(diag(Psisig))/meansig+1+n;
nusig = n+2;
Psisig = (nusig-n-1)*meansig; 

alp = reshape(A,k_alp,1);
invSig =  sparse(1:n,1:n,(mean(1./sig2,2)));
varA = reshape(varalp',k_alp,1);
B0 = eye(n);
Like_alp = zeros(n,1);
% VB iteration starts here
TOL =  10; % tolerance level
mllb = -1e10; delta = 1e10;  start_time = tic; iteration_idx = 0;
disp(['Starting VB Estimation...'])
ckappa(1) = Hyper.c0(1,1)-n*p/2; akappa(1) = 2*Hyper.c0(1,2); bkappa(1) =  sum(((alp(idx_kappa1).^2+varA(idx_kappa1))./C_alp(idx_kappa1)));
ckappa(2)= Hyper.c0(2,1)-(n-1)*n*p/2; akappa(2) = 2*Hyper.c0(2,2); bkappa(2) =  sum(((alp(idx_kappa2).^2+varA(idx_kappa2))./C_alp(idx_kappa2)));
Kalp_VB = zeros(k,k,n);
while delta > TOL

    for jj = 1:2
        [meankappa(jj), meaninvkappa(jj), meanlogkappa(jj)] = exp_GIG(akappa(jj),bkappa(jj),ckappa(jj));
    end
    [~,Hyper.iValp,logV_Minn] = update_iValp(p,meaninvkappa(1),meaninvkappa(2),100,Y0,Y,meanlogkappa(1),meanlogkappa(2));

    % update A
    Lambda_i = kron(chol(invSig,'lower'),sparse(1:T,1:T,1));
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
    YZA = (Y-Z*A)';
    tmpYZA = YZA*YZA';

    tmpXKX = diag(sum(temp1));
    tmp2 = tmpYZA+tmpXKX;
    Psisighat = Psisig + tmp2;
    %Sig = Psisighat/(nuhat-n-1);
    invSig = (nuhat+n-1)*(Psisighat\speye(n));
    %CiSigma = chol(invSig,'lower');

    % update VLB
    % Sig
    Like_Sig = -.5*(nusig+T)*ldet(Psisighat);

    % constant
    C = -.5*n*T*log(2*pi) +.5*(n^2*p+n) + .5*nusig*ldet(Psisig)-.5*n*nusig*log(2)-gammalnn(nusig,n)+...
        .5*(nusig+T)*n*log(2)+gammalnn((nusig+T),n);

    new_mllb = C + Like_kappa + sum(Like_alp) + Like_Sig;

    %delta = abs(new_mllb - mllb)
    delta = new_mllb - mllb
    %store_vlb(iteration_idx+1) = mllb;
    mllb = new_mllb;

    iteration_idx = iteration_idx + 1;
    if ( mod(iteration_idx, 100) == 0 )
        disp(  [ num2str(iteration_idx) ' iterations with delta = ', num2str(delta),'...' ] )
    end
    
end
time_spent = toc(start_time);


alp_VB = alp;A_VB = reshape(alp_VB,k,n)';
kappa_VB = meankappa;
%rho_VB = meanrho; s0sq_VB = s02mean; s1sq_VB = s12mean; s2_VB = means2;
mllbsv = mllb;
disp(['VB estimation takes ',num2str(time_spent), ' seconds.'])

%disp(['E(rho) = ', num2str(meanrho), '; E(s0) = ' , num2str(sqrt(s02mean)),'; E(s1) = ',num2str(sqrt(s12mean)),'; E(s2) = ',num2str(means2)])


lml_homo=nan;lmlstd_homo=nan;time_spent_ml=nan;
if ml_dummy
    M=10000;
    disp(['Starting ML Estimation...'])
    time_ml = tic;
    [lml_homo,lmlstd_homo] = ml_var_homo_redu(X,Y,Y0,y,Z,p,M,Hyper,alp_VB,Kalp_VB,...
    bkappa,nuhat,Psisighat,nusig,Psisig);

    time_spent_ml = toc(time_ml);
    disp(['log-ML estimation takes ',num2str(time_spent_ml), ' seconds.'])
end
%end

% cd E:\VB\0519\VB\Xuewen
% VBapprox_VARSVredu
% disp(['SVO: ' num2str(mllbsvo) '; SV: ' num2str(mllb)])