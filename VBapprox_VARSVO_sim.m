% This script estimates the BVAR in reduced-form with stationary autoregressive SV using VB
%clear all; %clc;
%function [A_VB,B0_VB,h_VB,h0_VB,isigh2_VB,kappa_VB,mllbsvo,po_VB,O_VB,lml_SVO,lmlstd_SVO,time_spent,time_spent_ml] = VBapprox_VARSVO_sim(Y,Y0,ml_dummy,p)
[T,n] = size(Y);
%y = reshape(Y',n*T,1);
tmpY = [Y0(end-p+1:end,:); Y];
X = zeros(T,n*p);
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];    
k_beta = n*(n-1)/2;       % dimension of B0
k_alp = n^2*p + n;        % dimension of A
k = k_alp/n;

% priors
Vh0 = 10*ones(n,1);
nuh = 3*ones(n,1); Sh = 0.1*ones(n,1).*(nuh-1);
%nuh = (n+3)*ones(n,1); Sh = nuh*0.15*12/p;
%nuh = nuh/2; Sh = Sh/2;
%kappa1 = .2^2; kappa4 = .2^2;kappa2=(.2^2)^2;kappa3=100;
kappa1 = .2^2; kappa3 = 100; kappa4 = .2^2;kappa2 = (.2^2)^2; 
[Hyper.alp0,Hyper.Valp,sig2_hat,U_hat] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
[Hyper.beta0,Hyper.Vbeta] = prior_B0(Y0,Y,kappa4,p);
Hyper.iValp = 1./Hyper.Valp;
Hyper.iVbeta = 1./Hyper.Vbeta;
Hyper.c0 = [1,1/.2^2; 1,1/.2^2; 1,1]; 
[C_alp,idx_kappa1,idx_kappa2] = get_C(n,p,sig2_hat);
apo = 10/4*ones(n,1);
bpo = 40*(1-1/16)*ones(n,1);
%  apo = 1 / ( 4 * 12) * 10 * 12*ones(n,1);
%  bpo = (10 * 12 - apo)*(1-1/48);
% % 10 years
%  apo = 1 / ( 10 * 12) * 10 * 12*ones(n,1);
%  bpo = (10 * 12 - apo)*(1-1/120);

ngrid = 20;
ogrid = [1;linspace(2,20,ngrid-1)'];

% initialize for storage
h_hat = zeros(T,n);
A = zeros(n,k); varalp = zeros(n,k);
akappa = zeros(3,1);bkappa = zeros(3,1);ckappa = zeros(3,1);
meaninvkappa = (.2^2)^(-1)*ones(3,1);meanlogkappa = log(.2^2)*ones(3,1);meankappa = .2^2*ones(3,1);
ohat = ones(T,n);
io2hat = ones(T,n);T1bar = T*ones(n,1);
%pohat = 1/48*ones(n,1);
TOL =  10; % tolerance level

% initialize
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T); HH = H'*H;
XX = X'*X; 
ih_and_d = h_hat;beta_approxvar = 0.1*ones(k_beta,1);
for ii=1:n
    iValpi = sparse(1:k,1:k,1./Hyper.Valp((ii-1)*k+1:ii*k));
    A(ii,:) = (XX + iValpi)\(X'*Y(:,ii));
    s2i = (Y(:,ii) - X*A(ii,:)').^2;
	h0i_tini = log(mean(s2i)); sig2hi_tini = .04; HiSH_h = HH/sig2hi_tini; 
    s2_ini = max(s2i,1e-4); Khi = HiSH_h + 1/4.9*speye(T); d_iKhi = band1inv(Khi);
    hi_hat = Khi\(HiSH_h*h0i_tini*ones(T,1) + (log(s2_ini) - 1.27)/4.9);
    h_hat(:,ii) = hi_hat; 	
	ih_and_d(:,ii) = -hi_hat + .5*d_iKhi;
    varalp(ii,:) = diag(var(A(ii,:)));
end
varA = reshape(varalp',k_alp,1);
nuh_hat = nuh + T/2;
Sh_hat = .1 * (nuh + T/2 - 1);
Esigh2 = nuh_hat./Sh_hat; 
Kh0_hat = ones(n,1)./Vh0 + Esigh2;
h0_hat = Kh0_hat.\(Esigh2.*h_hat(1,:));
alp = reshape(A',k_alp,1);  beta = 0.1*ones(k_beta,1);
XA = X*A'; 
B0_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
B0 = eye(n);B0(B0_id) = beta;
C_beta = Hyper.Vbeta/kappa4;
%logk = @(nu,zz) 1/2*log(pi./(2*nu))-nu.*log(exp(1)*zz./(2*nu));

% VB iteration starts here 
B0_Var = zeros(n,n); B0_CoVar = zeros(n,n,n); 
mllb = -1e10; delta = 1e10; start_time = tic; iteration_idx = 0;

ckappa(1) = Hyper.c0(1,1)-n*p/2; akappa(1) = 2*Hyper.c0(1,2); bkappa(1) =  sum(((alp(idx_kappa1).^2+varA(idx_kappa1))./C_alp(idx_kappa1)));
ckappa(2)= Hyper.c0(2,1)-(n-1)*n*p/2; akappa(2) = 2*Hyper.c0(2,2); bkappa(2) =  sum(((alp(idx_kappa2).^2+varA(idx_kappa2))./C_alp(idx_kappa2)));
ckappa(3) = Hyper.c0(3,1)-n*(n-1)/4; akappa(3) = 2*Hyper.c0(3,2); bkappa(3) =  sum((beta.^2+beta_approxvar)./C_beta);
% for jj = 1:3
%     [meankappa(jj), meaninvkappa(jj), meanlogkappa(jj)] = exp_GIG(akappa(jj),bkappa(jj),ckappa(jj),kappagrid);
% end
% [~,Hyper.iValp,logV_Minn] = update_iValp(p,meaninvkappa(1),meaninvkappa(2),100,Y0,Y,meanlogkappa(1),meanlogkappa(2));
% [~,Hyper.iVbeta,logVbeta] = update_iVbeta(Y0,Y,meaninvkappa(3),meanlogkappa(3));

%while delta > TOL * abs(mllb)
while delta > TOL
	% update kappa
    for jj = 1:3
        [meankappa(jj), meaninvkappa(jj), meanlogkappa(jj)] = exp_GIG(akappa(jj),bkappa(jj),ckappa(jj));
    end
    [~,Hyper.iValp,logV_Minn] = update_iValp(p,meaninvkappa(1),meaninvkappa(2),100,Y0,Y,meanlogkappa(1),meanlogkappa(2));
    [~,Hyper.iVbeta,logVbeta] = update_iVbeta(Y0,Y,meaninvkappa(3),meanlogkappa(3),p);

        % update alp
    Ytilde = Y*sparse(B0');
    G = zeros(n,n); Gtilde = zeros(T,n); 
	Like_alp = zeros(n,1); %Like_alp will be used in VLB calculation
  %  MiKalpi = zeros(k,k,n);
    Kalp_VB = zeros(k,k,n);
    for ii = 1:n            
        XAB0 = XA(:,[1:ii-1 ii+1:n])*sparse(B0(:,[1:ii-1 ii+1:n])');
        Lambda_i = vec(io2hat(:,ii:n).^(1/2).*exp(ih_and_d(:,ii:n)/2));
        Zi = vec(Ytilde(:,ii:n) - XAB0(:,ii:n)).*Lambda_i;    
        %iValpi = sparse(1:k,1:k,1./Hyper.Valp((ii-1)*k+1:ii*k));
        iValpi = sparse(1:k,1:k,Hyper.iValp((ii-1)*k+1:ii*k));
        alpi0 = Hyper.alp0((ii-1)*k+1:ii*k);
        Wi = kron(B0(ii:n,ii),X).*Lambda_i; 
		Wi_v = kron(B0_Var(ii:n,ii).^(1/2),X).*Lambda_i; 
        Kalpi = iValpi + Wi'*Wi+Wi_v'*Wi_v;		
        Kalp_VB(:,:,ii)=Kalpi;
        CKalpi = chol(Kalpi,'lower'); 

        Lambdaf_i = vec(io2hat.^(1/2).*exp(ih_and_d./2));
		%Wi_cov = kron(ones(n,1),X).*Lambdaf_i;
        Wi_cov = repmat(X,n,1).*Lambdaf_i;
		B0_CoVar_i = squeeze(B0_CoVar(ii,:,:)); Ai0 = A'; Ai0(:,ii) = zeros(k,1);
		Zi_cov = vec(((Y - X*Ai0)*B0_CoVar_i')).*Lambdaf_i;
        alpi_hat = (CKalpi')\(CKalpi\(iValpi*alpi0 + Wi'*Zi + Wi_cov'*Zi_cov));        
        A(ii,:) = alpi_hat;
        XA(:,ii) = X*alpi_hat;
		temp1 = diag(X/Kalpi*X'); 
		G(ii,:) = temp1'*(io2hat.*exp(ih_and_d)); 
        Gtilde(:,ii) = temp1; 
        Like_alp(ii) = -.5*(alpi_hat-alpi0)'*iValpi*(alpi_hat-alpi0) -.5*trace(CKalpi'\(CKalpi\iValpi)) ...
            -.5*sum(logV_Minn((ii-1)*k+1:ii*k)) - .5*ldet(Kalpi);
        % 		Like_alp(ii) = -.5*(alpi_hat-alpi0)'*iValpi*(alpi_hat-alpi0) -.5*trace(CKalpi'\(CKalpi\iValpi)) ...
        % 		               +.5*ldet(iValpi) - .5*ldet(Kalpi);

        tmp = Kalpi\speye(size(Kalpi,1));
        varalpi = diag(tmp);
        varalp(ii,:)=varalpi;
    end
    alp = reshape(A',k_alp,1); varA = reshape(varalp',k_alp,1);

        % update beta
    E = Y - XA;
    count_beta = 0; beta_approxvar = zeros(k_beta,1); B0_CoVar = zeros(n,n,n); Gbreve = zeros(T,n);Kbeta_VB = zeros(k_beta,k_beta);
    GiKbiinv=zeros(n,n);
    %B0_Var = zeros(n,n); 
	Like_beta = zeros(n,1); %Like_beta will be used in VLB calculation
    for ii=2:n
        X_betai = -E(:,1:ii-1);   
        %iD = sparse(1:T,1:T,exp(-h_and_d(:,ii)));
        iM = sparse(1:T,1:T,io2hat(:,ii).*exp(ih_and_d(:,ii)));
        %iVbetai = sparse(1:ii-1,1:ii-1,1./Hyper.Vbeta(count_beta+1:count_beta+ii-1)); 
        iVbetai = sparse(1:ii-1,1:ii-1,Hyper.iVbeta(count_beta+1:count_beta+ii-1));
		SGiS = diag(G(1:ii-1,ii));       
        Kbetai = iVbetai + X_betai'*iM*X_betai+SGiS; 
		iKbetai = Kbetai\speye(size(Kbetai,1)); 
        %betai_id = nonzeros(tril(reshape(1:(ii-1)^2,ii-1,ii-1))');
        B0_CoVar(1:(ii-1),ii,1:(ii-1)) = iKbetai;Kbeta_VB(count_beta+1:count_beta+ii-1,count_beta+1:count_beta+ii-1)=Kbetai;
        iKB0i = squeeze(B0_CoVar(:,ii,:)); Gbreve(:,ii) = diag(E*iKB0i*E');GiKbiinv(:,ii) = G(:,ii).*diag(iKB0i);
        % 		Gie = zeros(ii-1,1); Gie(ii-1) = G(ii-1,ii-1);
        %         betai_hat = iKbetai*(X_betai'*iM*E(:,ii)+Gie);
        %         Gie = zeros(ii-1,1); Gie(ii-1) = G(ii-1,ii);% Wei updated on 07/18
        %         betai_hat = Kbetai\(X_betai'*iM*E(:,ii)-Gie);% Wei updated on 07/18
        betai_hat = iKbetai*(X_betai'*iM*E(:,ii));
        beta(count_beta+1:count_beta+ii-1) = betai_hat;
        beta_approxvar(count_beta+1:count_beta+ii-1) = diag(iKbetai);
        Like_beta(ii) = -.5*betai_hat'*iVbetai*betai_hat -.5*trace(iKbetai*iVbetai) ...
            -.5*sum(logVbeta(count_beta+1:count_beta+ii-1)) + .5*ldet(iKbetai);
%         Like_beta(ii) = -.5*betai_hat'*iVbetai*betai_hat -.5*trace(iKbetai*iVbetai) ...
% 		                +.5*ldet(iVbetai) + .5*ldet(iKbetai);		
        count_beta = count_beta + ii-1;
    end
    B0(B0_id) = beta;
    B0_Var(B0_id) = beta_approxvar;
    
        % update h & sigh2
    B0E = E*sparse(B0'); 
    nuh_hat = nuh + T/2; Sh_hat = zeros(n,1);
	ih_and_d = zeros(T,n); %store -hi_hat+.5*d_iKhi
    Like_h = zeros(n,1); %Like_h will be used in VLB calculation	
	Like_h0 = zeros(n,1); %Like_h0 will be used in VLB calculation
    s2mat = zeros(T,n);Kh_VB = zeros(T,T,n);
    %     switch h_approx
    %         case 3
    for ii = 1:n
        % h
        s2_1 = B0E(:,ii).^2;
        s2_2 = Gtilde*((B0(ii,:)').^2);
        s2_3 = Gbreve(:,ii);
        s2 = s2_1 + s2_2 + s2_3;
        s2mat(:,ii)=s2;
        io2hats2 = io2hat(:,ii).*s2;
        Eisigh2 = Esigh2(ii); h0i_hat = h0_hat(ii);
        [hi_hat,Khi_hat,d_iKhi,d1_iKhi] = update_h_minKL(io2hats2,h0i_hat,Eisigh2);
        %[hi_hat,Khi_hat,d_iKhi,d1_iKhi] = update_h_LocalApprox(s2,h0i_hat,Eisigh2);
        %[hi_hat,Khi_hat,d_iKhi,d1_iKhi] = update_h_NormalApprox(s2,h0i_hat,Eisigh2);
        h_hat(:,ii) = hi_hat;Kh_VB(:,:,ii) = Khi_hat;
        ih_and_d(:,ii) = -hi_hat+.5*d_iKhi;
        Like_h(ii) = -.5*sum(hi_hat) -.5*ldet(Khi_hat);

        % sigh2
        Kh0i_hat = Kh0_hat(ii);
        trHHiKhi_hat = d_iKhi'*[2*ones(T-1,1);1] - 2*sum(d1_iKhi); % same as trace(HH/Khi_hat)
        Sh_hat(ii) = Sh(ii) + .5*((hi_hat-h0i_hat)'*HH*(hi_hat-h0i_hat) + 1/Kh0i_hat + trHHiKhi_hat);
        Like_h0(ii) = -.5*log(Vh0(ii)) -.5*log(Kh0i_hat);
    end
%     end 
	Esigh2 = nuh_hat./Sh_hat; 
	
       % update h0
    Kh0_hat = ones(n,1)./Vh0 + Esigh2; 
    h0_hat = 1./Kh0_hat.*(Esigh2.*h_hat(1,:)');   
    Like_sighi = -(h0_hat.^2+1./(Kh0_hat))./(2*Vh0) ...
	             + nuh.*log(Sh) - gammaln(nuh) ...
				 - (nuh+T/2).*log(Sh_hat) + gammaln(nuh+T/2);  %Like_h0 will be used in VLB calculation	
	
    % update po
    alppo = apo+(T-T1bar); betapo = bpo+T1bar;
    logpo = psi(alppo)-psi(alppo+betapo);
    log1mpo = psi(betapo)-psi(alppo+betapo);
   
    % update o2
    Co = zeros(T,n);T1bar = zeros(n,1);Like_o2po = zeros(n,1);opost_VB = zeros(T,ngrid,n);
    
    for ii = 1:n
        log1mpoi = log1mpo(ii);
        logpoi = logpo(ii);
        opri = [log1mpoi; repmat(logpoi-log(ngrid-1),ngrid-1,1)];
        uio = exp(ih_and_d(:,ii)).*s2mat(:,ii);
        ollikeli = -repmat(log(ogrid'),T,1)+repmat(-.5*uio,1,ngrid)./repmat(ogrid'.^2,T,1);
        opost = exp(ollikeli + repmat(opri',T,1)-repmat(max(ollikeli,[],2),1,ngrid));
        %opost = exp(ollikeli + repmat(opri',T,1)-repmat(median(ollikeli,2),1,ngrid));
        %opost = exp(ollikeli + repmat(opri',T,1));
        %Co(:,ii)=1./sum(opost,2);
        %Co(:,ii)=-log(sum(opost,2));
        Co(:,ii)=-log(sum(exp(ollikeli + repmat(opri',T,1)),2));
        opost = opost./repmat(sum(opost,2),1,ngrid);opost_VB(:,:,ii)=opost;
        io2hat(:,ii) = sum(repmat(ogrid'.^(-2),T,1).*opost,2);
        ohat(:,ii)= sum(repmat(ogrid',T,1).*opost,2);
        T1bar(ii) = sum(opost(:,1));
        Like_o2po(ii) = sum(-Co(:,ii)+.5*io2hat(:,ii).*exp(ih_and_d(:,ii)).*s2mat(:,ii))...
            -(T-T1bar(ii))*logpo(ii)-T1bar(ii)*log1mpo(ii)...
            +gammaln(apo(ii)+T-T1bar(ii))+gammaln(bpo(ii)+T1bar(ii))-gammaln(apo(ii)+bpo(ii)+T)...
            +gammaln(apo(ii)+bpo(ii))-gammaln(apo(ii))-gammaln(bpo(ii));
    end

%     % update kappa
    bkappa(1) =  sum(((alp(idx_kappa1).^2+varA(idx_kappa1))./C_alp(idx_kappa1)));
    bkappa(2) =  sum(((alp(idx_kappa2).^2+varA(idx_kappa2))./C_alp(idx_kappa2)));
    bkappa(3) =  sum((beta.^2+beta_approxvar)./C_beta);
%     for jj = 1:3
%         [meankappa(jj), meaninvkappa(jj), meanlogkappa(jj)] = exp_GIG(akappa(jj),bkappa(jj),ckappa(jj),kappagrid);
%     end
%     [~,Hyper.iValp,logV_Minn] = update_iValp(p,meaninvkappa(1),meaninvkappa(2),100,Y0,Y,meanlogkappa(1),meanlogkappa(2));
%     [~,Hyper.iVbeta,logVbeta] = update_iVbeta(Y0,Y,meaninvkappa(3),meanlogkappa(3));

    % update kappa likelihood
    % %    symbolic computation from MATLAB for log(besselk)
    %         Like_kappa = sum(Hyper.c0(:,1).*log(Hyper.c0(:,2))-gammaln(Hyper.c0(:,1))-ckappa./2.*(log(akappa./bkappa))+log(2)+double(log(besselk(sym(ckappa), sym(sqrt(akappa.*bkappa)))))+...
    %            (Hyper.c0(:,1)-ckappa).*meanlogkappa-(Hyper.c0(:,2)-.5*akappa).*meankappa+.5*bkappa.*meaninvkappa);

    % approximating log(besselk) based on equation 10.41.2 from the website
    % https://dlmf.nist.gov/10.41
    Like_kappa = sum(Hyper.c0(:,1).*log(Hyper.c0(:,2))-gammaln(Hyper.c0(:,1))-ckappa./2.*(log(akappa./bkappa))+log(2)+logBesselK(-ckappa,sqrt(akappa.*bkappa))+...
        (Hyper.c0(:,1)-ckappa).*meanlogkappa-(Hyper.c0(:,2)-.5*akappa).*meankappa+.5*bkappa.*meaninvkappa);

    % update the variational lower bound
    new_mllb = Like_kappa+sum(Like_alp) + sum(Like_beta) + sum(Like_h) + sum(Like_h0) + sum(Like_sighi) + sum(Like_o2po) - T*n/2*log(2*pi);
    %  new_mllb = sum(Like_alp) + sum(Like_beta) + sum(Like_h) + sum(Like_h0) + sum(Like_sighi)  - T*n/2*log(2*pi);

    for ii = 1:n
        iMbari = exp(ih_and_d(:,ii)).*io2hat(:,ii); B0Ei = B0E(:,ii); Gi = G(:,ii); B0i = B0(ii,:)'; Gbrevei = Gbreve(:,ii);
        %new_mllb = new_mllb  - .5*iMbari'*B0Ei.^2 - .5*Gi'*B0i.^2 - .5*iMbari'*Gbrevei + (T+n*p+ii+1)/2;
        new_mllb = new_mllb  - .5*iMbari'*B0Ei.^2 - .5*Gi'*B0i.^2 - .5*iMbari'*Gbrevei-.5*sum(GiKbiinv(:,ii)) + (T+n*p+ii+1)/2;
    end
    delta = abs(new_mllb - mllb);
    %store_vlb(iteration_idx+1) = mllb;
    mllb = new_mllb;
    
    iteration_idx = iteration_idx + 1;
    if ( mod(iteration_idx, 100) == 0 )
        disp(  [ num2str(iteration_idx) ' iterations with delta = ', num2str(delta),'...' ] )
    end
    
end
time_spent = toc(start_time);

po_VB = alppo./(alppo+betapo);
A_VB = A;alp_VB = alp;
B0_VB = B0;beta_VB=beta;
h0_VB = h0_hat;Kh0_VB = Kh0_hat;
h_VB = h_hat;
isigh2_VB = Esigh2;
kappa_VB = meankappa;
O_VB = ohat;
mllbsvo = mllb;
disp(['VB estimation takes ',num2str(time_spent), ' seconds.'])
lml_SVO=nan;lmlstd_SVO=nan;time_spent_ml=nan;
if ml_dummy
M=10000;
disp(['Starting ML Estimation...'])
time_ml = tic;
[lml_SVO,lmlstd_SVO] = ml_var_svo_redu(X,Y,Y0,M,Hyper,alp_VB,Kalp_VB,beta_VB,Kbeta_VB,...
    h_VB,Kh_VB,h0_VB,Kh0_VB,opost_VB,ogrid,alppo,betapo,bkappa,Vh0,nuh,Sh,apo,bpo);
time_spent_ml = toc(time_ml);
disp(['log-ML estimation takes ',num2str(time_spent_ml), ' seconds.'])
end
%end
% cd E:\VB\0519\VB\Xuewen
% VBapprox_VARSVredu
% disp(['SVO: ' num2str(mllbsvo) '; SV: ' num2str(mllb)])