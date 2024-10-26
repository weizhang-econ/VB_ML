%% Estimation for quarterly data

clear;
%rng(13242);
addpath('./utility')  
%load('./data_cleaning/data.mat');
%load('./data_cleaning/data_until032021.mat');
rng(154)
datasource = 1; % 1: own data; 2: Josh's data from JoE; 3: simulated data
ml_dummy=1;tstar=235;
CCMM=1;
for number_var = [104]
    switch datasource
        case 1
            load('data_annual.mat');
            %yt = yt(:,1:number_var);
            Y0 = yt(1:8,:);
            Y = yt(9:end,:);dates = dates(9:end);
            p=4;
        case 2
            data_all = load('macrodata_Q_2019Q4.csv');
            varid = [1,22,59,120,133,144,148,2,35,57,81,95,152,160,245]; % n = 15
            data = data_all(:,varid);
            Y0 = data(1:8,:);  % save the first 9 obs as the initial conditions
            Y = data(9:end,:);
            p=4;
        case 3
            outlier=true;qerror = false;lenza =false;tstar=80;
            [A_true,B0_true,h_true,sig2_true,Q2_true,O_true,probo_true,su_true,Sigma_true] = gen_para(T,n,p,outlier,qerror,lenza,tstar);
            %[Y0,Y] = gendata_VARSVOt_redu(T,n,p,A_true,B0_true,h_true,Q2_true,O_true,lenza,su_true,Sigma_true);
    end
    if CCMM
        varnames = [{'Real Income'},{'Real Consumption'},{'IP'},{'Capacity Utilization'},{'Unemployment'},...
            {'Nonfarm Payrolls'},{'Hours'},{'Hours Earnings'},{'PPI'},{'PCE'},{'Housing Starts'},{'S&P 500'},...
            {'USD/GBP FX Rate'},{'5-Year Yield'},{'10-Year Yield'},{'Baa Spead'}];
        ccmm_series = [{'GDPC1'},{'PCECC96'},{'INDPRO'},{'CUMFNS'},{'UNRATE'},...
            {'PAYEMS'},{'CES0600000007'},{'CES0600000008'},{'WPSFD49207'},{'PCEPILFE'},{'HOUST'},...
            {'S&P 500'},{'EXUSUKx'},{'GS5'},{'GS10'},{'BAAFFM'}];
        load_index = ismember(series,ccmm_series);
        yt_ccmm = yt(:,load_index);series_ccmm = series(load_index);
        load('Baaspread.mat')
        yt = [yt_ccmm,BAAFFM(2:end)];
        Y0 = yt(1:8,:);
        Y = yt(9:end,:);
        p=4;
    end
    for model = [1,2,3,9]
        switch model
            case 1
                model_name = 'VAR-SV';
                disp(['Starting VB for ' model_name '...']);
                [A_SV,B0_SV,h_SV,h0_SV,isigh2_SV,kappa_SV,~,mllb_SV,ml_SV,mlstd_SV,time_SV,time_SV_ml] = VBapprox_VARSVminn_redu(Y,Y0,ml_dummy,p);
                if ml_dummy ~= 0
                    fprintf(['log marginal likelihood of ' model_name ': %.1f \n'], ml_SV);
                end
            case 2
                model_name = 'VAR-SVO';
                disp(['Starting VB for ' model_name '...']);
                [A_SVO,B0_SVO,h_SVO,h0_SVO,isigh2_SVO,kappa_SVO,mllbsvo,po_SVO,O_SVO,ml_SVO,mlstd_SVO,time_SVO,time_SVO_ml] = VBapprox_VARSVO_sim(Y,Y0,ml_dummy,p);
                if ml_dummy ~= 0
                    fprintf(['log marginal likelihood of ' model_name ': %.1f \n'], ml_SVO);
                end
            case 3
                model_name = 'VAR-SVt';
                disp(['Starting VB for ' model_name '...']);
                [A_SVt,B0_SVt,h_SVt,h0_SVt,isigh2_SVt,kappa_SVt,q2_SVt,mllbsvt,ml_SVt,mlstd_SVt,time_SVt,time_SVt_ml] = VBapprox_VARSVt_sim(Y,Y0,ml_dummy,p);
                if ml_dummy ~= 0
                    fprintf(['log marginal likelihood of ' model_name ': %.1f \n'], ml_SVt);
                end
            case 4
                model_name = 'VAR-SVOt';
                disp(['Starting VB for ' model_name '...']);
                [A_SVOt,B0_SVOt,h_SVOt,h0_SVOt,isigh2_SVOt,kappa_SVOt,q2_SVOt,po_SVOt,O_SVOt,mllbsvot,ml_SVOt,mlstd_SVOt,time_SVOt] = VBapprox_VARSVOt_sim(Y,Y0,ml_dummy);
                if ml_dummy ~= 0
                    fprintf(['log marginal likelihood of ' model_name ': %.1f \n'], ml_SVOt);
                end
            case 5
                model_name = 'VAR-SV-fixed-kappa';
                disp(['Starting VB for ' model_name '...']);
                [A_SV,B0_SV,h_SV,h0_SV,isigh2_SV,sigh2_SV,mllb_SV,lml_SV,lmlstd_SV,time_SV,time_SV_ml] = VBapprox_VARSVminn_redu_fixedkappa(Y,Y0,ml_dummy);
            case 6
                model_name = 'VAR-SVO-fixed-kappa';
                disp(['Starting VB for ' model_name '...']);
                [A_SVO,B0_SVO,h_SVO,h0_SVO,isigh2_SVO,mllbsvo,po_SVO,O_SVO,lml_SVO,lmlstd_SVO,time_SVO,time_SVO_ml] = VBapprox_VARSVO_sim_fixkappa(Y,Y0,ml_dummy);
            case 7
                model_name = 'VAR-SVt-fixed-kappa';
                disp(['Starting VB for ' model_name '...']);
                [A_SVt,B0_SVt,h_SVt,h0_SVt,isigh2_SVt,q2_SVt,mllbsvt,lml_SVt,lmlstd_SVt,time_spent,time_spent_ml] = VBapprox_VARSVt_sim_fixkappa(Y,Y0,ml_dummy);
            case 8
                model_name = 'VAR-SVO-common-outliers';
                disp(['Starting VB for ' model_name '...']);
                [A_SVOc,B0_SVOc,h_SVOc,h0_SVOc,isigh2_SVOc,kappa_SVOc,mllbsvoc,po_SVOc,O_SVOc,ml_SVOc,mlstd_SVOc,time_SVOc,time_SVO_mlc] = VBapprox_VARSVO_sim(Y,Y0,ml_dummy,p);
                if ml_dummy ~= 0
                    fprintf(['log marginal likelihood of ' model_name ': %.1f \n'], ml_SVOc);
                end
            case 9
                model_name = 'VAR-lenza';
                disp(['Starting VB for ' model_name '...']);
                [A_VBl,kappa_VBl,rho_VBl, s0sq_VBl, s1sq_VBl, s2_VBl, mllbsvl,lml_lenza,lmlstd_lenza,time_spentl,time_spent_mll] = VBapprox_VARLenza_simv2(Y,Y0,ml_dummy,p,tstar);
                if ml_dummy ~= 0
                    fprintf(['log marginal likelihood of ' model_name ': %.1f \n'], lml_lenza);
                end
        end
    end

    filename = ['first',num2str(number_var),'_p12_noLenza.mat'];

    parsave(filename,Y, Y0,...
        A_SV,B0_SV,h_SV,h0_SV,isigh2_SV,kappa_SV,mllb_SV,ml_SV,mlstd_SV,time_SV,time_SV_ml,...
        A_SVO,B0_SVO,h_SVO,h0_SVO,isigh2_SVO,kappa_SVO,mllbsvo,po_SVO,O_SVO,ml_SVO,mlstd_SVO,time_SVO,time_SVO_ml,...
        A_SVt,B0_SVt,h_SVt,h0_SVt,isigh2_SVt,kappa_SVt,q2_SVt,mllbsvt,ml_SVt,mlstd_SVt,time_SVt,time_SVt_ml);

end
% tab = [ml_SV,ml_SVO,ml_SVt,lml_lenza;mlstd_SV,mlstd_SVO,mlstd_SVt,lmlstd_lenza];
% VLB = [mllb_SV;mllbsvo;mllbsvt;mllbsvl];

