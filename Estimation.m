%% Estimation

clear;
restoredefaultpath;
addpath('./utility')

rng(154)
ml_dummy=1; %1: estimate marginal likelihood; 2: not estimate marginal likelihood

datasource = 1; %1: 18-variable dataset from CCMM; 2: 180-variable dataset constructed from FRED-QD
switch datasource
    case 1
        load('data_ccmm.mat');
        yt = data;
        p=12;tstar = 724;
    case 2
        load('data.mat');
        p=4;tstar = 243;
end

%whole sample
Y0 = yt(1:8,:);
Y = yt(9:end,:);dates = dates(9:end);

model_name = 'VAR-SV';
disp(['Starting VB for ' model_name '...']);
VBapprox_VARSVminn_redu;
save('./Results/VARSV_Full.mat')

clearvars -except dates ml_dummy series Y0 Y p
model_name = 'VAR-no-SV';
disp(['Starting VB for ' model_name '...']);
VBapprox_VAR_homo;
save('./Results/VARnoSV_Full.mat')

clearvars -except dates ml_dummy series Y0 Y p
model_name = 'VAR-SVO';
disp(['Starting VB for ' model_name '...']);
VBapprox_VARSVO_sim;
save('./Results/VARSVO_Full.mat')

clearvars -except dates ml_dummy series Y0 Y p
model_name = 'VAR-SVt';
disp(['Starting VB for ' model_name '...']);
VBapprox_VARSVt_sim;
save('./Results/VARSVt_Full.mat')

clearvars -except dates ml_dummy series Y0 Y p
model_name = 'VAR-lenza';
disp(['Starting VB for ' model_name '...']);
VBapprox_VARLenza_simv2;
save('./Results/VARSVD_Full.mat')
