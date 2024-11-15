# Replication Code:
# "Bayesian Model Comparison for Large Bayesian VARs after the COVID-19 Pandemic"
by
- Joshua C. C. Chan (Purdue University)
- Xuewen Yu (Fudan University)
- Wei Zhang (Purdue University)

This paper is under Major Revision at Journal of Econometrics, Special Issue. It can be found on [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4913626).

# Overview
This repository provides replication codes and input data to replicate the results shown in the application section of our paper. All core scripts are in the main directory:
- `Estimation.m`: the master file
- `VBapprox_VARSVminn_redu.m`: for VAR-SV model
- `VBapprox_VARSVO_sim.m`: for VAR-SVO model
- `VBapprox_VARSVt_sim.m`: for VAR-SVt model
- `VBapprox_VARLenza_simv2.m`: for VAR-CVD model
- `VBapprox_VAR_homo.m`: for VAR with homoskedasticity
- `data.mat`: the 180-variable dataset constructed using FRED-QD dataset
- `ccmm_data.mat`: the 16-variable dataset used by Carriero et al. (2024)
  
In addition, this repository contains the following subdirectory:
- `Utility`: all the useful functions as well as functions for estimating marginal likelihood

# Disclaimer
This code is free to use for academic purposes only, provided that the paper is cited as:

Chan, J. C., Yu, X., & Zhang, W. (2024). Bayesian Model Comparison for Large Bayesian VARs after the COVID-19 Pandemic. Available at SSRN 4913626.

This code comes without technical support of any kind. Under no circumstances will the author be held responsible for any use (or misuse) of this code in any way.

# References
Carriero, A., Clark, T. E., Marcellino, M., & Mertens, E. (2024). Addressing COVID-19 outliers in BVARs with stochastic volatility. Review of Economics and Statistics, 106(5), 1403-1417.

