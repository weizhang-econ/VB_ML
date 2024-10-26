# Replication Code:
# "Bayesian Model Comparison for Large Bayesian VARs after the COVID-19 Pandemic"
by
- Joshua C. C. Chan (Purdue University)
- Xuewen Yu (Fudan University)
- Wei Zhang (Purdue University)

This paper is under review and it can be found on [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4913626).

# Overview
This repository provides replication codes and raw input data to replicate the results shown in our paper. All core scripts are in the main directory:
- `Estimation.m`: the master file
- `VBapprox_VARSVminn_redu.m`: VAR-SV
- `VBapprox_VARSVO_sim.m`: VAR-SVO
- `VBapprox_VARSVt_sim.m`: VAR-SVt
- `VBapprox_VARLenza_simv2.m`: VAR-CVD
- `VBapprox_VAR_homo.m`: VAR
- `data.mat`: dataset used
  
In addition, this repository contains the following subdirectory:
- `Utility`: all the useful functions as well as functions for estimating marginal likelihood

# Disclaimer
This code comes without technical support of any kind. Under no circumstances will the author be held responsible for any use (or misuse) of this code in any way.
