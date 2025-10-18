# HYPK-NatA-RNC Simulations
--------------------------------------------------------------------
This repository contains MATLAB scripts used for the kinetic modeling and quantitative analysis of HYPK–NatA-mediated cotranslational N-terminal acetylation (NTA), as described in the study: “HYPK promotes N-terminal acetylation through rapid ribosome exchange of NatA.”
This repository contains three MATLAB (.m) files:

1. HYPK_NatA_Ac_analytical.m:  Analytical model that computes the probability of different nascent protein states by solving the master as a function of NatA RNC binding kinteics, enzymatic rate and codon translation.
---------------------------------------------------------------------
2. HYPK_NatA_Ac_sim_model.m: Stochastic model that simulates cotranslational N-terminal acetylation whose input parameters are NatA RNC binding kinteics, enzymatic rate and initiation, codon translation and termination rate.
--------------------------------------------------------------------
3. HYPK_NatA_k23_fit_model.m:  This code fits the mathematical function to the k3 experimental data.
