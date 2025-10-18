# HYPK-NatA-Simulations
This repository contains MATLAB scripts used for the kinetic modeling and quantitative analysis of HYPK–NatA-mediated cotranslational N-terminal acetylation (NTA), as described in the study: “HYPK promotes N-terminal acetylation through rapid ribosome exchange of NatA.”
--------------------------------------------------------------------
This repository contains three MATLAB (.m) files:

1. HYPK_NatA_Ac_analytical.m:  Analytical model that computes the probability of N-terminal acetylation by solving the master equation for NatA–ribosome–HYPK kinetics as a function of elongation and binding rates.
---------------------------------------------------------------------
2. HYPK_NatA_Ac_sim_model.m: Stochastic model that simulates cotranslational N-terminal acetylation whose input parameters are NatA RNC binding kinteics, enzymatic rate and initiation, codon translation and termination rate.
--------------------------------------------------------------------
3. HYPK_NatA_k23_fit_model.m:  This code fits the mathematical function to the k3 experimental data.
