Scripts and data sets for Gerhardt et al 2021

I. TetR based IP:

fig1ddata 

> experimental data for the atc-TetR IP: inducer concentration, mean, CV

Tetgfpfit, tetgfperr, 0206.mat

> Deterministic model solution, fitting and parameter files

tetgfp2.sbproj, tetgfp_SSA

> Stochastic model SimBiology file and script to call and simulate the SimBiology model

Stats_0327.mat

> Final SSA runs for all 3 plasmid copy numbers. Can be unpacked to plot final mean and CV figures.

II. LuxR based IP:

luxrdata.m

> experimental data for all AHL concentrations

Luxsimplode, luxR_simplesteadystate

> Deterministic model solutions

directMethod.m
firstReactionMethod.m

> files for running Gillespie SSA (add source- https://github.com/nvictus/Gillespie)

luxrpropensities 

> Stochastic model; propensity function for Gillespie SSA

luxrsim_exn 

> Gillespie simulation script with extrinsic noise

revisions_*

> Scripts for comparing deterministic and stochastic model parameters for both IPs, LuxR and TetR based.