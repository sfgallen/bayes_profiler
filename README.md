# bayes_profiler
A Bayesian Markov chain Monte Carlo routine to invert river profiles and marine terrace data for rock uplift history and stream power parameters

## Author: Sean F. Gallen, 2020. [![DOI](https://zenodo.org/badge/299756636.svg)](https://zenodo.org/badge/latestdoi/299756636)

Department of Geosciences, Colorado State University

DOI: 10.5281/zenodo.4081849

# Introduction
This repository contains a number of MATLAB scripts and functions to run and plot the results of a Bayesian Markov chain Monte Carlo (MCMC) inversion (using the Metropolis-Hastings algorithm) of longitudinal river profiles and marine terrace-derived rock uplift rates developed by and presented in Gallen and Fernández-Blanco, 2021, "A New Data-driven Bayesian Inversion of Fluvial Topography Clarifies the Tectonic History of the Corinth Rift and Reveals a Channel Steepness Threshold", JGR-Earth Surface [https://doi.org/10.1029/2020JF005651].

The forward models assume a two-stage block-to-flexural (1D) uplift history, consistent with the evolution of the southern footwall of the Corinth Rift in Greece. The initial uplift is assumed to be spatially uniform and transitions to an uplift pattern consistent with flexure of broken elastic plate (footwall flexural uplift). The code uses the stream power incision model to evolve the river profiles through time. The stream power model states that:

    E = KA^mS^n

Where E is incision rate, A is upstream drainage area, S is local channel slope, K is a constant representing erodibility, and m and n are positive constants describe characteristics of hydrology, hydraulic scaling, and incision process (e.g. Whipple and Tucker, 1999).

Using marine terrace data from de Gelder et al. (2019) and river profiles extracted from a 20 m digital elevation model, the MCMC inverts for 7 unknown parameters: (1) the initial block uplift rate, (2) the final uplift rate at the broken plate segment (the fault trace), (3) the effective elastic thickness of the lithosphere (a proxy for lithospheric rigidity), (4) the timing of the transition from block to flexed uplift, (5) the erodibility constant, and (6) the m and (7) the n exponents in the stream power model. 
Gallen and Fernández-Blanco (2021) provide a comprehensive description of the forward and inverse models and the many assumptions underlying this approach. The scripts and functions contained within this repository allow the user to reproduce the general results of the inversion presented in the manuscript.

# Description of scripts and functions
The master script to run the inversion is block_2_flexed_profile_terrace_bayes_MCMC.m, which calls a series of functions to run the forward models (river_incision_forward_model.m, terrace_uplift_forward_model.m) and to help conduct the MCMC inversion (mod_loglikelihood.m, logprior_uniform.m, logproposal.m). For a nice introduction to inverse theory and Bayesian inversions, including MCMCs, the reader is referred to Aster et al. (2019) and for details on this specific model Gallen and Fernández-Blanco (2021). The script calls two data files, one related to marine terraces (corinth_terrace_elevation_age_de_Gelder_2019.xlsx) and the other related to river profiles (stream_data_east.mat). The river profile data is organized as a MATLAB structure based on Topotoolbox (Schwanghart and Scherler, 2014) and ChiProfiler (Gallen and Wegmann, 2017).

This repository also comes with a script to plot up the results after the inversion called “plot_basic_results.m”). This script will plot the chain paths, likelihood and acceptance rate history, and the best-fit model results compared to the observed data.

# Some practical notes
The inversion as presented by Gallen and Fernández-Blanco (2021) involved >3 million iterations that took more than three weeks to run on a single Intel Xeon E5-2680 v3 core. To put this another way, this model takes a long time to run in large part because of the 7 unknown parameters and the use of the relatively inefficient Metropolis-Hastings (MH) algorithm. I used the MH algorithm because it is easier to code from scratch relative to other fancier MCMC algorithms. 

Realizing how long the model takes to run, I set the number of total iterations in the main script to be 1.1e5 (1e4 for the burn in and 1e5 post-burn in). The script starts the model near the MAP solution, so even with the abbreviated number of iterations, the user can still get a good idea of how the model works. However, please note that the model still sometimes picks a "bad" initial starting position and struggles to start searching the parameter space efficiently with 1.1e5 iterations; therefore, it is good to run a couple of models to make sure things run smoothly.

Also note that, while MCMCs are great for their flexibility in application to various problems, they are difficult to generalize and thus require a significant amount of coding experience to be tailored for each application. This code is not sufficiently general to invert river profiles everywhere, rather is specifically designed for our application in Corinth. If you are interested in applying similar approaches but do not have the experience needed to revise the existing code, I am happy to collaborate and help you achieve your goals, so do not hesitate to email me.

# References:

Aster, Richard C., Brian Borchers, and Clifford H. Thurber. Parameter estimation and inverse problems. Elsevier, 2019.

Gallen and Fernández-Blanco, 2021, "A New Data-driven Bayesian Inversion of Fluvial Topography Clarifies the Tectonic History of the Corinth Rift and Reveals a Channel Steepness Threshold", JGR-Earth Surface. https://doi.org/10.1029/2020JF005651

Gallen, S. F., & Wegmann, K. W. (2017). River profile response to normal fault growth and linkage: An example from the Hellenic forearc of south-central Crete, Greece. Earth Surface Dynamics, 5(1). https://doi.org/10.5194/esurf-5-161-2017

de Gelder, G., Fernández-Blanco, D., Melnick, D., Duclaux, G., Bell, R. E., Jara-Muñoz, J., et al. (2019). Lithospheric flexure and rheology determined by climate cycle markers in the Corinth Rift. Scientific Reports, 9(1), 4260. https://doi.org/10.1038/s41598-018-36377-1

Schwanghart, W., & Scherler, D. (2014). Short Communication: TopoToolbox 2 – MATLAB-based software for topographic analysis and modeling in Earth surface sciences. Earth Surf. Dynam., 2(1), 1–7. https://doi.org/10.5194/esurf-2-1-2014

Whipple, K. X., & Tucker, G. E. (1999). Dynamics of the stream-power river incision model: Implications for height limits of mountain ranges, landscape response timescales, and research needs. Journal of Geophysical Research: Solid Earth, 104(B8), 17661–17674. https://doi.org/10.1029/1999JB900120
