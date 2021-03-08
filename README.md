# OPFOR
Code to implement Ordinal Probit Functional Outcome Regression as described by Meyer et al. (2021).

sample_script.m contains a sample script based off of a single simulated dataset. Sample code to implement all six models from the manuscript is in this file.

The main function, opfor contained in the file opfor.m, is a wrapper function that makes a call to the function opfm (in opfm.m). The function opfor was designed to run using defaults for the user's choice of O-Spline, B-Spline, or wavelet. Defaults are described in the manuscript.

N.B.: To run O-Spline models, the R file ``osullivan.R'' must be placed in the main MATLAB folder.

The Scripts folder contains scripts to run the simulations from manuscript including the naive simulation based on cumulative-link mixed-effects models (CLMM).

N.B.: To run the simulations for the CLMMs, the R file ``mixed.R'' must be placed in the main MATLAB folder.
