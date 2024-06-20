# Spatial Confounding in Joint Species Distribution Models

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11078613.svg)](https://doi.org/10.5281/zenodo.11078613)

<!-- badges: end -->

This repository code contain template code associated with the manuscript "Spatial Confounding in Joint Species Distribution Models" which is currently in review.

# Getting started

The following directories are currently available in this repository:

-   `simulation`, which contains the following folders/documents for reproducing the simulation results in the manuscript:
    -   A `simdatfnv0.R` script that contains a function for simulating spatial multivariate abundance data from a joint species distribution model assuming some observed and some unobserved covariates;
    -   `simulationexample_coarsescaleX.pdf` and `simulationexample_finescaleX.pdf`, which are the plots available in the manuscript demonstrating a coarse or fine-scaled observed covariate, along with fine or coarse-scaled spatial errors forming the residual component of the unobserved covariates.
    -   `binary_N400` folder, which contains two `R` scripts and two corresponding `cpp` scripts for reproducing the simulations involving binary responses, $N = 400$ spatial locations from the Great Britain butterfly data, and either a coarse or fine-scaled observed covariate. The models compared include independent factor analysis, and spatial factor analysis/restricted spatial factor analysis (SFA/RSFA);
    -   `negin_N400` folder, which is set up analogously to the `binary_N400` folder except involving negative binomial count responses instead;
    -   `tweedie_N400` folder, which is set up analogously to the `binary_N400` folder except involving non-negative continuous responses from the Tweedie distribution instead;
    -   A `analysis.R` script for gathering and visualizing the simulation results, after the `R` scripts in the above three folders have been run;
    -   A subfolder `plots` containing some of the plots available in the manuscript.
-   `application_butterflies`, which contains template scripts to for reproducing the case study results for the Great Britain butterfly data in the manuscript. Also contained are:
    -   Two `cpp` scripts implementing independent factor analysis and spatial factor analysis/restricted spatial factor analysis (SFA/RSFA); these scripts are subsequently passed into [TMB](https://cran.r-project.org/web/packages/TMB/index.html) for model fitting;
    -   A subfolder `plots` containing some of the plots available in the manuscript.
-   `application_SOCPR`, which contains template scripts to for reproducing the case study results for the Southern Ocean Continuous Plankton Recorder survey data in the manuscript. Also contained are:
    -   Two `cpp` scripts implementing independent factor analysis and spatial factor analysis/restricted spatial factor analysis (SFA/RSFA); these scripts are subsequently passed into [TMB](https://cran.r-project.org/web/packages/TMB/index.html) for model fitting;
    -   A `socpr_reduced.RData` file that contains contains synthetic version of the data;
    -   A subfolder `plots` containing some of the plots available in the manuscript.

The simulation results for $N = 1000$ can be produced by straightforwardly modifying the `R` scripts in the `simulation` folder.

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated.
3.  Required data files etc...

For all other issues, please contact the corresponding author at [francis.hui\@anu.edu.au](mailto:francis.hui@anu.edu.au)
