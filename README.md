# Spatial Confounding in Joint Species Distribution Models

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10460492.svg)](https://doi.org/10.5281/zenodo.10460492)


<!-- badges: end -->

This repository code contain template code associated with the manuscript "Spatial Confounding in Joint Species Distribution Models" which is soon to be/currently in review.


# Getting started

The following directories are currently available in this repository:

-   `models`, which contains `cpp` scripts implementing independent factor analysis and spatial factor analysis/restricted spatial factor analysis (SFA/RSFA). These scripts are subsequently passed into [TMB](https://cran.r-project.org/web/packages/TMB/index.html) for model fitting;

-   `application_butterflies`, which contains template scripts to for reproducing the case study results for the Great Britain butterfly data in the manuscript. Also contained in this directories is a subfolder `plots` which contains some of the plots available in the manuscript.


# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [francis.hui\@anu.edu.au](mailto:francis.hui@anu.edu.au)
