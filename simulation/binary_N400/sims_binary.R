#' ---
#' title: Simulations for spatial confounding in joint species distribution models -- Binary responses
#' author: xxx
#' date: Code started April 2024
#' ---


##-----------------------
#' # Load in packages and data
##-----------------------
rm(list = ls())
library(tidyverse)
library(geoR)
library(autoFRK)
library(mvtnorm)
library(TMB) 
library(foreach) 
library(doParallel) 
library(Matrix) 
library(INLA) 
library(vegan) 
registerDoParallel(cores = detectCores()-4)
here::i_am("simulation/binary_N400/sims_binary.R")
library(here)
source(here("simulation", "simdatfnv0.R"))


longlat <- read.csv(file = here("application_butterflies", "Butterfly_LatLon_scaled.csv"))
cov_dat <- longlat
str(cov_dat)
rm(longlat)
#' Note latitude and longitudes have been scaled and so do not make sense here; so treat them as Euclidean coordinates/distance metric, as per https://doi.org/10.1111/2041-210X.13106


set.seed(042024)
num_sites <- 400
cov_dat <- cov_dat %>%
    slice_sample(n = num_sites)


num_spp <- 10
set.seed(042024)
X_coefficients <- rmvnorm(num_spp, mean = c(0,1), sigma = diag(x = c(0.1,0.1^2)))
Z_coefficients <- rmvnorm(num_spp, mean = c(1,1), sigma = diag(x = c(0.1^2,0.1^2)))
big_basis_fn <- mrts(knot = cov_dat %>% select(x:y), k = 20)


##-----------------------
#' # Construct simulation function and run!
##-----------------------
#' @title Simulates multivariate abundance data from a joint species distribution model (with 10 species and 400 sites obtained from the Great Britain butterfly data application) involving up to one observed covariates and two unobserved covariates that are both spatially correlated with the observed covariate. Then fits independent factor analysis along with spatial factor analysis/restricted spatial factor analysis (SFA/RSFA), and compared results.
#'
#' @param NAI The seed for simulating multivariate abundance data.
#' @param scale_X A scalar that effectively controls the degree of correlation between the observed covariate and the two unobserved covariates. In particular, the latter are constructed as \code{z = -scale_X*obs_X + eps}, where \code{obs_X} denotes the simulated observed covariate, and \code{eps} denotes an additional spatial field that is generated independently of (and in this simulation, orthogonally to) the observed covariate.
#' 
#' @return A .RData file is saved containing the following elements:
#' \item{X_coefficients_dat}{A data frame containing the "true" unrestricted species-specific coefficients for the dataset, along with the estimates from the three factor analytic models fitted.}
#' \item{X_coefficients_ci_dat}{A data frame containing coverage and widths of the 95\% confidence intervals for each "true" unrestricted species-specific coefficients, and from the three factor analytic models fitted.}
#' \item{LV_results_dat}{A list containing the Procrustes error between the true residual ordination (i.e., portion of the two missing covariates not attributable to the observed covariates) and the two predicted latent variables from the three factor analytic models fitted.}
#' \item{LVcor_results_dat}{A list containing the Pearson correlations between the two predicted latent variables from the three factor analytic models fitted, and the observed covariate}


simfn <- function(NAI, scale_X = 0) {
    
    ##-------------------------
    #' ## Generate multivariate abundance data    
    ##-------------------------
    set.seed(NAI)
    
    #' ### Make observed and unobserved covariates
    X1 <- (big_basis_fn[,c(2:5)] %*% rnorm(4, sd = 0.25))
    cov_dat$climate <- 0.5*X1 + rnorm(num_sites, sd = 0.05)

    eps1 <- big_basis_fn[,c(8,12,14)] %*% rnorm(3, sd = 0.1)
    eps2 <- big_basis_fn[,c(9,13,16)] %*% rnorm(3, sd = 0.1)

    z1 <- (-scale_X*X1 - eps1) 
    z2 <- (-scale_X*X1 + eps2) 
    
    #' ### Construct "true" unrestricted species-specific coefficients for the dataset, as well as the true residual ordination.
    full_X_coefficients <- X_coefficients + t(lm(cbind(z1,z2) ~ climate, data = cov_dat)$coef %*% t(Z_coefficients))
    res_ZonX <- lm(cbind(z1,z2) ~ climate, data = cov_dat) %>% residuals 


    simdat <- simdat_jsdm_simple(family = "binomial",
                                        x1 = cov_dat$climate, X_coefficients = X_coefficients,
                                        z1 = z1, z2 = z2, Z_coefficients = Z_coefficients,
                                        seed = NAI)

    ##-------------------------
    #' ## Prepare objects and fit independent GLLVM 
    #' Fitting is done using Template Model Builder (TMB). 
    #' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood.
    ##-------------------------
    X <- model.matrix(~ climate, data = cov_dat) 
    num_lv <- 2
    
    filename <- "TMB_indFA.cpp"
    modelname <- strsplit(filename, "\\.")[[1]][1]
    compile(filename)
    dyn.load(dynlib(modelname))
    
    
    tidbits_data <- list(y = simdat %>% as.matrix, 
                         X = X, 
                         num_lv = num_lv)
    
    blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
    tidbits_parameters <- list(betas = matrix(0, num_spp, ncol(X)),
                               loadings = numeric(sum(lower.tri(blank_loadings, diag = TRUE))),
                               lvs = matrix(0, nrow = nrow(X), ncol = num_lv))
    rm(blank_loadings)
    
    tidbits_random <- c("lvs")
    
    objs <- MakeADFun(data = tidbits_data, 
                      parameters = tidbits_parameters, 
                      random = tidbits_random, 
                      DLL = modelname, 
                      hessian = FALSE, 
                      silent = TRUE)
    
    blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
    rownames(blank_loadings) <- colnames(tidbits_data$y)
    colnames(blank_loadings) <- paste0("LV", 1:num_lv)
    blank_loadings[lower.tri(blank_loadings,diag=TRUE)] <- 1:sum(lower.tri(blank_loadings,diag=TRUE))		
    upperlim <- rep(Inf, length(objs$par))
    lowerlim <- rep(-Inf, length(objs$par))
    lowerlim[grep("loadings",names(objs$par))][diag(blank_loadings)] <- 1e-3 ## Constraints on loading matrix    
    rm(blank_loadings)
    
    
    #' ### Fit and gather results
    fit_lvm <- nlminb(start = objs$par, 
                      objective = objs$fn, 
                      gradient = objs$gr, 
                      lower = lowerlim, 
                      upper = upperlim,
                      control = list(iter.max = 2000, eval.max = 2000, trace = 1))
    
    indlvm_results <- sdreport(objs) 
    pt_estimates <- as.list(indlvm_results, what = "Estimate", report = TRUE)
    all_results <- summary(indlvm_results, select = "report", p.value = TRUE)
    
    
    ci_alpha <- 0.05
    betaind_results <- data.frame(all_results[grep("betas$", rownames(all_results)),])
    betaind_results$lower <- betaind_results$Estimate - qnorm(1-ci_alpha/2) * betaind_results$Std
    betaind_results$upper <- betaind_results$Estimate + qnorm(1-ci_alpha/2) * betaind_results$Std
    rownames(betaind_results) <- paste0(rep(colnames(tidbits_data$y), ncol(X)), ":", rep(colnames(X), each = ncol(tidbits_data$y)))
    colnames(betaind_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 
    
    lvind_results <- pt_estimates$lvs
    loadingind_results <- pt_estimates$loadings_mat
    
    eta_ind <- list(Xbeta = tcrossprod(X, pt_estimates$betas), 
                    residual = tcrossprod(pt_estimates$lvs, pt_estimates$loadings_mat), 
                    residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
                    rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor)    
    
    dyn.unload(paste0(modelname,".so"))
    #file.remove(filename, paste0(modelname,".so"), paste0(modelname,".o"))
    gc()
    rm(all_results, pt_estimates, fit_lvm, indlvm_results)
    
    
    
    ##-------------------------
    #' ## Fit unrestricted spatial factor analytic (SFA) and restricted spatial factor analystic (RSFA) model.
    #' Fitting is done using Template Model Builder (TMB). Since the spatial domain is not very large (Grampians national park) and given the scaling that has taken place prior to accessing the data, then Euclidean coordinates are assumed. 
    #' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood, where the Matern covariance is approximated using a SPDE approach.
    #' Standard errors obtained using a generalized Delta method based on the joint covariance of the random effects and parameter estimates
    ##-------------------------
    filename <- "TMBINLA_SFA.cpp"
    modelname <- strsplit(filename, "\\.")[[1]][1]
    compile(filename)
    dyn.load(dynlib(modelname))
    
    
    #' ### Create SPDE mesh
    mesh <- inla.mesh.create(cov_dat %>% dplyr::select(x:y), plot.delay = NULL, refine = FALSE)
    spde <- inla.spde2.matern(mesh, alpha = 2)
    
    tidbits_data <- list(y = simdat %>% as.matrix, 
                         X = X, 
                         num_lv = num_lv,
                         OLSmatrix_transpose = X %*% solve(crossprod(X)), 
                         residualprojection = diag(nrow = nrow(X)) - X %*% tcrossprod(solve(crossprod(X)), X),
                         n_i = mesh$n,
                         meshidxloc=mesh$idx$loc - 1,
                         G0 = spde$param.inla$M0, 
                         G1 = spde$param.inla$M1, 
                         G2 = spde$param.inla$M2
                         )
    
    blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
    tidbits_parameters <- list(betas = matrix(0, nrow = num_spp, ncol = ncol(X)),
                               loadings = numeric(sum(lower.tri(blank_loadings, diag = TRUE))),
                               lvs = matrix(0, nrow = mesh$n, ncol = num_lv),
                               log_kappa = rep(0, num_lv))
    rm(blank_loadings, mesh, spde)
    
    tidbits_random <- c("lvs")
    
    
    objs <- MakeADFun(data = tidbits_data, 
                      parameters = tidbits_parameters, 
                      random = tidbits_random, 
                      DLL = modelname, 
                      hessian = FALSE, 
                      silent = TRUE)
    
    blank_loadings <- matrix(0, nrow = ncol(tidbits_data$y), ncol = num_lv)
    rownames(blank_loadings) <- colnames(tidbits_data$y)
    colnames(blank_loadings) <- paste0("LV", 1:num_lv)
    blank_loadings[lower.tri(blank_loadings,diag=TRUE)] <- 1:sum(lower.tri(blank_loadings,diag=TRUE))		
    upperlim <- rep(Inf, length(objs$par))
    lowerlim <- rep(-Inf, length(objs$par))
    lowerlim[grep("loadings",names(objs$par))][diag(blank_loadings)] <- 1e-3 ## Constraints on loading matrix    
    rm(blank_loadings)
    
    
    #' ### Fit and gather results
    fit_splvm <- nlminb(start = objs$par, 
                        objective = objs$fn, 
                        gradient = objs$gr, 
                        lower = lowerlim, 
                        upper = upperlim, 
                        control = list(iter.max = 2000, eval.max = 2000, trace = 1))
    
    
    splvm_results <- sdreport(objs)
    pt_estimates <- as.list(splvm_results, what = "Estimate", report = TRUE)
    all_results <- summary(splvm_results, select = "report", p.value = TRUE)
    
    
    ci_alpha <- 0.05
    betaSFA_results <- data.frame(all_results[grep("betas$", rownames(all_results)),])
    betaSFA_results$lower <- betaSFA_results$Estimate - qnorm(1-ci_alpha/2) * betaSFA_results$Std
    betaSFA_results$upper <- betaSFA_results$Estimate + qnorm(1-ci_alpha/2) * betaSFA_results$Std
    rownames(betaSFA_results) <- paste0(rep(colnames(tidbits_data$y), ncol(X)), ":", rep(colnames(X), each = ncol(tidbits_data$y)))
    colnames(betaSFA_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 
    
    betaRSFA_results <- data.frame(all_results[grep("betas_rsr", rownames(all_results)),])
    betaRSFA_results$lower <- betaRSFA_results$Estimate - qnorm(1-ci_alpha/2) * betaRSFA_results$Std
    betaRSFA_results$upper <- betaRSFA_results$Estimate + qnorm(1-ci_alpha/2) * betaRSFA_results$Std
    rownames(betaRSFA_results) <- paste0(rep(colnames(tidbits_data$y), ncol(X)), ":", rep(colnames(X), each = ncol(tidbits_data$y)))
    colnames(betaRSFA_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 
    
    lvSFA_results <- pt_estimates$lvs_units
    lvRSFA_results <- pt_estimates$lvs_rsr
    loadingSFA_results <- loadingRSFA_results <- pt_estimates$loadings_mat
    
    eta_SFA <- list(Xbeta = tcrossprod(X, pt_estimates$betas), 
                    residual = tcrossprod(pt_estimates$lvs_units, pt_estimates$loadings_mat),
                    residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
                    rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor)
    
    eta_RSFA <- list(Xbeta = tcrossprod(X, pt_estimates$betas_rsr), 
                     residual = tcrossprod(pt_estimates$lvs_rsr, pt_estimates$loadings_mat),
                     residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
                     rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor)
    
    
    
    dyn.unload(paste0(modelname,".so"))
    #file.remove(filename, paste0(modelname,".so"), paste0(modelname,".o"))
    
    rm(objs, tidbits_parameters, ci_alpha, filename, modelname, lowerlim, upperlim, tidbits_random, splvm_results, fit_splvm, all_results, pt_estimates)
    gc()
    
    
    ##-------------------------
    #' ## Compare results from the three models
    ##-------------------------
    #' ### Point performance of covariate effects, comparing against the "true" unrestricted species-specific coefficients for the dataset
    X_coefficients_dat <- data.frame(true = full_X_coefficients %>% unlist %>% as.vector(),
                                     indFA = betaind_results$Estimate,
                                     SFA = betaSFA_results$Estimate,
                                     RSFA = betaRSFA_results$Estimate) %>% 
        slice(-c(1:num_spp))
    
    #' ### Coverage probability of covariate effects
    X_coefficients_ci_dat <- rbind(
        bind_cols(betaind_results, true = full_X_coefficients %>% unlist %>% as.vector) %>% 
            mutate(coverage_probability = ifelse(lower < true & upper > true, 1, 0),
                   interval_length = upper - lower) %>% 
            slice(-(1:num_spp)) %>% 
            select(coverage_probability, interval_length),
        bind_cols(betaSFA_results, true = full_X_coefficients %>% unlist %>% as.vector) %>% 
            mutate(coverage_probability = ifelse(lower < true & upper > true, 1, 0),
                   interval_length = upper - lower) %>% 
            slice(-(1:num_spp)) %>% 
            select(coverage_probability, interval_length),
        bind_cols(betaRSFA_results, true = full_X_coefficients %>% unlist %>% as.vector) %>% 
            mutate(coverage_probability = ifelse(lower < true & upper > true, 1, 0),
                   interval_length = upper - lower) %>% 
            slice(-(1:num_spp)) %>% 
            select(coverage_probability, interval_length)) %>% 
        mutate(model = rep(c("Independent", "SFA", "RSFA"), each = (ncol(X)-1)*num_spp) %>% fct_inorder) 
    
    
    #' ### Recovery of unobserved component of unobserved covariates i.e., a scaled version of cbind(eps1,eps2)
    LV_results_dat <- list(Independent = procrustes(cbind(eps1, eps2), lvind_results),
                                 SFA = procrustes(cbind(eps1, eps2), lvSFA_results),
                                 RSFA = procrustes(cbind(eps1, eps2), lvRSFA_results))

    
    LVcor_results_dat <- list(Independent = cor(lvind_results %>% as.data.frame() %>% rename(LV1 = 1, LV2 = 2), X[,-1]),
                              SFA = cor(lvSFA_results %>% as.data.frame() %>% rename(LV1 = 1, LV2 = 2), X[,-1]),
                              RSFA = cor(lvRSFA_results %>% as.data.frame() %>% rename(LV1 = 1, LV2 = 2), X[,-1]))
    
    save(X_coefficients_dat, 
         X_coefficients_ci_dat,
         LV_results_dat, 
         LVcor_results_dat,
         file = paste0("simulation_binary_N", num_sites, "_scaleX", scale_X, "_dataset", NAI, ".RData"))
    }


num_datasets <- 200
foreach(k1 = 1:num_datasets) %dopar% simfn(NAI = k1, scale_X = 0)
foreach(k1 = 1:num_datasets) %dopar% simfn(NAI = k1, scale_X = 0.5)
foreach(k1 = 1:num_datasets) %dopar% simfn(NAI = k1, scale_X = 1)
foreach(k1 = 1:num_datasets) %dopar% simfn(NAI = k1, scale_X = 2)




##------------------------------
#' # Assessment of results
##------------------------------
num_datasets <- 200
scale_X_seq <- c(0,0.5,1,2)

all_diff_X_coefficients <- all_coverage_X_coefficients <- all_length_X_coefficients <- array(NA, 
                                                                                             dim = c(num_datasets, length(scale_X_seq), num_spp, 3),
                                                                                             dimnames = list(dataset = 1:num_datasets, scale_X = scale_X_seq, covariate = 1:(num_spp), method = c("Independent", "SFA", "RSFA"))) 
all_LVerror <- array(NA, 
                     dim = c(num_datasets, length(scale_X_seq), 3), 
                     dimnames = list(dataset = 1:num_datasets, scale_X = scale_X_seq, method = c("Independent", "SFA", "RSFA"))) 
all_LV_X_correlation <- array(NA, 
                              dim = c(num_datasets, length(scale_X_seq), 2, 3),
                              dimnames = list(dataset = 1:num_datasets, scale_X = scale_X_seq, LV = 1:2, method = c("Independent", "SFA", "RSFA"))) 


for(k0 in 1:length(scale_X_seq)) { for(k1 in 1:num_datasets) {
    make_filename <- here("simulation", "binary_N400", paste0("simulation_binary_N", num_sites, "_scaleX", scale_X_seq[k0], "_dataset", k1, ".RData"))
    
    if(!file.exists(make_filename))
        next;
    
    if(file.exists(make_filename))
        load(file = make_filename)
    
    all_diff_X_coefficients[k1,k0,,] <- X_coefficients_dat %>% mutate(across(everything(), ~(.x - true))) %>% select(-true) %>% as.matrix
    all_coverage_X_coefficients[k1,k0,,] <- matrix(X_coefficients_ci_dat$coverage_probability, ncol = 3, byrow = FALSE) 
    all_length_X_coefficients[k1,k0,,] <- matrix(X_coefficients_ci_dat$interval_length, ncol = 3, byrow = FALSE) 
    all_LVerror[k1,k0,] <- c(LV_results_dat$Independent$ss, LV_results_dat$SFA$ss, LV_results_dat$RSFA$ss)
    all_LV_X_correlation[k1,k0,,] <- cbind(LVcor_results_dat$Independent %>% abs,
                                           LVcor_results_dat$SFA %>% abs,
                                           LVcor_results_dat$RSFA %>% abs)
    
    } }



apply(all_diff_X_coefficients, c(2,3,4), function(x) mean(x, na.rm = TRUE))
apply(all_diff_X_coefficients, c(2,3,4), function(x) abs(mean(x, na.rm = TRUE)))
apply(all_diff_X_coefficients^2, c(2,3,4), function(x) sqrt(mean(x, na.rm = TRUE)))
apply(all_coverage_X_coefficients, c(2,3,4), mean, na.rm = TRUE)
apply(all_length_X_coefficients, c(2,3,4), mean, na.rm = TRUE)

apply(all_LVerror, c(2,3), mean, na.rm = TRUE)
apply(all_LV_X_correlation, c(2,3,4), mean, na.rm = TRUE)




##----------------------
sessioninfo::session_info()
##----------------------
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 (2021-11-01)
# os       Linux Mint 21.1
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2024-04-08
# rstudio  2023.12.0+369 Ocean Storm (desktop)
# pandoc   NA
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package           * version  date (UTC) lib source
# beeswarm            0.4.0    2021-06-01 [1] CRAN (R 4.1.2)
# cli                 3.6.1    2023-03-23 [1] CRAN (R 4.1.2)
# cluster             2.1.2    2021-04-17 [4] CRAN (R 4.1.1)
# codetools           0.2-18   2020-11-04 [4] CRAN (R 4.0.3)
# colorspace          2.1-0    2023-01-23 [1] CRAN (R 4.1.2)
# doParallel        * 1.0.17   2022-02-07 [1] CRAN (R 4.1.2)
# dplyr             * 1.1.2    2023-04-20 [1] CRAN (R 4.1.2)
# ellipsis            0.3.2    2021-04-29 [1] CRAN (R 4.1.2)
# fansi               1.0.4    2023-01-22 [1] CRAN (R 4.1.2)
# forcats           * 1.0.0    2023-01-29 [1] CRAN (R 4.1.2)
# foreach           * 1.5.2    2022-02-02 [1] CRAN (R 4.1.2)
# generics            0.1.3    2022-07-05 [1] CRAN (R 4.1.2)
# GGally              2.1.2    2021-06-21 [1] CRAN (R 4.1.2)
# ggbeeswarm          0.7.2    2023-04-29 [1] CRAN (R 4.1.2)
# ggplot2           * 3.4.2    2023-04-03 [1] CRAN (R 4.1.2)
# glue                1.6.2    2022-02-24 [1] CRAN (R 4.1.2)
# gtable              0.3.3    2023-03-21 [1] CRAN (R 4.1.2)
# here              * 1.0.1    2020-12-13 [1] CRAN (R 4.1.2)
# hms                 1.1.2    2022-08-19 [1] CRAN (R 4.1.2)
# INLA              * 22.12.16 2022-12-23 [1] local
# iterators         * 1.0.14   2022-02-05 [1] CRAN (R 4.1.2)
# lattice           * 0.20-45  2021-09-22 [4] CRAN (R 4.1.1)
# lifecycle           1.0.3    2022-10-07 [1] CRAN (R 4.1.2)
# lubridate         * 1.9.2    2023-02-10 [1] CRAN (R 4.1.2)
# magrittr            2.0.3    2022-03-30 [1] CRAN (R 4.1.2)
# MASS                7.3-55   2022-01-13 [4] CRAN (R 4.1.2)
# Matrix            * 1.6-4    2023-11-30 [1] CRAN (R 4.1.2)
# mgcv                1.9-0    2023-07-11 [1] CRAN (R 4.1.2)
# munsell             0.5.0    2018-06-12 [1] CRAN (R 4.1.2)
# nlme                3.1-155  2022-01-13 [4] CRAN (R 4.1.2)
# permute           * 0.9-7    2022-01-27 [1] CRAN (R 4.1.2)
# pillar              1.9.0    2023-03-22 [1] CRAN (R 4.1.2)
# pkgconfig           2.0.3    2019-09-22 [1] CRAN (R 4.1.2)
# plyr                1.8.8    2022-11-11 [1] CRAN (R 4.1.2)
# purrr             * 1.0.1    2023-01-10 [1] CRAN (R 4.1.2)
# R6                  2.5.1    2021-08-19 [1] CRAN (R 4.1.2)
# RandomFields      * 3.3.14   2022-01-18 [1] CRAN (R 4.1.2)
# RandomFieldsUtils * 1.2.5    2022-04-19 [1] CRAN (R 4.1.2)
# RColorBrewer        1.1-3    2022-04-03 [1] CRAN (R 4.1.2)
# Rcpp                1.0.10   2023-01-22 [1] CRAN (R 4.1.2)
# readr             * 2.1.4    2023-02-10 [1] CRAN (R 4.1.2)
# reshape             0.8.9    2022-04-12 [1] CRAN (R 4.1.2)
# rlang               1.1.1    2023-04-28 [1] CRAN (R 4.1.2)
# rprojroot           2.0.3    2022-04-02 [1] CRAN (R 4.1.2)
# rstudioapi          0.14     2022-08-22 [1] CRAN (R 4.1.2)
# scales              1.2.1    2022-08-20 [1] CRAN (R 4.1.2)
# sessioninfo         1.2.2    2021-12-06 [1] CRAN (R 4.1.2)
# sp                * 1.6-0    2023-01-19 [1] CRAN (R 4.1.2)
# stringi             1.7.12   2023-01-11 [1] CRAN (R 4.1.2)
# stringr           * 1.5.0    2022-12-02 [1] CRAN (R 4.1.2)
# tibble            * 3.2.1    2023-03-20 [1] CRAN (R 4.1.2)
# tidyr             * 1.3.0    2023-01-24 [1] CRAN (R 4.1.2)
# tidyselect          1.2.0    2022-10-10 [1] CRAN (R 4.1.2)
# tidyverse         * 2.0.0    2023-02-22 [1] CRAN (R 4.1.2)
# timechange          0.2.0    2023-01-11 [1] CRAN (R 4.1.2)
# TMB               * 1.9.10   2023-12-12 [1] CRAN (R 4.1.2)
# tweedie           * 2.3.5    2022-08-17 [1] CRAN (R 4.1.2)
# tzdb                0.3.0    2022-03-28 [1] CRAN (R 4.1.2)
# utf8                1.2.3    2023-01-31 [1] CRAN (R 4.1.2)
# vctrs               0.6.2    2023-04-19 [1] CRAN (R 4.1.2)
# vegan             * 2.6-4    2022-10-11 [1] CRAN (R 4.1.2)
# vipor               0.4.5    2017-03-22 [1] CRAN (R 4.1.2)
# withr               2.5.0    2022-03-03 [1] CRAN (R 4.1.2)
# 
# [1] /home/fkch/R/x86_64-pc-linux-gnu-library/4.1
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
