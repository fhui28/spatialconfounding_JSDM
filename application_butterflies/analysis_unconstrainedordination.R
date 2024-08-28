#' ---
#' title: Analysis of Butterflies data from Ovaskainen et al., (2016) using factor analytic models to explore the presence and impacts of spatial confounding on model-based unconstrained ordination only
#' author: Francis KC Hui
#' date: Code started December 2023
#' ---

#' # In preparation for takeoff
#' 
#' The data are sourced from [Wilkinson et al., (2019)](https://doi.org/10.1111/2041-210X.13106) and are publicly available from (https://zenodo.org/records/1452066) under `Datasets/Butterflies` with three files"
#' - Butterfly_PA.csv
#' - Butterfly_LatLon_scaled.csv
#' - Butterfly_Cov.csv
#' 
#' The script below assumes these datasets have been loaded into the working directory

##-----------------------
#' # Load in packages and explore data
##-----------------------
rm(list = ls())
library(tidyverse)
library(patchwork) 
library(TMB) 
library(INLA) 
library(DHARMa) 
# library(mgcv) 
library(foreach) 
library(doParallel) 
library(corrplot) 
#library(kableExtra) 
library(Matrix) 
registerDoParallel(cores = 6)
here::i_am("application_butterflies/analysis_unconstrainedordination.R")
library(here)


resp_dat <- read.csv(file = here("application_butterflies", "Butterfly_PA.csv"))
str(resp_dat)
colSums(resp_dat)
resp_dat <- resp_dat[,colSums(resp_dat) > 150] # Subset to most prevalence species. This step removes 10 of 55 species
str(resp_dat)
colSums(resp_dat)


cov_dat <- read.csv(file = here("application_butterflies", "Butterfly_Cov.csv"))
longlat <- read.csv(file = here("application_butterflies", "Butterfly_LatLon_scaled.csv"))
cov_dat <- cbind(longlat, cov_dat[,-1]) # Remove intercept
summary(cov_dat)
cov_dat <- cov_dat %>% 
    mutate(blwood = scale(blwood) %>% as.vector) %>% 
    mutate(conwood = scale(cov_dat$conwood) %>% as.vector) # Note growing degree days (climate) has already been standardized
str(cov_dat) 

summary(cov_dat) # Note latitude and longitudes have been scaled and so do not make sense here; so treat them as Euclidean coordinates/distance metric, as per https://doi.org/10.1111/2041-210X.13106
boxplot(cov_dat)


#' ## Some simple spatial plots
ggplot(cov_dat %>% pivot_longer(-(x:y), names_to = "predictors"), aes(x = x, y = y, color = value)) +
    geom_point() +
    facet_wrap(. ~ predictors, nrow = 4, scales = "free") +
    scale_color_viridis_c() +
    theme_bw()

cov_datlong <- cov_dat %>% 
    pivot_longer(-(x:y), names_to = "predictors") %>% 
    filter(predictors != "chalk_limestone") %>% 
    mutate(predictors = fct_inorder(predictors)) 
levels(cov_datlong$predictors) <- c("GDD", "Broadleaved woodland", "Coniferous woodland")

ggplot(cov_datlong, aes(x = x, y = y, color = value)) +
    geom_point() +
    facet_wrap(. ~ predictors, nrow = 1, scales = "free") +
    scale_color_viridis_c() +
    labs(color = "Value") +
    theme_bw() #+
ggsave(file = here("application_butterflies", "plots", "covariates.pdf"), width = 8, height = 4)


ggplot(cbind(longlat, resp_dat)  %>% pivot_longer(-(x:y), names_to = "species"), aes(x = x, y = y, color = value %>% factor)) +
    geom_point() +
    facet_wrap(. ~ species, nrow = 4, scales = "free") +
    scale_color_viridis_d() +
    theme_bw() # Not really too useful given the number of species...



#' ## Is there evidence of residual spatial correlation?
num_spp <- resp_dat %>% 
    ncol
for(k0 in 1:num_spp) {
    fit_cw <- glm(response ~ climate + blwood + conwood, family = binomial(), 
                  data = data.frame(response = resp_dat[,k0], cov_dat))
    dotestSAC <- simulateResiduals(fit_cw) %>% testSpatialAutocorrelation(., x = cov_dat$x, y = cov_dat$y, plot = FALSE)
    print(dotestSAC)
}
rm(fit_cw, dotestSAC)
#' All species exhibit evidence of residual spatial correlation, after accounting for the three covariates, although note the responses are binary!

dev.off()




##-------------------------
#' # Prepare objects and fit independent GLLVM 
#' Fitting is done using Template Model Builder (TMB) -- Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood.
##-------------------------
X <- model.matrix(~ 1, data = cov_dat) 
num_lv <- 2

filename <- "TMB_indFA.cpp"
modelname <- strsplit(filename, "\\.")[[1]][1]
compile(filename)
dyn.load(dynlib(modelname))


tidbits_data <- list(y = resp_dat %>% as.matrix, 
                     X = X, 
                     num_lv = num_lv)

blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
tidbits_parameters <- list(betas = matrix(0, nrow = num_spp, ncol = ncol(X)),
                           alpha = numeric(ncol(X)),
                           log_sds = numeric(ncol(X)),
                           loadings = numeric(sum(lower.tri(blank_loadings, diag = TRUE))),
                           lvs = matrix(0, nrow = nrow(X), ncol = num_lv))
rm(blank_loadings)

tidbits_random <- c("lvs", "betas")


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


#' ## Fit and gather results
fit_lvm <- nlminb(start = objs$par, 
                  objective = objs$fn, 
                  gradient = objs$gr, 
                  lower = lowerlim, 
                  upper = upperlim,
                  control = list(iter.max = 10000, eval.max = 10000, trace = 1))

indlvm_results <- sdreport(objs) 
pt_estimates <- as.list(indlvm_results, what = "Estimate", report = TRUE)
all_results <- summary(indlvm_results, select = "report", p.value = TRUE)


ci_alpha <- 0.05
betaind_results <- data.frame(all_results[grep("betas$", rownames(all_results)),])
betaind_results$lower <- betaind_results$Estimate - qnorm(1-ci_alpha/2) * betaind_results$Std
betaind_results$upper <- betaind_results$Estimate + qnorm(1-ci_alpha/2) * betaind_results$Std
rownames(betaind_results) <- paste0(rep(colnames(tidbits_data$y), ncol(X)), ":", rep(colnames(X), each = ncol(tidbits_data$y)))
colnames(betaind_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 

alphaind_results <- data.frame(all_results[grep("alpha$", rownames(all_results)),,drop=FALSE])
alphaind_results$lower <- alphaind_results$Estimate - qnorm(1-ci_alpha/2) * alphaind_results$Std
alphaind_results$upper <- alphaind_results$Estimate + qnorm(1-ci_alpha/2) * alphaind_results$Std
rownames(alphaind_results) <- colnames(X)
colnames(alphaind_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper")

lvind_results <- pt_estimates$lvs
loadingind_results <- pt_estimates$loadings_mat

eta_ind <- list(Xbeta = tcrossprod(X, pt_estimates$betas), 
                residual = tcrossprod(pt_estimates$lvs, pt_estimates$loadings_mat), 
                residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
                rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor)    

dyn.unload(paste0(modelname,".so"))
gc()
rm(all_results, pt_estimates, fit_lvm, indlvm_results)


save(betaind_results, 
     alphaind_results,
     eta_ind, 
     loadingind_results, 
     lvind_results, 
     file = here("application_butterflies", "unconstrainedordination_independentspatialfit.RData"))


##-------------------------
#' # Prepare objects and fit spatial factor analytic (SFA) model. 
#' Fitting is done using Template Model Builder (TMB). Since the spatial domain is not very large (Grampians national park) and given the scaling that has taken place prior to accessing the data, then Euclidean coordinates are assumed. 
#' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood, where the Matern covariance is approximated using a SPDE approach.
##-------------------------
X <- model.matrix(~ 1, data = cov_dat) 
num_lv <- 2

filename <- "TMB_SFA.cpp"
modelname <- strsplit(filename, "\\.")[[1]][1]
compile(filename)
dyn.load(dynlib(modelname))


# Create SPDE mesh
Loc <- cov_dat %>% dplyr::select(x:y)
mesh <- inla.mesh.create(Loc, plot.delay = NULL, refine = FALSE)
spde <- inla.spde2.matern(mesh, alpha = 2)


tidbits_data <- list(y = resp_dat %>% as.matrix, 
                     X = X, 
                     num_lv = num_lv,
                     OLSmatrix_transpose = X %*% solve(crossprod(X)), 
                     residualprojection = diag(nrow = nrow(X)) - X %*% tcrossprod(solve(crossprod(X)), X),
                     n_i = mesh$n,
                     meshidxloc=mesh$idx$loc - 1,
                     G0 = spde$param.inla$M0, 
                     G1 = spde$param.inla$M1, 
                     G2 = spde$param.inla$M2)

blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
tidbits_parameters <- list(betas = matrix(0, nrow = num_spp, ncol = ncol(X)),
                           alpha = numeric(ncol(X)),
                           log_sds = numeric(ncol(X)),
                           loadings = numeric(sum(lower.tri(blank_loadings, diag = TRUE))),
                           lvs = matrix(0, nrow = mesh$n, ncol = num_lv),
                           log_kappa = rep(0, num_lv))
rm(blank_loadings, mesh, spde)

tidbits_random <- c("lvs", "betas")


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


#' ## Fit and gather results
fit_splvm <- nlminb(start = objs$par, 
                    objective = objs$fn, 
                    gradient = objs$gr, 
                    lower = lowerlim, 
                    upper = upperlim, 
                    control = list(iter.max = 10000, eval.max = 10000, trace = 1))


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


alphaSFA_results <- data.frame(all_results[grep("alpha$", rownames(all_results)),])
alphaSFA_results$lower <- alphaSFA_results$Estimate - qnorm(1-ci_alpha/2) * alphaSFA_results$Std
alphaSFA_results$upper <- alphaSFA_results$Estimate + qnorm(1-ci_alpha/2) * alphaSFA_results$Std
rownames(alphaSFA_results) <- colnames(X)
colnames(alphaSFA_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper")

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

rm(Loc, objs, tidbits_parameters, ci_alpha, filename, modelname, lowerlim, upperlim, tidbits_random, splvm_results, fit_splvm, all_results, pt_estimates)
gc()

save(betaRSFA_results, betaSFA_results,
     alphaSFA_results,
     lvSFA_results, lvRSFA_results,
     eta_RSFA, eta_SFA,
     loadingSFA_results, loadingRSFA_results,
     file = here("application_butterflies", "unconstrainedordination_spatialfit.RData"))



##------------------------------
#' # Explore results for unconstrained ordination
##------------------------------
load(file = here("application_butterflies", "unconstrainedordination_independentspatialfit.RData"))
load(file = here("application_butterflies", "unconstrainedordination_spatialfit.RData"))

X <- model.matrix(~ 1, data = cov_dat) 
num_lv <- 2

tidbits_data <- list(y = resp_dat %>% as.matrix, 
                     X = X, 
                     num_lv = num_lv)


spp_names <- c("small.tortoiseshell", "orange.tip", "ringlet", "dgreen.fritillary",  
               "swashed.fritillary", "brown.argus", "nbrown.argus", 
               "pbordered.fritillary", "spbordered.fritillary", 
               "green.hairstreak", "holly.blue", "small.heath", 
               "large.heath", "clouded.yellow", "small.blue", "scotch.argus", 
               "dingy.skipper", "marsh.fritillary", "brimstone", 
               "doburgundy.fritillary", "grayling", "peacock", 
               "wall.brown", "small.copper", "adonis.blue", 
               "chill.blue", "meadow.brown", "marble.white", "purple.hairstreak", 
               "camberwell.beauty", "large.skipper", "speckled.wood", "large.white", 
               "gveined.white", "small.white", "sstudded.blue", "comma", 
               "common.blue", "grizzled.skipper", "gatekeeper", "brown.hairstreak", 
               "essex.skipper", "small.skipper", "red.admiral", "painted.lady") %>% 
    str_replace("\\.", " ") %>% 
    str_to_title()



#' ## Plots of point estimates for latent variables
results_dat <- rbind(lvind_results, lvSFA_results, lvRSFA_results) %>% #lvspplus_results
    as.data.frame()
colnames(results_dat)[1:num_lv] <- paste("Factor", 1:num_lv)
results_dat$model <- rep(c("independent", "SFA", "RSFA"), each = nrow(X)) %>% fct_inorder #"SFA+"
results_dat$x <- rep(cov_dat$x, 3)
results_dat$y <- rep(cov_dat$y, 3)
levels(results_dat$model)[1] <- "Independent"


ggplot(results_dat %>% rename(F1 = "Factor 1", F2 = "Factor 2"), aes(x = F1, y = F2)) +
    geom_point(color = "lightgrey") +
    facet_wrap(. ~ model) +
    theme_bw()
#ggsave(file = here("application_butterflies", "plots", "unconstrainedordination.pdf"), width = 8, height = 8)


ggplot(results_dat %>% pivot_longer("Factor 1":"Factor 2", names_to = "Factor"), aes(x = model, y = value)) +
    #geom_jitter(color = "lightgrey") +
    geom_boxplot(notch = TRUE) +
    ggbeeswarm::geom_quasirandom(color = "lightgrey", method = "pseudorandom", alpha = 0.2) +
    geom_hline(yintercept = 0, color = "darkblue", linetype = 2) +
    facet_wrap(. ~ Factor) +
    theme_bw()
ggsave(file = here("application_butterflies", "plots", "unconstrainedordination_boxplotLV.pdf"), width = 8, height = 8)


results_dat %>% 
    pivot_longer("Factor 1":"Factor 2", names_to = "Factor") %>%
    group_by(model, Factor) %>%
    summarise(mean = mean(value) %>% round(3), median = median(value) %>% round(3))


ggplot(results_dat %>% pivot_longer("Factor 1":"Factor 2", names_to = "Factor"), aes(x = x, y = y, color = value)) + #%>% filter(model != "SFA+")
    geom_point() +
    scale_color_viridis_c() +
    facet_wrap(Factor ~ model, nrow = 2) +
    labs(color = "Value") +
    theme_bw()
ggsave(file = here("application_butterflies", "plots", "unconstrainedordination_LVwithlongitudelatitude.pdf"), width = 8, height = 8)



##----------------------
sessioninfo::session_info()
##----------------------
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 (2021-11-01)
# os       Linux Mint 21.1
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2024-04-25
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version     date (UTC) lib source
# boot          1.3-28      2021-05-03 [4] CRAN (R 4.1.1)
# class         7.3-20      2022-01-13 [4] CRAN (R 4.1.2)
# classInt      0.4-8       2022-09-29 [1] CRAN (R 4.1.2)
# cli           3.6.2       2023-12-11 [1] CRAN (R 4.1.2)
# codetools     0.2-18      2020-11-04 [4] CRAN (R 4.0.3)
# colorspace    2.1-0       2023-01-23 [1] CRAN (R 4.1.2)
# corrplot    * 0.92        2021-11-18 [1] CRAN (R 4.1.2)
# DBI           1.1.3       2022-06-18 [1] CRAN (R 4.1.2)
# DHARMa      * 0.4.6       2022-09-08 [1] CRAN (R 4.1.2)
# doParallel  * 1.0.17      2022-02-07 [1] CRAN (R 4.1.2)
# dplyr         1.1.2       2023-04-20 [1] CRAN (R 4.1.2)
# e1071         1.7-13      2023-02-01 [1] CRAN (R 4.1.2)
# fansi         1.0.6       2023-12-08 [1] CRAN (R 4.1.2)
# fmesher       0.1.5       2023-12-20 [1] CRAN (R 4.1.2)
# foreach     * 1.5.2       2022-02-02 [1] CRAN (R 4.1.2)
# generics      0.1.3       2022-07-05 [1] CRAN (R 4.1.2)
# ggplot2       3.4.3       2023-08-14 [1] CRAN (R 4.1.2)
# glue          1.7.0       2024-01-09 [1] CRAN (R 4.1.2)
# gtable        0.3.4       2023-08-21 [1] CRAN (R 4.1.2)
# here        * 1.0.1       2020-12-13 [1] CRAN (R 4.1.2)
# INLA        * 24.02.09    2024-02-09 [1] local
# iterators   * 1.0.14      2022-02-05 [1] CRAN (R 4.1.2)
# KernSmooth    2.23-20     2021-05-03 [4] CRAN (R 4.0.4)
# lattice       0.20-45     2021-09-22 [4] CRAN (R 4.1.1)
# lifecycle     1.0.4       2023-11-07 [1] CRAN (R 4.1.2)
# lme4          1.1-34.9000 2023-08-25 [1] Github (lme4/lme4@ff5ecd0)
# magrittr      2.0.3       2022-03-30 [1] CRAN (R 4.1.2)
# MASS          7.3-55      2022-01-13 [4] CRAN (R 4.1.2)
# Matrix      * 1.6-1       2023-08-14 [1] CRAN (R 4.1.2)
# minqa         1.2.5       2022-10-19 [1] CRAN (R 4.1.2)
# munsell       0.5.0       2018-06-12 [1] CRAN (R 4.1.2)
# nlme          3.1-162     2023-01-31 [1] CRAN (R 4.1.2)
# nloptr        2.0.3       2022-05-26 [1] CRAN (R 4.1.2)
# patchwork   * 1.1.2       2022-08-19 [1] CRAN (R 4.1.2)
# pillar        1.9.0       2023-03-22 [1] CRAN (R 4.1.2)
# pkgconfig     2.0.3       2019-09-22 [1] CRAN (R 4.1.2)
# proxy         0.4-27      2022-06-09 [1] CRAN (R 4.1.2)
# R6            2.5.1       2021-08-19 [1] CRAN (R 4.1.2)
# Rcpp          1.0.12      2024-01-09 [1] CRAN (R 4.1.2)
# rlang         1.1.3       2024-01-10 [1] CRAN (R 4.1.2)
# rprojroot     2.0.4       2023-11-05 [1] CRAN (R 4.1.2)
# rstudioapi    0.15.0      2023-07-07 [1] CRAN (R 4.1.2)
# scales        1.2.1       2022-08-20 [1] CRAN (R 4.1.2)
# sessioninfo   1.2.2       2021-12-06 [1] CRAN (R 4.1.2)
# sf            1.0-9       2022-11-08 [1] CRAN (R 4.1.2)
# sp          * 2.1-3       2024-01-30 [1] CRAN (R 4.1.2)
# tibble        3.2.1       2023-03-20 [1] CRAN (R 4.1.2)
# tidyselect    1.2.0       2022-10-10 [1] CRAN (R 4.1.2)
# TMB         * 1.9.9       2023-11-28 [1] CRAN (R 4.1.2)
# units         0.8-1       2022-12-10 [1] CRAN (R 4.1.2)
# utf8          1.2.4       2023-10-22 [1] CRAN (R 4.1.2)
# vctrs         0.6.5       2023-12-01 [1] CRAN (R 4.1.2)
# withr         3.0.0       2024-01-16 [1] CRAN (R 4.1.2)
# 
# [1] /home/fkch/R/x86_64-pc-linux-gnu-library/4.1
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────



