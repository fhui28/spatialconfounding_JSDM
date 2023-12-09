#' ---
#' title: Analysis of Butterflies data from Ovaskainen et al., (2016) using spatial factor analytic models
#' abstract: Idea is to explore whether spatial confounding is present in the model, and compare how different models handle SC. 
#' Note the models are set up such that the *species-specific coefficients to the measured environmental predictors are treated as random slopes* ala community-level pool. This is done to bring some stability into the model, otherwise treating them as fixed effects has the tendency to cause issues with standard error calculation.
#' author: Francis KC Hui
#' date: Code started July 2023
#' ---

##-----------------------
#' # Load in packages and explore data
##-----------------------
rm(list = ls())
library(tidyverse)
library(patchwork) 
library(TMB) 
library(INLA) 
library(DHARMa) 
library(mgcv) 
library(foreach) 
library(doParallel) 
library(corrplot) 
library(kableExtra) 
library(Matrix) 
library(mustashe) 
registerDoParallel(cores = 6)
here::i_am("application_butterflies/analysis.R")
library(here)


resp_dat <- read.csv(file = "Butterfly_PA.csv")
str(resp_dat)
colSums(resp_dat)
resp_dat <- resp_dat[,colSums(resp_dat) > 150] 
str(resp_dat)
colSums(resp_dat)

cov_dat <- read.csv(file = "Butterfly_Cov.csv")
longlat <- read.csv(file = "Butterfly_LatLon_scaled.csv")
cov_dat <- cbind(longlat, cov_dat[,-1]) # Remove intercept
cov_dat$blwood <- scale(cov_dat$blwood) %>% as.vector
cov_dat$conwood <- scale(cov_dat$conwood) %>% as.vector
str(cov_dat) 

summary(cov_dat) # Note latitude and longitudes have been scaled and so do not make sense here; so treat them as Euclidean coordinates/distance metric, as per https://doi.org/10.1111/2041-210X.13106
boxplot(cov_dat)


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
ggsave(file = "plots/covariates.pdf", width = 8, height = 4)

#GGally::ggpairs(cov_dat)

ggplot(cbind(longlat, resp_dat)  %>% pivot_longer(-(x:y), names_to = "species"), aes(x = x, y = y, color = value %>% factor)) +
    geom_point() +
    facet_wrap(. ~ species, nrow = 4, scales = "free") +
    scale_color_viridis_d() +
    theme_bw()


#' ## Is there evidence of residual spatial correlation?
num_spp <- resp_dat %>% ncol
# for(k0 in 1:num_spp) {
#     fit_cw <- glm(response ~ climate + blwood + conwood, family = binomial(), #+ chalk_limestone
#                   data = data.frame(response = resp_dat[,k0], cov_dat))
#     dotestSAC <- simulateResiduals(fit_cw) %>% testSpatialAutocorrelation(., x = cov_dat$x, y = cov_dat$y, plot = FALSE)
#     print(dotestSAC)
#     }
# rm(fit_cw, dotestSAC)
# All species exhibit evidence of residual spatial correlation, after accounting for all four covariates!

dev.off()


##-------------------------
#' # Prepare objects and fit independent GLLVM 
#' Relevant output saved in an RData file, so **do not need to run**
##-------------------------
#' Fitting is done using Template Model Builder (TMB). 
#' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood.

X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2

filename <- "TMB_indFA.cpp"
modelname <- strsplit(filename, "\\.")[[1]][1]
compile(filename)
dyn.load(dynlib(modelname))

tidbits_data <- list(y = resp_dat %>% as.matrix, 
                     X = X, 
                     num_lv = num_lv
                     )

blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
tidbits_parameters <- list(betas = matrix(0, nrow = num_spp, ncol = ncol(X)),
                           alpha = numeric(ncol(X)),
                           log_sds = numeric(ncol(X)),
                           loadings = numeric(sum(lower.tri(blank_loadings, diag = TRUE))),
                           lvs = matrix(0, nrow = nrow(X), ncol = num_lv)
                           )
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
lowerlim[grep("loadings",names(objs$par))][diag(blank_loadings)] <- 1e-6 ## Constraints on loading matrix    
rm(blank_loadings)


#' ## Fit and gather results
fit_lvm <- nlminb(start = objs$par, 
                 objective = objs$fn, 
                 gradient = objs$gr, 
                 lower = lowerlim, 
                 upper = upperlim,
                 #method = "L-BFGS-B",
                 control = list(iter.max = 100000, eval.max = 10000, trace = 1))

indlvm_results <- sdreport(objs) 
pt_estimates <- as.list(indlvm_results, what = "Estimate", report = TRUE)
all_results <- summary(indlvm_results, select = "report", p.value = TRUE)


ci_alpha <- 0.05
betaind_results <- data.frame(all_results[grep("betas$", rownames(all_results)),])
betaind_results$lower <- betaind_results$Estimate - qnorm(1-ci_alpha/2) * betaind_results$Std
betaind_results$upper <- betaind_results$Estimate + qnorm(1-ci_alpha/2) * betaind_results$Std
rownames(betaind_results) <- paste0(rep(colnames(tidbits_data$y), ncol(X)), ":", rep(colnames(X), each = ncol(tidbits_data$y)))
colnames(betaind_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 

alphaind_results <- data.frame(all_results[grep("alpha$", rownames(all_results)),])
alphaind_results$lower <- alphaind_results$Estimate - qnorm(1-ci_alpha/2) * alphaind_results$Std
alphaind_results$upper <- alphaind_results$Estimate + qnorm(1-ci_alpha/2) * alphaind_results$Std
rownames(alphaind_results) <- colnames(X)
colnames(alphaind_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper")

lvind_results <- pt_estimates$lvs
loadingind_results <- pt_estimates$loadings_mat

eta_ind <- list(Xbeta = tcrossprod(X, pt_estimates$betas), 
    residual = tcrossprod(pt_estimates$lvs, pt_estimates$loadings_mat), 
    residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
    rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor
    )    
    
dyn.unload(paste0(modelname,".so"))
#file.remove(filename, paste0(modelname,".so"), paste0(modelname,".o"))
gc()
rm(all_results, pt_estimates, fit_lvm, indlvm_results)


save(betaind_results, 
     eta_ind, 
     loadingind_results, 
     lvind_results, 
     file = "independentspatialfit.RData")



##-------------------------
#' # Prepare objects and fit unrestricted spatial GLLVM/factor analytic (SFA) model. 
#' Relevant objects saved in an RData file at the end, so **do not need to run**
##-------------------------
#' Fitting is done using Template Model Builder (TMB). Since the spatial domain is not very large (Grampians national park) and given the scaling that has taken place prior to accessing the data, then Euclidean coordinates are assumed. 
#' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood, where the Matern covariance is approximated using a SPDE approach.
#' Standard errors obtained using a generalized Delta method based on the joint covariance of the random effects and parameter estimates

X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2

filename <- "TMBINLA_SFA.cpp"
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
                     G2 = spde$param.inla$M2
                     )

blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
tidbits_parameters <- list(betas = matrix(0, nrow = num_spp, ncol = ncol(X)),
                           alpha = numeric(ncol(X)),
                           log_sds = numeric(ncol(X)),
                           loadings = numeric(sum(lower.tri(blank_loadings, diag = TRUE))),
                           lvs = matrix(0, nrow = mesh$n, ncol = num_lv),
                           log_kappa = rep(0, num_lv)
                           )
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
lowerlim[grep("loadings",names(objs$par))][diag(blank_loadings)] <- 1e-6 ## Constraints on loading matrix    
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
               rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor
               )

eta_RSFA <- list(Xbeta = tcrossprod(X, pt_estimates$betas_rsr), 
                 residual = tcrossprod(pt_estimates$lvs_rsr, pt_estimates$loadings_mat),
                 residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
                 rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor
                 )



dyn.unload(paste0(modelname,".so"))
#file.remove(filename, paste0(modelname,".so"), paste0(modelname,".o"))

rm(Loc, objs, tidbits_parameters, ci_alpha, filename, modelname, lowerlim, upperlim, tidbits_random, splvm_results, fit_splvm, all_results, pt_estimates)
gc()

save(betaRSFA_results, betaSFA_results,
     alphaSFA_results,
     lvSFA_results, lvRSFA_results,
     eta_RSFA, eta_SFA,
     loadingSFA_results, loadingRSFA_results,
     file = "spatialfit.RData")



##-------------------------
#' # Prepare objects and fit restricted spatial spatial GLLVM/factor analytic (RSFA) model. 
#' **Not run since it encounters computational issues, plus alphas (which this model was specifically designed to run to get) are not of interest**
##-------------------------
#' Fitting is done using Template Model Builder (TMB). Since the spatial domain is not very large (Grampians national park) and given the scaling that has taken place prior to accessing the data, then Euclidean coordinates are assumed. 
#' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood, where the Matern covariance is approximated using a SPDE approach.
#' Standard errors obtained using a generalized Delta method based on the joint covariance of the random effects and parameter estimates

X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2

filename <- "TMBINLA_RSFA.cpp"
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
                     residualprojection = diag(nrow = nrow(X)) - X %*% tcrossprod(solve(crossprod(X)), X),
                     n_i = mesh$n,
                     meshidxloc=mesh$idx$loc - 1,
                     G0 = spde$param.inla$M0, 
                     G1 = spde$param.inla$M1, 
                     G2 = spde$param.inla$M2
)

blank_loadings <- matrix(0, nrow = num_spp, ncol = num_lv)
tidbits_parameters <- list(betas = matrix(0, nrow = num_spp, ncol = ncol(X)),
                           alpha = numeric(ncol(X)),
                           log_sds = numeric(ncol(X)),
                           loadings = numeric(sum(lower.tri(blank_loadings, diag = TRUE))),
                           lvs = matrix(0, nrow = mesh$n, ncol = num_lv),
                           log_kappa = rep(0, num_lv)
)
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
upperlim <- rep(100, length(objs$par))
lowerlim <- rep(-100, length(objs$par))
lowerlim[grep("loadings",names(objs$par))][diag(blank_loadings)] <- 1e-6 ## Constraints on loading matrix    
rm(blank_loadings)


#' ## Fit and gather results
fit_restrictedsplvm <- nlminb(start = objs$par, 
                    objective = objs$fn, 
                    gradient = objs$gr, 
                    lower = lowerlim, 
                    upper = upperlim, 
                    control = list(iter.max = 10000, eval.max = 10000, trace = 1))


restrictedsplvm_results <- sdreport(objs)
pt_estimates <- as.list(restrictedsplvm_results, what = "Estimate", report = TRUE)
all_results <- summary(restrictedsplvm_results, select = "report", p.value = TRUE)


ci_alpha <- 0.05
betaRSFA_results <- data.frame(all_results[grep("betas$", rownames(all_results)),])
betaRSFA_results$lower <- betaRSFA_results$Estimate - qnorm(1-ci_alpha/2) * betaRSFA_results$Std
betaRSFA_results$upper <- betaRSFA_results$Estimate + qnorm(1-ci_alpha/2) * betaRSFA_results$Std
rownames(betaRSFA_results) <- paste0(rep(colnames(tidbits_data$y), ncol(X)), ":", rep(colnames(X), each = ncol(tidbits_data$y)))
colnames(betaRSFA_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 

alphaRSFA_results <- data.frame(all_results[grep("alpha$", rownames(all_results)),])
alphaRSFA_results$lower <- alphaRSFA_results$Estimate - qnorm(1-ci_alpha/2) * alphaRSFA_results$Std
alphaRSFA_results$upper <- alphaRSFA_results$Estimate + qnorm(1-ci_alpha/2) * alphaRSFA_results$Std
rownames(alphaRSFA_results) <- colnames(X)
colnames(alphaRSFA_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 

lvRSFA_results <- pt_estimates$lvs_rsr
loadingRSFA_results <- pt_estimates$loadings_mat

eta_RSFA <- list(Xbeta = tcrossprod(X, pt_estimates$betas), 
                 residual = tcrossprod(pt_estimates$lvs_rsr, pt_estimates$loadings_mat),
                 residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
                 rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor
                 )


dyn.unload(paste0(modelname,".so"))
#file.remove(filename, paste0(modelname,".so"), paste0(modelname,".o"))

rm(Loc, objs, tidbits_parameters, ci_alpha, filename, modelname, lowerlim, upperlim, tidbits_random, restrictedsplvm_results, fit_restrictedsplvm, all_results, pt_estimates)
gc()

save(betaRSFA_results,
     alphaRSFA_results,
     lvRSFA_results,
     eta_SFA,
     loadingRSFA_results,
     file = "restrictedspatialfit.RData")



##-------------------------
#' # Prepare objects and fit Spatial+ GLLVM: **NOT USED**
#' Note it is not doing Spatial+ properly for non-Gaussian responses as in the paper. The first step of Spatial+ involves fitting a weighted spatial model to each covariate, with weights from a spatial model of the response itself. However, in a JSDM this causes because the residuals of the covariates from this first stage Spatial+ will then be species-specific. So in the second stage of Spatial+, the X's themselves will be response specific. *So for now our implementation ignores this and is experimental*
#' Using thin plate splines with k = 300 everywhere, as per [https://doi.org/10.1111/biom.13656]. This lead to very LVs to the SFA model i.e., they were still correlated with measured covariates. On the other hand, the effects of the covariate were more diminished with many no longer statistically significant.
#' Now trying default number of basis functions. 
#' Relevant objects saved in an RData file at the end, so **do not need to run**
##-------------------------
X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2

#' ## Construct residuals from a spatial model for the strongly spatial/more continuous covariate
Xfits <- foreach(k0 = 2:ncol(X)) %dopar% gam(covariate ~ s(x, y), data = data.frame(covariate = X[,k0], cov_dat))
lapply(Xfits, summary) # The covariates are not entirely spatial, even climate. This is good, as it suggests there is potential residual signal for spatial+ to work with?!
lapply(Xfits, gam.check) # Still needs more according to mgcv

resX <- cbind(X[,1], sapply(Xfits, residuals, type = "response"))
colnames(resX) <- colnames(X)
summary(resX)     
GGally::ggpairs(resX[,-1] %>% as.data.frame())
ggplot(bind_cols(longlat, resX[,-1]) %>% pivot_longer(-(x:y), names_to = "predictors"), aes(x = x, y = y, color = value)) +
    geom_point() +
    facet_wrap(. ~ predictors, nrow = 4, scales = "free") +
    scale_color_viridis_c() +
    theme_bw()


#' ## Now fit the response model      
filename <- "TMBINLA_SFA.cpp"
modelname <- strsplit(filename, "\\.")[[1]][1]
compile(filename)
dyn.load(dynlib(modelname))


# Create SPDE mesh
Loc <- cov_dat %>% dplyr::select(x:y)
mesh <- inla.mesh.create(Loc, plot.delay = NULL, refine = FALSE)
spde <- inla.spde2.matern(mesh, alpha = 2)


tidbits_data <- list(y = resp_dat %>% as.matrix, 
                     X = resX, 
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
lowerlim[grep("loadings",names(objs$par))][diag(blank_loadings)] <- 1e-6 ## Constraints on loading matrix    
rm(blank_loadings)


#' ## Fit and gather results
fit_spplus <- nlminb(start = objs$par, 
                    objective = objs$fn, 
                    gradient = objs$gr, 
                    lower = lowerlim, 
                    upper = upperlim, 
                    control = list(iter.max = 10000, eval.max = 10000, trace = 1))


spplus_results <- sdreport(objs)
pt_estimates <- as.list(spplus_results, what = "Estimate", report = TRUE)
all_results <- summary(spplus_results, select = "report", p.value = TRUE)


ci_alpha <- 0.05
betaspplus_results <- data.frame(all_results[grep("betas$", rownames(all_results)),])
betaspplus_results$lower <- betaspplus_results$Estimate - qnorm(1-ci_alpha/2) * betaspplus_results$Std
betaspplus_results$upper <- betaspplus_results$Estimate + qnorm(1-ci_alpha/2) * betaspplus_results$Std
rownames(betaspplus_results) <- paste0(rep(colnames(tidbits_data$y), ncol(X)), ":", rep(colnames(X), each = num_spp))
colnames(betaspplus_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 

alphaspplus_results <- data.frame(all_results[grep("alpha$", rownames(all_results)),])
alphaspplus_results$lower <- alphaspplus_results$Estimate - qnorm(1-ci_alpha/2) * alphaspplus_results$Std
alphaspplus_results$upper <- alphaspplus_results$Estimate + qnorm(1-ci_alpha/2) * alphaspplus_results$Std
rownames(alphaspplus_results) <- colnames(X)
colnames(alphaspplus_results) <- c("Estimate", "StdErr", "z_value", "P_val", "lower", "upper") 


lvspplus_results <- pt_estimates$lvs_units
loadingspplus_results <- pt_estimates$loadings_mat

eta_spplus <- list(Xbeta = tcrossprod(resX, pt_estimates$betas), 
               residual = tcrossprod(pt_estimates$lvs_units, pt_estimates$loadings_mat),
               residual2 = diag(tcrossprod(pt_estimates$loadings_mat)),
               rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor
               )

               
dyn.unload(paste0(modelname,".so"))
#file.remove(filename, paste0(modelname,".so"), paste0(modelname,".o"))

rm(Loc, tidbits_parameters, ci_alpha, filename, modelname, lowerlim, upperlim, tidbits_random, fit_spplus, all_results)
gc()


save(betaspplus_results, 
     eta_spplus,
     loadingspplus_results,
     lvspplus_results, 
     file = "spatialplusfit_kdefault.RData")

     

##------------------------------
#' # Explore results for covariate effects and LVs
##------------------------------
load(file = "independentspatialfit.RData")
load(file = "spatialfit.RData")
#load(file = "spatialplusfit_kdefault.RData")
#load(file = "spatialplusfit_k300.RData")

X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2

tidbits_data <- list(y = resp_dat %>% as.matrix, 
                     X = X, 
                     num_lv = num_lv
                     )

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


#' ## Plots of 95% Wald intervals for covariate effects
results_dat <- rbind(
                     betaind_results[-(1:num_spp),], 
                     betaSFA_results[-(1:num_spp),], 
                     betaRSFA_results[-(1:num_spp),]
                     #betaspplus_results[-(1:num_spp),]
                     ) %>%
    rownames_to_column() %>%
    mutate(sig = P_val < 0.05) %>%
    mutate(species = rep(rep(colnames(tidbits_data$y), ncol(X)-1), 3) %>% fct_inorder) %>%
    mutate(predictors = rep(rep(colnames(X)[-1], each = num_spp), 3) %>% fct_inorder) %>%
    mutate(model = rep(c("Independent", "SFA", "RSFA"), each = (ncol(X)-1)*num_spp) %>% fct_inorder) 
levels(results_dat$predictors) <- c("GDD", "Broadleaved woodland", "Coniferous woodland")
levels(results_dat$species) <- spp_names
    
    
myColors <- c("slateblue","coral3")
names(myColors) <- levels(results_dat$sig)

ggplot(data = results_dat, 
       aes(x = species, y = Estimate, col = sig, Lower = lower, Upper = upper, group = model, shape = model)) +
   geom_hline(yintercept = 0, linetype = 2) + 
   geom_errorbar(aes(ymin = lower, ymax = upper, col = sig), width = 0.1, position = position_dodge(width = 0.4)) + 
   geom_point(fill = "white", position = position_dodge(width = 0.4)) + 
   facet_grid(. ~ predictors, scales = "free") + 
   scale_colour_manual(values = myColors) +
   scale_shape_manual(values=c(21:24)) +
   labs(color = "Significant", x = "Species", y = "Estimated regression coefficients", shape = "Effect type") + 
   theme_bw() +
   theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "bottom") +
   coord_flip() +
   guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave(file = "plots/caterpillarcoefficients.pdf", width = 8, height = 10)


#' ## Tile plot of regression  coefficients only
tile_climate <- ggplot(results_dat %>% filter(predictors == "GDD"), #%>% mutate(Estimate2 = ifelse(sig == 1, Estimate, 0)), 
       aes(x = model, y = species, fill = Estimate)) +
   geom_tile(linetype = 0) + 
   geom_tile() +
   #facet_wrap(. ~ predictors, nrow = 1) +
   labs(x = "Model", y = "Species", fill = "Estimated \ncoefficients", title = "GDD") +
   scale_fill_viridis_c() +
   coord_flip() +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

tile_bwoodland <- ggplot(results_dat %>% filter(predictors == "Broadleaved woodland"), #%>% mutate(Estimate2 = ifelse(sig == 1, Estimate, 0)), 
       aes(x = model, y = species, fill = Estimate)) +
   geom_tile(linetype = 0) + 
   geom_tile() +
   #facet_wrap(. ~ predictors, nrow = 1) +
   labs(x = "Model", y = "Species", fill = "Estimated \ncoefficients", title = "Broadleaved woodland") +
   scale_fill_viridis_c() +
   coord_flip() +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

tile_cwoodland <- ggplot(results_dat %>% filter(predictors == "Coniferous woodland"), #%>% mutate(Estimate2 = ifelse(sig == 1, Estimate, 0)), 
       aes(x = model, y = species, fill = Estimate)) +
   geom_tile(linetype = 0) + 
   geom_tile() +
   #facet_wrap(. ~ predictors, nrow = 1) +
   labs(x = "Model", y = "Species", fill = "Estimated \ncoefficients", title = "Coniferous woodland") +
   scale_fill_viridis_c() +
   coord_flip() +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- (tile_climate / tile_bwoodland / tile_cwoodland)
p

ggsave(file = "plots/tilecoefficients.pdf", width = 10, height = 12)

rm(p)


results_dat %>% 
    group_by(predictors, model) %>% 
    summarise(mean = mean(sig))

results_dat %>% 
    group_by(predictors, model) %>% 
    summarise(mean = mean(sig)) %>% 
   kbl(format = "latex", caption = "", digits = 3, 
       col.names = c("Predictor", "Model", "Proportion of significant results"), 
       align = "lrr",
       label = "tab:butterfly_propsigresults",
       booktabs = TRUE,
       centering = TRUE)
#' There are some but not too many differences between SFA, RSFA, and independence in terms of whether or not there is statistically clear evidence of covariate effects for different species. Independent and RSFA models behave more similarly to each other compared to FA, with statistically clear evidence that climate and blwood are important for many species, and less so for conwood. SFA however declares much less species significant for all three predictors. 
#' For conwood, RSFAS declares the most number of species as statistically significant, many more than SFA and independence model. 
# # Finally SFA+ is probably the most different out of the four methods, as it tends to declare many things to be not statistically significant. Arguably it behaves closest to SFA, especially for blwood and which is to be expected, but it is also more conservative than all other methods, especially for climate.


ggplot(results_dat, aes(x = model, y = StdErr)) +
   geom_boxplot() +
   geom_point() +
   ggbeeswarm::geom_beeswarm() +
   labs(x = "Model", y = "Standard Errors") +
   facet_wrap(predictors ~ ., nrow = 1, scales = "free") +
   theme_bw()
ggsave(file = "plots/intervalwidths.pdf", width = 8, height = 5)
#' Interestingly the widths are SFA intervals are generally narrower for climate and blwood. Interval widths are more similar for independent and RSFA models, and if anything the standard errors are largest for the independence model.  SFA+ has the least variable of the CI widths, much smaller than all three other methods.

ggplot(results_dat %>% 
           mutate(width = upper - lower) %>% 
           dplyr::select(species:width) %>% 
           pivot_wider(names_from = model, values_from = width) %>% 
           mutate(ratio = RSFA/SFA),
       aes(x = predictors, y = ratio)) +
    geom_boxplot() +
    theme_bw() +
    labs(x = "Predictors", y = "Ratio (RSFA/SFA)") +
    geom_hline(yintercept = 1, linetype = 2, col = "darkgrey")

results_dat %>% 
           dplyr::select(StdErr, species:model) %>% 
           pivot_wider(names_from = model, values_from = StdErr) %>% 
           mutate(ratioRSFA = RSFA/SFA, ratioindependent = Independent/SFA) %>%
           group_by(predictors) %>%
           summarise(meanRSFA = mean(ratioRSFA), meanIndependent = mean(ratioindependent))
# Interestingly the intervals are quite similar in this dataset between RSFA and SFA. However these results are not necessarily generalizable, especially given the known results in the existing literature. 



#' ## Plots of point estimates for latent variables
results_dat <- rbind(lvind_results, lvSFA_results, lvRSFA_results) %>% #lvspplus_results
    as.data.frame()
colnames(results_dat)[1:num_lv] <- paste("Factor", 1:num_lv)
results_dat$model <- rep(c("independent", "SFA", "RSFA"), each = nrow(X)) %>% fct_inorder #"SFA+"
results_dat$x <- rep(cov_dat$x, 3)
results_dat$y <- rep(cov_dat$y, 3)

allcors <- rbind(
   data.frame(Model = "Independent", LV = paste("Factor", 1:2), cor(lvind_results, X[,-1])), 
   data.frame(Model = "SFA", LV = paste("Factor", 1:2), cor(lvSFA_results, X[,-1])), 
   data.frame(Model = "RSFA", LV = paste("Factor", 1:2), cor(lvRSFA_results, X[,-1]))
   #data.frame(Model = "Spatial+", LV = paste("Factor", 1:2), cor(lvspplus_results, X[,-1]))
    ) %>% 
   pivot_longer(climate:conwood, names_to = "Predictors") %>% 
   mutate(Predictors = fct_inorder(Predictors), Model = fct_inorder(Model))
levels(allcors$Predictors) <- c("Climate", "Broadleaved woodland", "Coniferous woodland")
allcors <- allcors %>% 
   arrange(Predictors, Model, LV) %>% 
   relocate(Predictors, Model, LV, value)

allcors %>% 
    as.data.frame()


kbl(allcors, format = "latex", caption = "", digits = 3, 
    col.names = c("Predictors", "Model", "Factor", "Correlation"), 
    align = "lllr",
    label = "tab:lvcorrelation",
    booktabs = TRUE,
    centering = TRUE)
#' Notice the SFA and SFA+ predicted LVs are strongly correlated with climate and to a lesser extent with blwood, and the correlations are not small (suggesting confounding is present to some extent). RSFA forces orthgonality. The predicted LVs from the independence model also exhibit fairly much weaker correlation with the covariates. 
#' Finally, the predicted LVs from the SFA+ model plus exhibits strong correlations with the three predictors, even more so than the SFA model. This is perhaps not surprising, given it SFA+ is actually designed to make the LVs orthogonal to the residual projected predictors, which is in fact what we see. This reflects what SFA+ in essence tries to do i.e., move all the spatial variation in X to the spatial random effect...


ggplot(results_dat %>% pivot_longer("Factor 1":"Factor 2", names_to = "Factor"), aes(x = model, y = value)) +
    geom_violin() +
    #ggbeeswarm::geom_beeswarm(color = "darkgrey") +
    geom_hline(yintercept = 0, color = "darkblue", linetype = 2) +
    facet_wrap(. ~ Factor) +
    theme_bw()

results_dat %>% 
    pivot_longer("Factor 1":"Factor 2", names_to = "Factor") %>%
    group_by(model, Factor) %>%
    summarise(mean = mean(value) %>% round(3), median = median(value) %>% round(3))
# Shows that RSR tends to produce more centered LVs for residual ordination, as expected. By contrast, the SFA and SFA+ models tend to exhibit a overall strong non-zero location effect. We can come back to this issue for basic unconstrained ordination...

levels(results_dat$model)[1] <- "Independent"
ggplot(results_dat %>% pivot_longer("Factor 1":"Factor 2", names_to = "Factor"), aes(x = x, y = y, color = value)) + #%>% filter(model != "SFA+")
    geom_point() +
    scale_color_viridis_c() +
    facet_wrap(Factor ~ model, nrow = 2) +
   labs(color = "Value") +
    theme_bw()
ggsave(file = "plots/LVwithlongitudelatitude.pdf", width = 8, height = 8)
# The north-south spatial pattern of the predicted LVs from SFA quite obvious, especially for the first LV. The similarlity to the pattern to climate is clear. This pattern also occurs but to a much lesser with the SFA+ plus. 
# This north-south pattern is not present in the RSFA or independence model, although RSFA exhibiting a noticeable east-west spatial pattern for LV1 (and analogously but to a lesser extent for LV2 i nthe independence model). Finally, given the above correlations, we do see that the LVs in the SFA+ model exhibits a spatial pattern that contains the spatial patterns of the predictors as well!


ordinationp_climate <- ggplot(results_dat %>% 
           rename(f1 = "Factor 1", f2 = "Factor 2") %>% 
           left_join(., bind_cols(X[,-1], longlat), by = c("x","y")), aes(x = f1, y = f2, color = climate)) + 
    geom_point(alpha = 0.8) +
    scale_color_viridis_c() +
    labs(color = "GDD", x = "Factor 1", y = "Factor 2") +
    facet_wrap(. ~ model, nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom")

ordinationp_blwood <- ggplot(results_dat %>% 
                                  rename(f1 = "Factor 1", f2 = "Factor 2") %>% 
                                  left_join(., bind_cols(X[,-1], longlat), by = c("x","y")), aes(x = f1, y = f2, color = blwood)) + 
    geom_point(alpha = 0.8) +
    scale_color_viridis_c() +
    labs(color = "Broadleaved woodland", x = "Factor 1", y = "Factor 2") +
    facet_wrap(. ~ model, nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom")

ordinationp_conwood <- ggplot(results_dat %>% 
                                  rename(f1 = "Factor 1", f2 = "Factor 2") %>% 
                                  left_join(., bind_cols(X[,-1], longlat), by = c("x","y")), aes(x = f1, y = f2, color = conwood)) + 
    geom_point(alpha = 0.8) +
    scale_color_viridis_c() +
    labs(color = "Coniferous woodland", x = "Factor 1", y = "Factor 2") +
    facet_wrap(. ~ model, nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom")


p <- ordinationp_climate / ordinationp_blwood / ordinationp_conwood
ggsave(p, file = "plots/residualordination.pdf", width = 12, height = 15)

# ggplot(results_dat %>% 
#            rename(f1 = "Factor 1", f2 = "Factor 2") %>% 
#            left_join(., bind_cols(X[,-1], longlat), by = c("x","y")), aes(x = f1, y = f2, z = climate)) + 
#     stat_summary_hex(bins = 50) +
#     scale_color_viridis_c() +
#     facet_wrap(. ~ model, nrow = 1) +
#     theme_bw()



##------------------------------
#' # Explore results for residual covariance and variance partitioning
##------------------------------
summary((eta_SFA$Xbeta + eta_SFA$residual) - (eta_RSFA$Xbeta + eta_RSFA$residual))
# Checks the linear predictor is the same. Good!


varpart_ind <- data.frame(Xbeta = apply(eta_ind$Xbeta, 2, var), LV = eta_ind$residual2) %>% #LV = apply(eta_ind$residual, 2, var)
    apply(., 1, prop.table) %>% 
    data.frame
varpart_SFA <- data.frame(Xbeta = apply(eta_SFA$Xbeta, 2, var), LV = eta_SFA$residual2) %>% #LV = apply(eta_SFA$residual, 2, var)
    apply(., 1, prop.table) %>% 
    data.frame
varpart_RSFA <- data.frame(Xbeta = apply(eta_RSFA$Xbeta, 2, var), LV = eta_RSFA$residual2) %>% #LV = apply(eta_RSFA$residual, 2, var)
    apply(., 1, prop.table) %>% 
    data.frame
# varpart_spplus <- data.frame(Xbeta = apply(eta_spplus$Xbeta, 2, var), LV = eta_spplus$residual2) %>%
#     apply(., 1, prop.table) %>%
#     data.frame
colnames(varpart_ind) <- colnames(varpart_SFA) <- colnames(varpart_RSFA) <- colnames(resp_dat)

round(varpart_ind, 3)
round(varpart_SFA, 3)
round(varpart_RSFA, 3)
#round(varpart_spplus, 3)


# Construct a nice stacked barplot representing the variance partitioning
v_pretty <- rbind(varpart_ind[,order(varpart_SFA[2,])] %>% rownames_to_column(var = "component"), 
                  varpart_SFA[,order(varpart_SFA[2,])] %>% rownames_to_column(var = "component"), 
                  varpart_RSFA[,order(varpart_SFA[2,])] %>% rownames_to_column(var = "component")
                  #varpart_spplus[,order(varpart_sp[2,])] %>% rownames_to_column(var = "component")
                  ) %>% 
    as.data.frame %>% 
    mutate(model = rep(c("Independent","SFA","RSFA"), each = 2)) %>% #"SFA+"
    pivot_longer(-c("component","model"), names_to = "Species") %>% 
    mutate(Species = fct_inorder(Species), model = fct_inorder(model), component = fct_inorder(component))
levels(v_pretty$Species) <- spp_names
levels(v_pretty$component) <- c("Covariates", "Latent factors")


ggplot(v_pretty %>% filter(model != "SFA+"), aes(x = Species, y = value, fill = component)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(x = "Species", y = "Proportion", fill = "Model component") + #Percentage of variance explained
    theme_bw() +
    facet_wrap(. ~ model, nrow = 4) +
    scale_fill_viridis_d() +
   labs(y = "Variance") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 70, hjust = 1))
ggsave(file = "plots/variancepartition.pdf", width = 8, height = 8)
# For almost all species, the independence and RSFA models pushes more of the variance to be explained into measured predictors. In fact the results between these two models are fairly similar in terms of variance partitioning. By contrast, the SFA+ model has basically all the variation in the data explained by the LVs component. That being said, variance partitioning in the SFA+ model does not make much sense though since you are partitioning the model into a part explained by the residual projection of the X. This is not comparable to the other three models, let alone whether it makes sense or not?!


par(mfrow = c(2,2))
corrplot(eta_ind$rescov, type = "lower", diag = FALSE, title = "Independent", mar = c(2,5,2,2))
corrplot(eta_SFA$rescov, type = "lower", diag = FALSE, title = "SFA", mar = c(2,5,2,2))
corrplot(eta_RSFA$rescov, type = "lower", diag = FALSE, title = "RSFA", mar = c(2,5,2,2))
#corrplot(eta_spplus$rescov, type = "lower", diag = FALSE, title = "SFA+", mar = c(2,5,2,2))
# These empirical residual covariances between all three models are fairly similar though. In fact, and defining the residual correlation as based solely on the loading matrix, SFA and RSFA should produce basically the set of results. However, it is important to acknowledge that such a correlation construct may not make much sense in a *spatial* factor analysis and variation thereof though...The SFA+ model produces the most different pattern in terms of ordinations



##----------------------
sessionInfo()
##----------------------
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 21.1
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
# [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C               LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8     LC_MONETARY=en_AU.UTF-8    LC_MESSAGES=en_AU.UTF-8   
# [7] LC_PAPER=en_AU.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#     [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] here_1.0.1        mustashe_0.1.4    kableExtra_1.3.4  corrplot_0.92     doParallel_1.0.17 iterators_1.0.14  mgcv_1.9-0        nlme_3.1-155      DHARMa_0.4.6     
# [10] INLA_22.12.16     sp_1.6-0          foreach_1.5.2     Matrix_1.5-3      TMB_1.9.2         patchwork_1.1.2   lubridate_1.9.2   forcats_1.0.0     stringr_1.5.0    
# [19] dplyr_1.1.2       purrr_1.0.1       readr_2.1.4       tidyr_1.3.0       tibble_3.2.1      ggplot2_3.4.2     tidyverse_2.0.0  
# 
# loaded via a namespace (and not attached):
# [1] httr_1.4.5         viridisLite_0.4.2  splines_4.1.2      vipor_0.4.5        pillar_1.9.0       lattice_0.20-45    glue_1.6.2         digest_0.6.31     
# [9] RColorBrewer_1.1-3 rvest_1.0.3        minqa_1.2.5        colorspace_2.1-0   htmltools_0.5.4    plyr_1.8.8         pkgconfig_2.0.3    scales_1.2.1      
# [17] webshot_0.5.5      svglite_2.1.1      tzdb_0.3.0         lme4_1.1-33        timechange_0.2.0   generics_0.1.3     farver_2.1.1       ellipsis_0.3.2    
# [25] withr_2.5.0        cli_3.6.1          magrittr_2.0.3     evaluate_0.21      GGally_2.1.2       fansi_1.0.4        MASS_7.3-55        xml2_1.3.3        
# [33] beeswarm_0.4.0     textshaping_0.3.6  tools_4.1.2        hms_1.1.2          lifecycle_1.0.3    munsell_0.5.0      compiler_4.1.2     systemfonts_1.0.4 
# [41] rlang_1.1.1        grid_4.1.2         nloptr_2.0.3       rstudioapi_0.14    labeling_0.4.2     rmarkdown_2.20     boot_1.3-28        gtable_0.3.3      
# [49] codetools_0.2-18   reshape_0.8.9      R6_2.5.1           knitr_1.42         fastmap_1.1.1      utf8_1.2.3         rprojroot_2.0.3    ragg_1.2.5        
# [57] stringi_1.7.12     ggbeeswarm_0.7.2   Rcpp_1.0.10        vctrs_0.6.2        tidyselect_1.2.0   xfun_0.37         
