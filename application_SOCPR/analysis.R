#' ---
#' title: Analysis of SO-CPR count data from Hui (2020) using spatial factor analytic models. 
#' author: xxx
#' date: Code started April 2024
#' ---

#' # In preparation for takeoff...
#' 
#' The data are publicly available from (https://data.aad.gov.au/aadc/cpr/).
#' 
#' The script below assumes these data has been tidied, subsetted, and saved into an .RData file called \code{socpr_reduced.RData} which has been loaded into the working directory. The .RData file is assumed to contain the following two objects: 
#' - \code{covariate_dat}, which is a data frame containing the spatial coordinates (longitude and latitude) and three covariates (bathymetry, slope, salinity). For the application below, we will treat the longitude and latitude as Euclidean coordinates. 
#' - \code{resp_dat}, which a data frame containing the recorded species counts.
#' 
#' For illustration purposes, we have provide \code{socpr_reduced.RData} file in the Github repository that **contains synthetic versions** of \code{covariate_dat} and \code{resp_dat}.


##-----------------------
#' # Load in packages and explore data
##-----------------------
rm(list = ls())
library(tidyverse)
library(patchwork) 
library(sp)
library(TMB) 
library(INLA) 
library(DHARMa) 
#library(mgcv) 
library(foreach) 
library(doParallel) 
library(corrplot) 
#library(kableExtra) 
library(Matrix) 
#library(mustashe) 
registerDoParallel(cores = 6)
here::i_am("application_SOCPR/analysis.R")
library(here)

load(file = "socpr_reduced.RData")

summary(covariate_dat)
covariate_dat <- covariate_dat %>% 
    as.data.frame %>% 
    mutate(Bathymetry = scale(Bathymetry) %>% as.vector,
           Slope = scale(Slope) %>% as.vector,
           Salinity = scale(Salinity) %>% as.vector) 
str(covariate_dat) 

summary(covariate_dat) #' For the application below, we will treat the longitude and latitude as Euclidean coordinates.
boxplot(covariate_dat %>% select(Bathymetry:Salinity))


#' ## Some simple spatial plots
ggplot(covariate_dat %>% pivot_longer(-(Longitude:Latitude), names_to = "predictors"), aes(x = Longitude, y = Latitude, color = value)) +
    geom_point() +
    facet_wrap(. ~ predictors, nrow = 4, scales = "free") +
    scale_color_viridis_c() +
    theme_bw()


covariate_datlatlong <- covariate_dat %>% 
    pivot_longer(-(Longitude:Latitude), names_to = "predictors") %>% 
    mutate(predictors = fct_inorder(predictors)) 

ggplot(covariate_datlatlong, aes(x = Longitude, y = Latitude, color = value)) +
    geom_point() +
    facet_wrap(. ~ predictors, nrow = 1, scales = "free") +
    scale_color_viridis_c() +
    labs(color = "Value") +
    theme_bw()
ggsave(file = here("application_SOCPR", "plots", "covariates.pdf"), width = 8, height = 4)


ggplot(cbind(covariate_dat %>% select(Longitude:Latitude), resp_dat) %>% 
           pivot_longer(-(Longitude:Latitude), names_to = "species"), aes(x = Longitude, y = Latitude, color = value)) +
    geom_point() +
    facet_wrap(. ~ species, nrow = 4, scales = "free") +
    scale_color_viridis_c(trans = "log1p") +
    theme_bw() # Not really too useful given the number of species...


#' ## Is there evidence of residual spatial correlation?
num_spp <- resp_dat %>% ncol
for(k0 in 1:num_spp) {
    fit_cw <- MASS::glm.nb(response ~ Bathymetry + Slope + Salinity,
                  data = data.frame(response = resp_dat[,k0], covariate_dat))
    dotestSAC <- simulateResiduals(fit_cw) %>% testSpatialAutocorrelation(., x = covariate_dat$Longitude, y = covariate_dat$Latitude, plot = FALSE)
    print(dotestSAC)
    }
rm(fit_cw, dotestSAC, k0)
# In the real dataset, all species exhibit evidence of residual spatial correlation, after accounting for all three covariates!

dev.off()


##-------------------------
#' # Prepare objects and fit independent GLLVM 
#' Fitting is done using Template Model Builder (TMB) -- Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood.
##-------------------------
X <- model.matrix(~ Bathymetry + Slope + Salinity, data = covariate_dat) 
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
                           log_dispparam = numeric(num_spp),
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
lowerlim[grep("loadings",names(objs$par))][diag(blank_loadings)] <- 1e-6 ## Constraints on loading matrix    
rm(blank_loadings)


#' ## Fit and gather results
fit_lvm <- nlminb(start = objs$par, 
                  objective = objs$fn, 
                  gradient = objs$gr, 
                  lower = lowerlim, 
                  upper = upperlim,
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
     eta_ind, 
     loadingind_results, 
     lvind_results, 
     file = "independentspatialfit.RData")



##-------------------------
#' # Prepare objects and fit spatial factor analytic (SFA) model. 
#' Fitting is done using Template Model Builder (TMB). Since the spatial domain is not very large (Grampians national park) and given the scaling that has taken place prior to accessing the data, then Euclidean coordinates are assumed. 
#' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood, where the Matern covariance is approximated using a SPDE approach.
##-------------------------
X <- model.matrix(~ Bathymetry + Slope + Salinity, data = covariate_dat) 
num_lv <- 2

filename <- "TMBINLA_SFA.cpp"
modelname <- strsplit(filename, "\\.")[[1]][1]
compile(filename)
dyn.load(dynlib(modelname))


# Create SPDE mesh
Loc <- covariate_dat %>% dplyr::select(Longitude:Latitude)
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
                           log_dispparam = numeric(num_spp),
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
     lvSFA_results, lvRSFA_results,
     eta_RSFA, eta_SFA,
     loadingSFA_results, loadingRSFA_results,
     file = "spatialfit.RData")


##-------------------------
#' # Explore results for covariate effects
##-------------------------
load(file = here("application_SOCPR", "independentspatialfit.RData"))
load(file = here("application_SOCPR", "spatialfit.RData"))

X <- model.matrix(~ Bathymetry + Slope + Salinity, data = covariate_dat) 
num_lv <- 2
num_spp <- ncol(resp_dat)

tidbits_data <- list(y = resp_dat %>% as.matrix, 
                     X = X, 
                     num_lv = num_lv)

spp_names <- colnames(resp_dat)



#' ## Plots of 95% Wald intervals for covariate effects
#' ## Plots of 95% Wald intervals for covariate effects
results_dat <- rbind(
    betaind_results[-(1:num_spp),], 
    betaSFA_results[-(1:num_spp),], 
    betaRSFA_results[-(1:num_spp),]) %>%
    rownames_to_column() %>%
    mutate(sig = P_val < 0.05) %>%
    mutate(species = rep(rep(colnames(tidbits_data$y), ncol(X)-1), 3) %>% fct_inorder) %>%
    mutate(predictors = rep(rep(colnames(X)[-1], each = num_spp), 3) %>% fct_inorder) %>%
    mutate(model = rep(c("Independent", "SFA", "RSFA"), each = (ncol(X)-1)*num_spp) %>% fct_inorder) 
levels(results_dat$predictors) <- c("Bathymetry", "Slope", "Salinity")
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
tile_bathy <- ggplot(results_dat %>% filter(predictors == "Bathymetry"), #%>% mutate(Estimate2 = ifelse(sig == 1, Estimate, 0)), 
                       aes(x = model, y = species, fill = Estimate)) +
    geom_tile(linetype = 0) + 
    geom_tile() +
    #facet_wrap(. ~ predictors, nrow = 1) +
    labs(x = "Model", y = "Species", fill = "Estimated \ncoefficients", title = "Bathymetry") +
    scale_fill_viridis_c() +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

tile_slope <- ggplot(results_dat %>% filter(predictors == "Slope"), #%>% mutate(Estimate2 = ifelse(sig == 1, Estimate, 0)), 
                         aes(x = model, y = species, fill = Estimate)) +
    geom_tile(linetype = 0) + 
    geom_tile() +
    #facet_wrap(. ~ predictors, nrow = 1) +
    labs(x = "Model", y = "Species", fill = "Estimated \ncoefficients", title = "Slope") +
    scale_fill_viridis_c() +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

tile_salinity <- ggplot(results_dat %>% filter(predictors == "Salinity"), #%>% mutate(Estimate2 = ifelse(sig == 1, Estimate, 0)), 
                         aes(x = model, y = species, fill = Estimate)) +
    geom_tile(linetype = 0) + 
    geom_tile() +
    #facet_wrap(. ~ predictors, nrow = 1) +
    labs(x = "Model", y = "Species", fill = "Estimated \ncoefficients", title = "Salinity") +
    scale_fill_viridis_c() +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- (tile_bathy / tile_slope / tile_salinity)
p

ggsave(file = here("application_SOCPR", "plots", "tilecoefficients.pdf"), width = 10, height = 12)
rm(p)


results_dat %>% 
    group_by(predictors, model) %>% 
    summarise(mean = mean(sig))


##-------------------------
#' # Explore results for residual between-species correlation and variance partitioning
##-------------------------
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
colnames(varpart_ind) <- colnames(varpart_SFA) <- colnames(varpart_RSFA) <- colnames(resp_dat)

round(varpart_ind, 3)
round(varpart_SFA, 3)
round(varpart_RSFA, 3)


# Construct a nice stacked barplot representing the variance partitioning
v_pretty <- rbind(varpart_ind[,order(varpart_SFA[2,])] %>% rownames_to_column(var = "component"), 
                  varpart_SFA[,order(varpart_SFA[2,])] %>% rownames_to_column(var = "component"), 
                  varpart_RSFA[,order(varpart_SFA[2,])] %>% rownames_to_column(var = "component")) %>% 
    as.data.frame %>% 
    mutate(model = rep(c("Independent","SFA","RSFA"), each = 2)) %>% 
    pivot_longer(-c("component","model"), names_to = "Species") %>% 
    mutate(Species = fct_inorder(Species), model = fct_inorder(model), component = fct_inorder(component))
levels(v_pretty$Species) <- spp_names
levels(v_pretty$component) <- c("Covariates", "Latent factors")


ggplot(v_pretty, aes(x = Species, y = value, fill = component)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(x = "Species", y = "Proportion", fill = "Model component") + 
    theme_bw() +
    facet_wrap(. ~ model, nrow = 4) +
    scale_fill_viridis_d() +
    labs(y = "Variance") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 70, hjust = 1))
ggsave(file = here("application_SOCPR", "plots", "variancepartition.pdf"), width = 8, height = 8)

par(mfrow = c(2,2))
corrplot(eta_ind$rescov, type = "lower", diag = FALSE, title = "Independent", mar = c(2,5,2,2))
corrplot(eta_SFA$rescov, type = "lower", diag = FALSE, title = "SFA", mar = c(2,5,2,2))
corrplot(eta_RSFA$rescov, type = "lower", diag = FALSE, title = "RSFA", mar = c(2,5,2,2))
# These empirical residual covariances between all three models are fairly similar though. In fact, and defining the residual correlation as based solely on the loading matrix, SFA and RSFA should produce basically the set of results. However, it is important to acknowledge that such a correlation construct may not make much sense in a *spatial* factor analysis and variation thereof though.


##-------------------------
#' # Explore results for model-based residual ordination
##-------------------------
results_dat <- rbind(lvind_results, lvSFA_results, lvRSFA_results) %>%
    as.data.frame()
colnames(results_dat)[1:num_lv] <- paste("Factor", 1:num_lv)
results_dat$model <- rep(c("independent", "SFA", "RSFA"), each = nrow(X)) %>% fct_inorder 
results_dat$Longitude <- rep(covariate_dat$Longitude, 3)
results_dat$Latitude <- rep(covariate_dat$Latitude, 3)

allcors <- rbind(
    data.frame(Model = "Independent", LV = paste("Factor", 1:2), cor(lvind_results, X[,-1])), 
    data.frame(Model = "SFA", LV = paste("Factor", 1:2), cor(lvSFA_results, X[,-1])), 
    data.frame(Model = "RSFA", LV = paste("Factor", 1:2), cor(lvRSFA_results, X[,-1]))) %>% 
    pivot_longer(Bathymetry:Salinity, names_to = "Predictors") %>% 
    mutate(Predictors = fct_inorder(Predictors), Model = fct_inorder(Model))

allcors <- allcors %>% 
    arrange(Predictors, Model, LV) %>% 
    relocate(Predictors, Model, LV, value) %>% 
    as.data.frame()
allcors    


results_dat %>% 
    pivot_longer("Factor 1":"Factor 2", names_to = "Factor") %>%
    group_by(model, Factor) %>%
    summarise(mean = mean(value) %>% round(3), median = median(value) %>% round(3))

levels(results_dat$model)[1] <- "Independent"
ggplot(results_dat %>% pivot_longer("Factor 1":"Factor 2", names_to = "Factor"), aes(x = Longitude, y = Latitude, color = value)) + 
    geom_point() +
    scale_color_viridis_c() +
    facet_wrap(Factor ~ model, nrow = 2) +
    labs(color = "Value") +
    theme_bw()
ggsave(file = here("application_SOCPR", "plots", "LVwithlongitudelatitude.pdf"), width = 8, height = 8)


ordinationp_bathymetry <- ggplot(results_dat %>% 
                                  rename(f1 = "Factor 1", f2 = "Factor 2") %>% 
                                  left_join(., bind_cols(X[,-1], covariate_dat %>% select(Longitude:Latitude)), by = c("Longitude","Latitude")), 
                                 aes(x = f1, y = f2, color = Bathymetry)) + 
    geom_point(alpha = 0.8) +
    scale_color_viridis_c() +
    labs(color = "Bathymetry", x = "Factor 1", y = "Factor 2") +
    facet_wrap(. ~ model, nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom")

ordinationp_slope <- ggplot(results_dat %>% 
                                 rename(f1 = "Factor 1", f2 = "Factor 2") %>% 
                                left_join(., bind_cols(X[,-1], covariate_dat %>% select(Longitude:Latitude)), by = c("Longitude","Latitude")), 
                            aes(x = f1, y = f2, color = Slope)) + 
    geom_point(alpha = 0.8) +
    scale_color_viridis_c() +
    labs(color = "Slope", x = "Factor 1", y = "Factor 2") +
    facet_wrap(. ~ model, nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom")

ordinationp_salinity <- ggplot(results_dat %>% 
                                  rename(f1 = "Factor 1", f2 = "Factor 2") %>% 
                                   left_join(., bind_cols(X[,-1], covariate_dat %>% select(Longitude:Latitude)), by = c("Longitude","Latitude")), 
                               aes(x = f1, y = f2, color = Salinity)) + 
    geom_point(alpha = 0.8) +
    scale_color_viridis_c() +
    labs(color = "Salinity", x = "Factor 1", y = "Factor 2") +
    facet_wrap(. ~ model, nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom")


p <- ordinationp_bathymetry / ordinationp_slope / ordinationp_salinity
ggsave(p, file = here("application_SOCPR", "plots", "residualordination.pdf"), width = 12, height = 15)


##-------------------------
sessioninfo::session_info()
##-------------------------
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
