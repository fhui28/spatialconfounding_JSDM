#' ---
#' title: Analysis of Butterflies data from Ovaskainen et al., (2016) using factor analytic models to explore the presence and impacts of spatial confounding.
#' author: xxx
#' date: Code started July 2023
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
here::i_am("application_butterflies/analysis.R")
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
X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2

filename <- here("models", "TMB_indFA.cpp")
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
                  DLL = "TMB_indFA", 
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
    rescov = tcrossprod(pt_estimates$loadings_mat) %>% cov2cor)    
    
dyn.unload(paste0(modelname,".so"))
gc()
rm(all_results, pt_estimates, fit_lvm, indlvm_results)


save(betaind_results, 
     alphaind_results,
     eta_ind, 
     loadingind_results, 
     lvind_results, 
     file = here("application_butterflies", "independentspatialfit.RData"))



##-------------------------
#' # Prepare objects and fit spatial factor analytic (SFA) model. 
#' Fitting is done using Template Model Builder (TMB). Since the spatial domain is not very large (Grampians national park) and given the scaling that has taken place prior to accessing the data, then Euclidean coordinates are assumed. 
#' Maximum likelihood estimation via TMB of the Laplace approximated log-likelihood, where the Matern covariance is approximated using a SPDE approach.
##-------------------------
X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2

filename <- here("models", "TMB_SFA.cpp")
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
                  DLL = "TMB_SFA", 
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
     file = here("application_butterflies", "spatialfit.RData"))



##-------------------------
#' # Explore results for covariate effects
##-------------------------
load(file = here("application_butterflies", "independentspatialfit.RData"))
load(file = here("application_butterflies", "spatialfit.RData"))

X <- model.matrix(~ climate + blwood + conwood, data = cov_dat) 
num_lv <- 2
num_spp <- ncol(resp_dat)

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
ggsave(file = here("application_butterflies", "plots", "caterpillarcoefficients.pdf"), width = 8, height = 10)


#' ## Tile plot of regression coefficients only
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

ggsave(file = here("application_butterflies", "plots", "tilecoefficients.pdf"), width = 10, height = 12)
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
ggsave(file = here("application_butterflies", "plots", "variancepartition.pdf"), width = 8, height = 8)


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
results_dat$x <- rep(cov_dat$x, 3)
results_dat$y <- rep(cov_dat$y, 3)

allcors <- rbind(
   data.frame(Model = "Independent", LV = paste("Factor", 1:2), cor(lvind_results, X[,-1])), 
   data.frame(Model = "SFA", LV = paste("Factor", 1:2), cor(lvSFA_results, X[,-1])), 
   data.frame(Model = "RSFA", LV = paste("Factor", 1:2), cor(lvRSFA_results, X[,-1]))) %>% 
   pivot_longer(climate:conwood, names_to = "Predictors") %>% 
   mutate(Predictors = fct_inorder(Predictors), Model = fct_inorder(Model))
levels(allcors$Predictors) <- c("Climate", "Broadleaved woodland", "Coniferous woodland")

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
ggplot(results_dat %>% pivot_longer("Factor 1":"Factor 2", names_to = "Factor"), aes(x = x, y = y, color = value)) + 
    geom_point() +
    scale_color_viridis_c() +
    facet_wrap(Factor ~ model, nrow = 2) +
   labs(color = "Value") +
    theme_bw()
ggsave(file = here("application_butterflies", "plots", "LVwithlongitudelatitude.pdf"), width = 8, height = 8)


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
ggsave(p, file = here("application_butterflies", "plots", "residualordination.pdf"), width = 12, height = 15)



##-------------------------
sessionInfo()
##-------------------------
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
