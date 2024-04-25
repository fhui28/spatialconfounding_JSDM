#' @title Simulate multivariate abundance data from a joint species distribution model involving up to two observed covariates (x1,x2) and two unobserved covariates (z1, z2)
#'
#' @param family The distribution of the responses for the model. 
#' @param x1/x2 Vectors of the two observed covariates. The second one can be optional
#' @param X_coefficients A matrix of species-specific coefficients corresponding to the observed covariates, where the number of rows must be equal to the number of species and the number of columns must be equal to two or three (corresponding to intercept, x1, and x2 if provided).
#' @param z1/z2 Vectors of the two unobserved covariates.
#' @param Z_coefficients A matrix of species-specific coefficients corresponding to the unobserved covariates, where the number of rows must be equal to the number of species and the number of columns must be equal to two (z1, z2).
#' @param dispparam Vector of species-specific dispersion parameters, if required.
#' @param trial_size The trial size used for simulating binomial responses.
#' @param powerparam Vector of species-specific power parameters, if required.
#' @param seed Seed for simulating data.
#' 
#' @return A matrix of multivariate abundance data.

library(tweedie)

simdat_jsdm_simple <- function(family = gaussian(), 
                                      x1, x2 = NULL, X_coefficients, 
                                      z1, z2, Z_coefficients, 
                                      dispparam = NULL, trial_size = 1, powerparam = NULL, 
                                      seed = NULL) {
    
    #' # Set up 
    set.seed(seed)
    
    if(!(family %in% c("gaussian","binomial","poisson","negative_binomial","tweedie")))
        stop("Family not supported. Sorry!")
    if(nrow(X_coefficients) != nrow(Z_coefficients))
        stop("Number of rows in X_coefficients and Z_coefficients should be the same, both equal to the number of species.")
    if(ncol(Z_coefficients) != 2)
        stop("Number of columns in Z_coefficients must be two.")
    if(family %in% c("gaussian", "negative_binomial", "tweedie")) {
        if(length(dispparam) != nrow(X_coefficients))
            stop("Length of dispparam should be equal to the number of rows in X_coefficients i.e., the number of species.")
        }
    if(family %in% c("tweedie")) {
        if(length(powerparam) != nrow(X_coefficients))
            stop("Length of powerparam should be equal to the number of rows in X_coefficients i.e., the number of species.")
        }
    
    if(is.null(rownames(X_coefficients)))
        rownames(X_coefficients) <- rownames(Z_coefficients) <- paste0("spp", 1:nrow(X_coefficients))
    
    X <- as.matrix(cbind(1, x1, x2))
    if(ncol(X_coefficients) != ncol(X))
        stop("Number of columns in X_coefficients is not compatible with the number of x1 and x2 supplied.")
    num_units <- nrow(X)
    num_spp <- nrow(X_coefficients)
         
    sim_y <- matrix(0, nrow = num_units, ncol = num_spp)
    colnames(sim_y) <- rownames(X_coefficients)

     
    #' # Make linear predictor
    missing_predictors <- as.matrix(cbind(z1, z2))
    all_eta <- tcrossprod(X, as.matrix(X_coefficients)) + tcrossprod(missing_predictors, as.matrix(Z_coefficients))

    #' ## Generate responses
    for(j in 1:num_spp) {
        if(family == "gaussian")
            sim_y[,j] <- rnorm(n = num_units, mean = all_eta[,j], sd = sqrt(dispparam[j]))
        if(family == "binomial") 
            sim_y[,j] <- rbinom(n = num_units, size = trial_size, prob = plogis(all_eta[,j]))
        if(family == "poisson") 
            sim_y[,j] <- rpois(n = num_units, lambda = exp(all_eta[,j]))
        if(family == "negative_binomial") 
            sim_y[,j] <- rnbinom(n = num_units, mu = exp(all_eta[,j]), size = 1/dispparam[j])
        if(family == "tweedie") 
            sim_y[,j] <- rtweedie(n = num_units, mu = exp(all_eta[,j]), phi = dispparam[j], p = powerparam[j])
        }
    

    set.seed(NULL)
    return(sim_y)
    }
	

