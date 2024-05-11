#' ---
#' title: Analysis of simulation results 
#' author: xxx
#' date: Code started April 2024
#' ---


##-----------------------
#' # Load in packages, data, and fits
##-----------------------
rm(list = ls())
library(tidyverse)
library(colorspace)
library(patchwork)
here::i_am("simulation/analysis.R")
library(here)


#' # In preparation for takeoff...
#' 
#' It is assumed the R scripts in the folders binary_N400, negbin_N400, and tweedie_N400 has been run, thus producing 1200 .RData files in each folder: 200 simulated datasets for each combination of scale of the observed covariate (2 choices, coarse-scaled or fine-scaled) and \code{scale_X} (3 values of of 0, 0.5, 1).


##------------------------------
#' # Gather results
##------------------------------
num_sites <- 400
num_spp <- 10
num_datasets <- 200
scale_X_seq <- c(0,0.5,1)
response_seq <- c("binary", "negbin", "tweedie")
scalecovariate_seq <- c("Coarse-scale observed covariate", "Fine-scale observed covariate")

bias_results <- RMSE_results <- ci_results <- intervallength_results <- array(NA, 
                                                                              dim = c(3, 2, length(scale_X_seq), num_spp, 3),
                                                                              dimnames = list(response = response_seq, 
                                                                                              scalecovariate = scalecovariate_seq,
                                                                                              scale_X = scale_X_seq, 
                                                                                              species = paste0("species", 1:num_spp), 
                                                                                              methods = c("Independent", "SFA", "RSFA")))

LVerror_results <- array(NA,
                         dim = c(3, 2, length(scale_X_seq), num_datasets, 3),
                         dimnames = list(response = response_seq, 
                                         scalecovariate = scalecovariate_seq,
                                         scale_X = scale_X_seq, 
                                         datasets = 1:num_datasets,
                                         methods = c("Independent", "SFA", "RSFA")))
LV_X_correlation_results <- array(NA, 
                                  dim = c(3, 2, length(scale_X_seq), num_datasets, 3),
                                  dimnames = list(response = response_seq, 
                                             scalecovariate = scalecovariate_seq,
                                             scale_X = scale_X_seq, 
                                             datasets = 1:num_datasets,
                                             methods = c("Independent", "SFA", "RSFA")))



for(k0 in 1:length(response_seq)) { for(k1 in 1:length(scale_X_seq)) { 
    #' ## Coarse scale X
    all_diff_X_coefficients <- all_coverage_X_coefficients <- all_length_X_coefficients <- array(NA, dim = c(num_datasets, num_spp, 3)) 
    all_LVerror <- array(NA, dim = c(num_datasets, 3)) 
    all_LV_X_correlation <- array(NA, dim = c(num_datasets, 3)) 
        
    for(k2 in 1:num_datasets) {
        make_filename <- here("simulation", paste0(response_seq[k0], "_N", num_sites), paste0("simulation_", response_seq[k0], "_N", num_sites, "_scaleX", scale_X_seq[k1], "_dataset", k2, ".RData"))

        if(!file.exists(make_filename))
            next;
        
        if(file.exists(make_filename)) {
            message("Loading ", make_filename)
            load(file = make_filename)
            }

        all_diff_X_coefficients[k2,,] <- X_coefficients_dat %>% mutate(across(everything(), ~(.x - true))) %>% select(-true) %>% as.matrix
        all_coverage_X_coefficients[k2,,] <- matrix(X_coefficients_ci_dat$coverage_probability, ncol = 3, byrow = FALSE) 
        all_length_X_coefficients[k2,,] <- matrix(X_coefficients_ci_dat$interval_length, ncol = 3, byrow = FALSE) 
        all_LVerror[k2,] <- c(LV_results_dat$Independent$ss, LV_results_dat$SFA$ss, LV_results_dat$RSFA$ss)
        all_LV_X_correlation[k2,] <- c(LVcor_results_dat$Independent %>% abs %>% max,
                                               LVcor_results_dat$SFA %>% abs %>% max,
                                               LVcor_results_dat$RSFA %>% abs %>% max)
        
        rm(X_coefficients_ci_dat, X_coefficients_dat, LV_results_dat, LVcor_results_dat)
        }
    
    bias_results[k0,1,k1,,] <- apply(all_diff_X_coefficients, c(2,3), function(x) mean(x, na.rm = TRUE))
    RMSE_results[k0,1,k1,,] <- apply(all_diff_X_coefficients^2, c(2,3), function(x) sqrt(mean(x, na.rm = TRUE)))
    ci_results[k0,1,k1,,] <- apply(all_coverage_X_coefficients, c(2,3), mean, na.rm = TRUE)
    intervallength_results[k0,1,k1,,] <- apply(all_length_X_coefficients, c(2,3), mean, na.rm = TRUE)

    LVerror_results[k0,1,k1,,] <- all_LVerror
    LV_X_correlation_results[k0,1,k1,,] <- all_LV_X_correlation
    
    rm(all_diff_X_coefficients, all_coverage_X_coefficients, all_length_X_coefficients, all_LVerror, all_LV_X_correlation)
    rm(list = ls(pattern = "beta"))
    rm(list = ls(pattern = "eta"))
    rm(list = ls(pattern = "loading"))
    
    
    #' ## Fine scale X
    all_diff_X_coefficients <- all_coverage_X_coefficients <- all_length_X_coefficients <- array(NA, dim = c(num_datasets, num_spp, 3)) 
    all_LVerror <- array(NA, dim = c(num_datasets, 3)) 
    all_LV_X_correlation <- array(NA, dim = c(num_datasets, 3)) 
    
    for(k2 in 1:num_datasets) {
        make_filename <- here("simulation", paste0(response_seq[k0], "_N", num_sites), "finescaleX", paste0("simulation_", response_seq[k0], "_finescaleX_N", num_sites, "_scaleX", scale_X_seq[k1], "_dataset", k2, ".RData"))
        
        if(!file.exists(make_filename))
            next;
        
        if(file.exists(make_filename)) {
            message("Loading ", make_filename)
            load(file = make_filename)
        }
        
        all_diff_X_coefficients[k2,,] <- X_coefficients_dat %>% mutate(across(everything(), ~(.x - true))) %>% select(-true) %>% as.matrix
        all_coverage_X_coefficients[k2,,] <- matrix(X_coefficients_ci_dat$coverage_probability, ncol = 3, byrow = FALSE) 
        all_length_X_coefficients[k2,,] <- matrix(X_coefficients_ci_dat$interval_length, ncol = 3, byrow = FALSE) 
        all_LVerror[k2,] <- c(LV_results_dat$Independent$ss, LV_results_dat$SFA$ss, LV_results_dat$RSFA$ss)
        all_LV_X_correlation[k2,] <- c(LVcor_results_dat$Independent %>% abs %>% max,
                                       LVcor_results_dat$SFA %>% abs %>% max,
                                       LVcor_results_dat$RSFA %>% abs %>% max)
        
        rm(X_coefficients_ci_dat, X_coefficients_dat, LV_results_dat, LVcor_results_dat)
    }
    
    bias_results[k0,2,k1,,] <- apply(all_diff_X_coefficients, c(2,3), function(x) mean(x, na.rm = TRUE))
    RMSE_results[k0,2,k1,,] <- apply(all_diff_X_coefficients^2, c(2,3), function(x) sqrt(mean(x, na.rm = TRUE)))
    ci_results[k0,2,k1,,] <- apply(all_coverage_X_coefficients, c(2,3), mean, na.rm = TRUE)
    intervallength_results[k0,2,k1,,] <- apply(all_length_X_coefficients, c(2,3), mean, na.rm = TRUE)
    
    LVerror_results[k0,2,k1,,] <- all_LVerror
    LV_X_correlation_results[k0,2,k1,,] <- all_LV_X_correlation
    
    rm(all_diff_X_coefficients, all_coverage_X_coefficients, all_length_X_coefficients, all_LVerror, all_LV_X_correlation)
    rm(list = ls(pattern = "beta"))
    rm(list = ls(pattern = "eta"))
    rm(list = ls(pattern = "loading"))
    } }


##------------------------------
#' # Present results
##------------------------------
choose_response <- "binary"

p_bias <- ggplot(bias_results %>% 
                     as.data.frame.table %>% 
                     mutate(methods = fct_inorder(methods),
                            value = Freq %>% abs,
                            Freq = NULL) %>% 
                     filter(response == choose_response),
                 aes(x = scale_X, y = value, fill = methods)) +
    geom_boxplot(color = "blueviolet") +
    facet_wrap(. ~ scalecovariate, nrow = 1) +
    scale_fill_viridis_d() +
    scale_y_log10() +
    labs(x = expression(Constant~C[0]), y = "Absolute bias", fill = "Methods") +
    theme_bw()


p_rmse <- ggplot(RMSE_results %>% 
                     as.data.frame.table %>% 
                     mutate(methods = fct_inorder(methods),
                            value = Freq,
                            Freq = NULL) %>% 
                     filter(response == choose_response),
                 aes(x = scale_X, y = value, fill = methods)) +
    geom_boxplot(color = "blueviolet") +
    facet_wrap(. ~ scalecovariate, nrow = 1) +
    scale_fill_viridis_d() +
    labs(x = expression(Constant~C[0]), y = "RMSE", linetype = "Methods", fill = "Methods") +
    scale_y_log10() +
    theme_bw()


p_ci <- ggplot(ci_results %>% 
                     as.data.frame.table %>% 
                     mutate(methods = fct_inorder(methods),
                            value = Freq,
                            Freq = NULL) %>% 
                     filter(response == choose_response),
               aes(x = scale_X, y = value, color = methods, fill = methods)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_boxplot(color = "blueviolet") +
    facet_wrap(. ~ scalecovariate, nrow = 1) +
    scale_fill_viridis_d() +
    labs(x = expression(Constant~C[0]), y = "Coverage probability", fill = "Methods") +
    theme_bw()

p_intervalwidth <- ggplot(intervallength_results %>% 
                              as.data.frame.table %>% 
                              mutate(methods = fct_inorder(methods),
                                     value = Freq,
                                     Freq = NULL) %>% 
                              filter(response == choose_response),
                          aes(x = scale_X, y = value, color = methods, fill = methods)) +
    geom_boxplot(color = "blueviolet") +
    facet_wrap(. ~ scalecovariate, nrow = 1, scales = "free_y") +
    scale_fill_viridis_d() +
    scale_y_log10() +
    labs(x = expression(Constant~C[0]), y = "Mean interval width", fill = "Methods") +
    theme_bw()


p_lverr <- ggplot(LVerror_results %>% 
                     as.data.frame.table %>% 
                     mutate(methods = fct_inorder(methods),
                            value = Freq,
                            Freq = NULL) %>% 
                     filter(response == choose_response),
                  aes(x = scale_X, y = value, color = methods, fill = methods)) +
    geom_boxplot(color = "blueviolet") +
    facet_wrap(. ~ scalecovariate, nrow = 1) +
    scale_fill_viridis_d() +
    scale_y_log10() +
    labs(x = expression(Constant~C[0]), y = "Procrustes error", fill = "Methods") +
    theme_bw()


p_lvXcor <- ggplot(LV_X_correlation_results %>% 
                      as.data.frame.table %>% 
                      mutate(methods = fct_inorder(methods),
                             value = Freq,
                             Freq = NULL) %>% 
                      filter(response == choose_response),
                   aes(x = scale_X, y = value, color = methods, fill = methods)) +
    geom_boxplot(color = "blueviolet") +
    facet_wrap(. ~ scalecovariate, nrow = 1) +
    scale_fill_viridis_d() +
    labs(x = expression(Constant~C[0]), y = "Max. correlation between LV and covariate", fill = "Methods") +
    theme_bw()



(p_bias / p_rmse / p_ci) + plot_layout(guides = "collect", axis_titles =  "collect_x") & theme(legend.position = "bottom")

ggsave(file = here("simulation", "plots", paste0(choose_response, "_N400_coefficientresults.pdf")), width = 8, height = 12)


p_intervalwidth & theme(legend.position = "bottom")

ggsave(file = here("simulation", "plots", paste0(choose_response, "_N400_intervalwidth.pdf")), width = 8, height = 6)


(p_lverr | p_lvXcor) & theme(legend.position = "bottom")

ggsave(file = here("simulation", "plots", paste0(choose_response, "_LVresults.pdf")), width = 12, height = 6)



##----------------------
sessioninfo::session_info()
##----------------------
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 (2021-11-01)
# os       Linux Mint 21.1
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2024-04-26
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# cli           3.6.2   2023-12-11 [1] CRAN (R 4.1.2)
# colorspace  * 2.1-0   2023-01-23 [1] CRAN (R 4.1.2)
# dplyr       * 1.1.2   2023-04-20 [1] CRAN (R 4.1.2)
# fansi         1.0.6   2023-12-08 [1] CRAN (R 4.1.2)
# forcats     * 1.0.0   2023-01-29 [1] CRAN (R 4.1.2)
# generics      0.1.3   2022-07-05 [1] CRAN (R 4.1.2)
# ggplot2     * 3.4.3   2023-08-14 [1] CRAN (R 4.1.2)
# glue          1.7.0   2024-01-09 [1] CRAN (R 4.1.2)
# gtable        0.3.4   2023-08-21 [1] CRAN (R 4.1.2)
# here        * 1.0.1   2020-12-13 [1] CRAN (R 4.1.2)
# hms           1.1.3   2023-03-21 [1] CRAN (R 4.1.2)
# lifecycle     1.0.4   2023-11-07 [1] CRAN (R 4.1.2)
# lubridate   * 1.9.2   2023-02-10 [1] CRAN (R 4.1.2)
# magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.1.2)
# munsell       0.5.0   2018-06-12 [1] CRAN (R 4.1.2)
# patchwork   * 1.2.0   2024-01-08 [1] CRAN (R 4.1.2)
# pillar        1.9.0   2023-03-22 [1] CRAN (R 4.1.2)
# pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.1.2)
# purrr       * 1.0.2   2023-08-10 [1] CRAN (R 4.1.2)
# R6            2.5.1   2021-08-19 [1] CRAN (R 4.1.2)
# readr       * 2.1.4   2023-02-10 [1] CRAN (R 4.1.2)
# rlang         1.1.3   2024-01-10 [1] CRAN (R 4.1.2)
# rprojroot     2.0.4   2023-11-05 [1] CRAN (R 4.1.2)
# rstudioapi    0.15.0  2023-07-07 [1] CRAN (R 4.1.2)
# scales        1.2.1   2022-08-20 [1] CRAN (R 4.1.2)
# sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.1.2)
# stringi       1.7.12  2023-01-11 [1] CRAN (R 4.1.2)
# stringr     * 1.5.0   2022-12-02 [1] CRAN (R 4.1.2)
# tibble      * 3.2.1   2023-03-20 [1] CRAN (R 4.1.2)
# tidyr       * 1.3.0   2023-01-24 [1] CRAN (R 4.1.2)
# tidyselect    1.2.0   2022-10-10 [1] CRAN (R 4.1.2)
# tidyverse   * 2.0.0   2023-02-22 [1] CRAN (R 4.1.2)
# timechange    0.2.0   2023-01-11 [1] CRAN (R 4.1.2)
# tzdb          0.3.0   2022-03-28 [1] CRAN (R 4.1.2)
# utf8          1.2.4   2023-10-22 [1] CRAN (R 4.1.2)
# vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.1.2)
# withr         3.0.0   2024-01-16 [1] CRAN (R 4.1.2)
# 
# [1] /home/fkch/R/x86_64-pc-linux-gnu-library/4.1
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
