library(tidyverse)
library(MASS)

#-----------------------------------------------------------------------------------------------------
# Example of simulating data based on:
# Gasparini, A. (2018). Rsimsum: Summarise results from monte carlo simulation studies. Journal of Open
# Source Software, 3, 739. https://doi.org/10.21105/joss.00739
#-----------------------------------------------------------------------------------------------------

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SIMULATE 1 DATASET FOR DGP 1
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(1234)
nsim <- 1
simlist <- list()
states <- matrix(ncol = 626, nrow = nsim)
for (k in 1:nsim){
  simdata <- DGP_FUNCTION_3(N = 1440,
                            n = 25, D11 = 3.3, D21 = -0.05, D22 = 0.05, true.betas = c(0.39, 0.19, 0.21), sigma.e = 1.0977,
                            alpha = c("alpha.12" = 0.12, "alpha.23" = 0.06), 
                            pwc_knot1 = 4.8, pwc_knot2 = 17.4182, 
                            lambda = c("lambda.121" = -3.5, "lambda.122" = -5, "lambda.123" = -3.9, 
                                       "lambda.231" = -5.7, "lambda.232" = -5.9, "lambda.233" = -6.5), 
                            CA = 25, RandomCensoring = FALSE, prop.censored = 0.12,
                            rho = 1,
                            X1mean = 9.07,
                            X1sd = 1.6,
                            X2min = 1, X2max = 15, X3p <- 0.5,
                            gammas <- c("X1cont" = 0.08, "X3cat" = -0.32, 
                                        "X1cont" = 0.13)
  )
  simdata$sim_number <- k
  simlist[[k]] <- simdata
  states[k, ] <- .Random.seed
}
sasdata <- do.call("rbind", simlist)

# SAVE THE DATASET (OPTIONAL SINCE IT IS ALREADY INCLUDED IN THE FILES)
setwd('M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER4')
saveRDS(sasdata,'sasdata_DGP1_n1.rds')

# OPTIONAL CODE TO STORE ADDITIONAL RESULTS, E.G. WHEN SETTING nsim EQUAL TO 200
# saveRDS(states,'states_DGP1_n200.rds')
# saveRDS(sasdata,'sasdata_DGP1_n200.rds')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SETTINGS FOR DGP 2
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(2234)
DGP_FUNCTION_3(N = 1440,
               n = 25, D11 = 3.3, D21 = -0.05, D22 = 0.05, true.betas = c(0.39, 0.19, 0.21), sigma.e = 1.0977,
               alpha = c("alpha.12" = 0.12, "alpha.23" = 0.06), 
               pwc_knot1 = 4.8, pwc_knot2 = 17.4182, 
               lambda = c("lambda.121" = -3.5, "lambda.122" = -5, "lambda.123" = -3.9, 
                          "lambda.231" = -5.7, "lambda.232" = -5.9, "lambda.233" = -6.5), 
               CA = 25, RandomCensoring = FALSE, prop.censored = 0.12,
               rho = 0.5,
               X1mean = 9.07,
               X1sd = 1.6,
               X2min = 1, X2max = 15, X3p <- 0.5,
               gammas <- c("X1cont" = 0.08, "X3cat" = -0.32, 
                           "X1cont" = 0.13)
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SETTINGS FOR DGP 3
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(3234)
DGP_FUNCTION_3(N = 1440,
               n = 25, D11 = 3.3, D21 = -0.05, D22 = 0.05, true.betas = c(0.39, 0.19, 0.21), sigma.e = 1.0977,
               alpha = c("alpha.12" = 0.12, "alpha.23" = 0.09), 
               pwc_knot1 = 4.8, pwc_knot2 = 17.4182, 
               lambda = c("lambda.121" = -3.4, "lambda.122" = -4.9, "lambda.123" = -3.8, 
                          "lambda.231" = -5.1, "lambda.232" = -5.3, "lambda.233" = -5.9), 
               CA = 25, RandomCensoring = FALSE, prop.censored = 0.12,
               rho = 1,
               X1mean = 9.07,
               X1sd = 1.6,
               X2min = 1, X2max = 15, X3p = 0.5,
               gammas = c("X1cont" = 0.08, "X3cat" = -0.32, 
                          "X1cont" = 0.13)
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SETTINGS FOR DGP 4
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(4234)
DGP_FUNCTION_3(N = 1440,
               n = 25, D11 = 3.3, D21 = -0.05, D22 = 0.05, true.betas = c(0.39, 0.19, 0.21), sigma.e = 1.0977,
               alpha = c("alpha.12" = 0.12, "alpha.23" = 0.06), 
               pwc_knot1 = 4.8, pwc_knot2 = 17.4182, 
               lambda = c("lambda.121" = -3.5, "lambda.122" = -5, "lambda.123" = -3.9, 
                          "lambda.231" = -5.7, "lambda.232" = -5.9, "lambda.233" = -6.5), 
               CA = 25, RandomCensoring = TRUE, prop.censored = 0.12,
               rho = 1,
               X1mean = 9.07,
               X1sd = 1.6,
               X2min = 1, X2max = 15, X3p = 0.5,
               gammas = c("X1cont" = 0.08, "X3cat" = -0.32, 
                          "X1cont" = 0.13)
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SETTINGS FOR DGP 5
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THIS REQUIRES A DIFFERENT FUNCTION THAN THAT USED FOR DGP 1, 2, 3, AND 4
set.seed(5234)
DGP_FUNCTION_4(N = 1440,
               n = 25, D11 = 3.3, D21 = -0.05, D22 = 0.05, true.betas = c(0.39, 0.19, 0.21), sigma.e = 1.0977,
               alpha = c("alpha.12" = 0.12, "alpha.23" = 0.06), 
               pwc_knot1 = 4.8, pwc_knot2 = 17.4182, 
               lambda = c("lambda.121" = -3.5, "lambda.122" = -5, "lambda.123" = -3.9, 
                          "lambda.231" = -5.7, "lambda.232" = -5.9, "lambda.233" = -6.5), 
               CA = 25, RandomCensoring = FALSE, prop.censored = 0.12,
               rho = 1,
               X1mean = 9.07,
               X1sd = 1.6,
               X2min = 1, X2max = 15, X3p = 0.5,
               gammas = c("X1cont" = 0.08, "X3cat" = -0.32, 
                          "X1cont" = 0.13)
)
