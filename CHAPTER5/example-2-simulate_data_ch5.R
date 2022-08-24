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
  simdata <- DGP_TRI_FUNCTION( # [1] "PROCESS" VARIABLES
                                N = 1440, # SAMPLE SIZE
                                n = 25, # NUMBER OF LONGITUDINAL MEASUREMENTS
                                CA = 25, # ADMINISTRATIVE CENOSRING
                                RandomCensoring = FALSE, # IS RANDOM CENSORING ALLOWED?
                                prop.censored = 0.12, # THE PROPORTION OF PEOPLE CENSORED
                                rho = 1, # THE MULITSTATE VISIT FREQUENCY PARAMETER
                                
                                # [2] RANDOM EFFECTS - WE ARE ASSUMING 6, ALL OF WHICH ARE CORRELATED
                                D11 = 1.4017*1.4017, 
                                D21 = -0.04346, 
                                D22 = 0.2248*0.2248,
                                D31 = 0.2565,
                                D32 = 0.01517,
                                D33 = 0.6314*0.6314,
                                D41 = -0.01313,
                                D42 = 0.008311,
                                D43 = -0.00639,
                                D44 = 0.07677*0.07677,
                                D51 = -1.8694,
                                D52 = 0.7317,
                                D53 = 1.2070,
                                D54 = 0.2288,
                                D55 = 11.4293*11.4293,
                                D61 = -0.00610,
                                D62 = -0.09324,
                                D63 = -0.05309,
                                D64 = -0.05461,
                                D65 = -4.2248,
                                D66 = 0.9635*0.9635,
                                
                                # [3] COVARIATES
                                X1mean = 9.07,
                                X1sd = 1.6,
                                X2p <- 0.5,
                                X3p <- 0.5,
                                
                                # [4] LINEAR MIXED-EFFECTS MODEL PARAMETERS
                                sigma.ETDRS = 1.1362,
                                sigma.AER = 0.6687,
                                sigma.EGFR = 8.3159,
                                true.betas.ETDRS = c(2.0799, 0.1920, 0.1527, -2.3148),
                                true.betas.AER = c(2.0727, 0.02547, 0.05528, -0.3228),
                                true.betas.EGFR = c(120.27, -1.4380, 0.4161, 1.2858),
                                
                                # [5] BASELINE INTENSITY FUNCTION PARAMETERS/CONSTANTS
                                pwc_knot1 = 4.8,
                                pwc_knot2 = 17.4182,
                                lambda = c("lambda.121" = -2.8, "lambda.122" = -3.6, "lambda.123" = -3.0,
                                           "lambda.231" = -4.5, "lambda.232" = -5.2, "lambda.233" = -5.7),
                                # [6] ADDITIONAL TRANSITION INTENSITY PARAMETERS
                                gammas <- c("X1cont" = 0.07, "X3cat" = -0.27, 
                                            "X1cont" = 0.10),
                                # [7] ASSOCIATION PARAMETERS
                                alpha12 = c("alpha.121" = 0.11, "alpha.122" = 0.23, "alpha.123" = -0.005), 
                                alpha23 = c("alpha.231" = 0.05, "alpha.232" = 0.15, "alpha.233" = -0.01)
  )
  simdata$sim_number <- k
  simlist[[k]] <- simdata
  states[k, ] <- .Random.seed
}
simdata <- do.call("rbind", simlist)

# SAVE THE DATASET (OPTIONAL SINCE IT IS ALREADY INCLUDED IN THE FILES)
setwd('M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5')
saveRDS(simdata,'simdata_DGP1_n1.rds')

# OPTIONAL CODE TO STORE ADDITIONAL RESULTS, E.G. WHEN SETTING nsim EQUAL TO 200
# saveRDS(states,'states_DGP1_n100.rds')
# saveRDS(simdata,'simdata_DGP1_n100.rds')




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SETTINGS FOR DGP 5
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THIS REQUIRES A DIFFERENT FUNCTION THAN THAT USED FOR DGP 1
set.seed(1234) 
DGP_TRI_FUNCTION_2( # [1] "PROCESS" VARIABLES
                    N = 1440, # SAMPLE SIZE
                    n = 25, # NUMBER OF LONGITUDINAL MEASUREMENTS
                    CA = 25, # ADMINISTRATIVE CENOSRING
                    RandomCensoring = FALSE, # IS RANDOM CENSORING ALLOWED?
                    prop.censored = 0.12, # THE PROPORTION OF PEOPLE CENSORED
                    rho = 1, # THE MULITSTATE VISIT FREQUENCY PARAMETER
                    
                    # [2] RANDOM EFFECTS - WE ARE ASSUMING 6, ALL OF WHICH ARE CORRELATED
                    D11 = 1.8275*1.8275,
                    D21 = -0.045,
                    D22 = 0.225,
                    D31 = 0.47,
                    D32 = 0.015,
                    D33 = 0.6566*0.6566,
                    D41 = -0.019,
                    D42 = 0.008,
                    D43 = -0.007,
                    D44 = 0.07663*0.07663,
                    D51 = -2.56,
                    D52 = 0.736,
                    D53 = 1.10,
                    D54 = 0.23,
                    D55 = 11.4447*11.4447,
                    D61 = -0.005,
                    D62 = -0.093,
                    D63 = -0.053,
                    D64 = -0.054,
                    D65 = -4.23,
                    D66 = 0.9639*0.9639,
                    
                    # [3] COVARIATES
                    X1mean = 9.07,
                    X1sd = 1.6,
                    X2p <- 0.5,
                    X3p <- 0.5,
                    
                    # [4] LINEAR MIXED-EFFECTS MODEL PARAMETERS
                    sigma.ETDRS = 1.10,
                    sigma.AER = 0.65,
                    sigma.EGFR = 8.16,
                    true.betas.ETDRS = c(0.37, 0.19 , 0.21) ,
                    true.betas.AER = c(1.51, 0.02, 0.10),
                    true.betas.EGFR = c(118.59, -1.39, 0.65),
                    
                    # [5] BASELINE INTENSITY FUNCTION PARAMETERS/CONSTANTS
                    pwc_knot1 = 4.8,
                    pwc_knot2 = 17.4182,
                    lambda = c("lambda.121" = -2.8, "lambda.122" = -3.6, "lambda.123" = -3.0,
                               "lambda.231" = -4.5, "lambda.232" = -5.2, "lambda.233" = -5.7),
                    # [6] ADDITIONAL TRANSITION INTENSITY PARAMETERS
                    gammas <- c("X1cont" = 0.07, "X3cat" = -0.27, 
                                "X1cont" = 0.10),
                    # [7] ASSOCIATION PARAMETERS
                    alpha12 = c("alpha.121" = 0.11, "alpha.122" = 0.23, "alpha.123" = -0.005), 
                    alpha23 = c("alpha.231" = 0.05, "alpha.232" = 0.15, "alpha.233" = -0.01)
)