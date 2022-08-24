library(rstanarm)
library(rstan)
library(tidyverse)
library(nlme)


# ----------------------------------------------------------------------------------------------------
# THIS PROGRAM CREATES THE DATASET REQUIRE TO RUN STAN (IN THE PROPER FORMAT), CALLED standata_paper3.rds
# IT ALSO CREATES THE INITIAL VALUES, CALLED inits_paper3.rds
# BOTH OF THESE FILES ARE INCLUDED IN THE CHPATER 5 FILES, SO THIS PROGRAM DOES NOT NEED TO BE RUN
# ----------------------------------------------------------------------------------------------------

# WE REQUIRE SOME FUNCTIONS FROM THE RSTANARM PACKAGE
# THE SOURCE FILES ARE LOCATED HERE: https://github.com/stan-dev/rstanarm
# THESE FUNCTIONS HAVE BEEN COPIED TO THE CHAPTER 5 FOLDER
# WE USED VERSION 2.21.3

setwd("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5")
simdata <- readRDS('simdata_DGP1_n1.rds')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# [1] RUN THE MVMER FUNCTION WITH RSTANARM, TO INITIATE THE DATA
# THIS FUNCTION RUNS A MULTIVARIATE GENERALIZED LINEAR MIXED-EFFECTS MODEL OR OUR 3 LONGITUDINAL OUTCOMES
# WE WILL USE THE OUTPUT AND STANDATA FOR THIS MODEL IN OUR JOINT MULTISTATE MODEL
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df.long <- subset(simdata, outcome_number==2)
mlmm <- stan_mvmer(
  formula = list(
    y1 ~ year + X1cont + X2cat + (year | id), 
    y2 ~ year + X1cont + X2cat + (year | id),
    y3 ~ year + X1cont + X2cat + (year | id)),
  data = df.long,
  chains = 1, cores = 1, seed = 12345, iter = 10) # WE ONLY NEED TO INITIALIZE THE DATA

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# [2] BUILD THE DATA FOR FUTURE USE IN STAN, USING THE RSTANARM FUNCTIONS
# THE CODE MUST BE RUN LINE BY LINE, WITH THE END RSEULTS BEING THE STAN DATA SET
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/AUXILLARY FILES/rstanarm-master/R/jm_data_block.R")
source("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/AUXILLARY FILES/rstanarm-master/R/data_block.R")
source("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/AUXILLARY FILES/rstanarm-master/R/misc.R")
source("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/AUXILLARY FILES/rstanarm-master/R/stan_glm.fit.R")
source("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/AUXILLARY FILES/rstanarm-master/R/stanreg-methods.R")
source("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/AUXILLARY FILES/rstanarm-master/R/stanmodels.R")

mlmm$call
mlmm$algorithm
# ENTER THE FUNCTION PARAMETERS
# formula = list(y ~ year + X1cont + (year | id))
formula = list(
  y1 ~ year + X1cont + X2cat + (year | id), 
  y2 ~ year + X1cont + X2cat + (year | id),
  y3 ~ year + X1cont + X2cat + (year | id))
data = df.long
family = gaussian
weights <- NULL		          
prior = normal(autoscale=TRUE)
prior_intercept = normal(autoscale=TRUE)
prior_aux = cauchy(0, 5, autoscale=TRUE)
prior_covariance = lkj(autoscale=TRUE)
prior_PD = FALSE
# algorithm = c("sampling", "meanfield", "fullrank")
# algorithm = c("sampling")
algorithm = mlmm$algorithm
adapt_delta = NULL
max_treedepth = 10L
init = "random"
QR = FALSE
sparse = FALSE

#-----------------------------
# Pre-processing of arguments
#-----------------------------  

# algorithm <- match.arg(algorithm) # ?????????

if (missing(weights)) weights <- NULL

if (!is.null(weights)) 
  stop("'weights' are not yet implemented.")
if (QR)               
  stop("'QR' decomposition is not yet implemented.")
if (sparse)
  stop("'sparse' option is not yet implemented.")

# Formula
formula <- validate_arg(formula, "formula"); M <- length(formula)
if (M > 3L)
  stop("'stan_mvmer' is currently limited to a maximum of 3 outcomes.")

# Data
data <- validate_arg(data, "data.frame", validate_length = M)  
data <- xapply(formula, data, FUN = get_all_vars) # drop additional vars

# Family
ok_classes <- c("function", "family", "character")
ok_families <- c("binomial", "gaussian", "Gamma", 
                 "inverse.gaussian", "poisson", "neg_binomial_2")
family <- validate_arg(family, ok_classes, validate_length = M)
family <- lapply(family, validate_famlink, ok_families)

# Observation weights
if (!is.null(weights)) {
  if (!is(weights, "list")) 
    weights <- rep(list(weights), M)
  weights <- lapply(weights, validate_weights)
}

# Is prior* already a list?
prior <- broadcast_prior(prior, M)
prior_intercept <- broadcast_prior(prior_intercept, M)
prior_aux <- broadcast_prior(prior_aux, M)

#-----------
# Fit model
#----------- 
# stanfit <- stan_jm.fit(formulaLong = formula, dataLong = data, family = family,
#                        weights = weights, priorLong = prior, 
#                        priorLong_intercept = prior_intercept, priorLong_aux = prior_aux, 
#                        prior_covariance = prior_covariance, prior_PD = prior_PD, 
#                        algorithm = algorithm, adapt_delta = adapt_delta, 
#                        max_treedepth = max_treedepth, init = init, 
#                        QR = QR, sparse = sparse, ...)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# NOW MOVE TO STAN_JM.FIT.R
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
formulaLong = formula
dataLong = data
family = family
weights = weights
priorLong = prior
priorLong_intercept = prior_intercept
priorLong_aux = prior_aux
prior_covariance = prior_covariance
prior_PD = prior_PD
algorithm = algorithm
adapt_delta = adapt_delta
max_treedepth = max_treedepth
init = init
QR = QR
sparse = sparse

formulaEvent <- NULL
dataEvent <- NULL
#-----------------------------
# Pre-processing of arguments
#-----------------------------  

if (!requireNamespace("survival"))
  stop("the 'survival' package must be installed to use this function.")

# Set seed if specified
# dots <- list(...)
dots <- list(seed=12345) # I ALTERED THIS LINE....MAYBE SHOULD HAVE INCLUDED SEED IN THE TOP CALL
if ("seed" %in% names(dots))
  set.seed(dots$seed)

# algorithm <- match.arg(algorithm) # ?????????????
# basehaz   <- match.arg(basehaz) # ????????????????

# if (missing(basehaz_ops)) basehaz_ops <- NULL
# if (missing(weights))     weights     <- NULL
# if (missing(id_var))      id_var      <- NULL
# if (missing(time_var))    time_var    <- NULL
# if (missing(grp_assoc))   grp_assoc   <- NULL

basehaz_ops <- NULL
weights     <- NULL
id_var      <- NULL
time_var    <- NULL
grp_assoc   <- NULL

if (!is.null(weights)) 
  stop("'weights' are not yet implemented.")
if (QR)               
  stop("'QR' decomposition is not yet implemented.")
if (sparse)
  stop("'sparse' option is not yet implemented.")

# Error if args not supplied together
supplied_together(formulaLong, dataLong, error = TRUE)
supplied_together(formulaEvent, dataEvent, error = TRUE)

# Determine whether a joint longitudinal-survival model was specified
is_jm <- supplied_together(formulaLong, formulaEvent)
stub <- if (is_jm) "Long" else "y"

if (is_jm && is.null(time_var))
  stop("'time_var' must be specified.")

# Formula
formulaLong <- validate_arg(formulaLong, "formula"); M <- length(formulaLong)

# Data
dataLong <- validate_arg(dataLong, "data.frame", validate_length = M)  
if (is_jm)
  dataEvent <- as.data.frame(dataEvent)

# Family
ok_classes <- c("function", "family", "character")
ok_families <- c("binomial", "gaussian", "Gamma", 
                 "inverse.gaussian", "poisson", "neg_binomial_2")
family <- validate_arg(family, ok_classes, validate_length = M)
family <- lapply(family, validate_famlink, ok_families)
family <- lapply(family, append_mvmer_famlink)

# Observation weights
has_weights <- !is.null(weights)

# Priors
priorLong <- broadcast_prior(priorLong, M)
priorLong_intercept <- broadcast_prior(priorLong_intercept, M)
priorLong_aux <- broadcast_prior(priorLong_aux, M)

#--------------------------
# Longitudinal submodel(s)
#--------------------------

# Info for separate longitudinal submodels
y_mod <- xapply(formulaLong, dataLong, family, FUN = handle_y_mod)

# Construct single cnms list for all longitudinal submodels
y_cnms  <- fetch(y_mod, "z", "group_cnms")
cnms <- get_common_cnms(y_cnms, stub = stub)
cnms_nms <- names(cnms)
if (length(cnms_nms) > 2L)
  stop("A maximum of 2 grouping factors are allowed.")

# Construct single list with unique levels for each grouping factor
y_flist <- fetch(y_mod, "z", "group_list")
flevels <- get_common_flevels(y_flist)

# Ensure id_var is a valid grouping factor in all submodels
if (is_jm) {
  id_var <- check_id_var(id_var, y_cnms, y_flist)
  id_list <- check_id_list(id_var, y_flist)
  if (!is.null(weights))
    weights <- check_weights(weights, id_var)
}

# Observation weights
y_weights <- lapply(y_mod, handle_weights, weights, id_var)

#----------- Prior distributions -----------# 

# Valid prior distributions
ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                  "laplace", "lasso")  # disallow product normal
ok_intercept_dists <- ok_dists[1:3]
ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
ok_covariance_dists <- c("decov", "lkj")

y_vecs <- fetch(y_mod, "y", "y")     # used in autoscaling
x_mats <- fetch(y_mod, "x", "xtemp") # used in autoscaling

# Note: *_user_prior_*_stuff objects are stored unchanged for constructing 
# prior_summary, while *_prior_*_stuff objects are autoscaled

# Priors for longitudinal submodels
y_links <- fetch(y_mod, "family", "link")
y_user_prior_stuff <- y_prior_stuff <- 
  xapply(priorLong, nvars = fetch(y_mod, "x", "K"), link = y_links,
         FUN = handle_glm_prior, 
         args = list(default_scale = 2.5, ok_dists = ok_dists))

y_user_prior_intercept_stuff <- y_prior_intercept_stuff <- 
  xapply(priorLong_intercept, link = y_links, 
         FUN = handle_glm_prior,
         args = list(nvars = 1, default_scale = 10, 
                     ok_dists = ok_intercept_dists))

y_user_prior_aux_stuff <- y_prior_aux_stuff <- 
  xapply(priorLong_aux, FUN = handle_glm_prior, 
         args = list(nvars = 1, default_scale = 5, link = NULL, 
                     ok_dists = ok_aux_dists))  

b_user_prior_stuff <- b_prior_stuff <- handle_cov_prior(
  prior_covariance, cnms = cnms, ok_dists = ok_covariance_dists)

# Autoscaling of priors
y_prior_stuff <- 
  xapply(y_prior_stuff, response = y_vecs, predictors = x_mats, 
         family = family, FUN = autoscale_prior)
y_prior_intercept_stuff <- 
  xapply(y_prior_intercept_stuff, response = y_vecs,
         family = family, FUN = autoscale_prior)
y_prior_aux_stuff <- 
  xapply(y_prior_aux_stuff, response = y_vecs,
         family = family, FUN = autoscale_prior)
if (b_prior_stuff$prior_dist_name == "lkj") { # autoscale priors for ranef sds
  b_prior_stuff <- split_cov_prior(b_prior_stuff, cnms = cnms, submodel_cnms = y_cnms)
  b_prior_stuff <- xapply(
    cnms_nms, FUN = function(nm) {
      z_mats <- fetch(y_mod, "z", "z", nm)
      xapply(b_prior_stuff[[nm]], response = y_vecs, predictors = z_mats, 
             family = family, FUN = autoscale_prior)
    })
} 

#----------- Data for export to Stan -----------# 

standata <- list(
  M = as.integer(M), 
  has_weights  = as.integer(!all(lapply(weights, is.null))),
  family = fetch_array(y_mod, "family", "mvmer_family"),
  link   = fetch_array(y_mod, "family", "mvmer_link"),
  weights = as.array(numeric(0)), # not yet implemented
  prior_PD = as.integer(prior_PD)
)  

# Offset
Y_offset <- fetch(y_mod, "offset", pad_length = 3)
standata$has_offset <- has_offset <-
  fetch_array(y_mod, "has_offset", pad_length = 3)
standata$y1_offset <- if (has_offset[1]) Y_offset[[1]] else as.array(integer(0))  
standata$y2_offset <- if (has_offset[2]) Y_offset[[2]] else as.array(integer(0))  
standata$y3_offset <- if (has_offset[3]) Y_offset[[3]] else as.array(integer(0)) 

# Dimensions
standata$has_aux <- 
  fetch_array(y_mod, "has_aux", pad_length = 3)
standata$resp_type <- 
  fetch_array(y_mod, "y", "resp_type", pad_length = 3)
standata$intercept_type <- 
  fetch_array(y_mod, "intercept_type", "number", pad_length = 3)
standata$yNobs <- 
  fetch_array(y_mod, "x", "N", pad_length = 3)
standata$yNeta <- 
  fetch_array(y_mod, "x", "N", pad_length = 3) # same as Nobs for stan_mvmer
standata$yK <- 
  fetch_array(y_mod, "x", "K", pad_length = 3)

# Response vectors
Y_integer <- fetch(y_mod, "y", "integer")
standata$yInt1 <- if (M > 0) Y_integer[[1]] else as.array(integer(0))  
standata$yInt2 <- if (M > 1) Y_integer[[2]] else as.array(integer(0))  
standata$yInt3 <- if (M > 2) Y_integer[[3]] else as.array(integer(0)) 

Y_real <- fetch(y_mod, "y", "real")
standata$yReal1 <- if (M > 0) Y_real[[1]] else as.array(double(0)) 
standata$yReal2 <- if (M > 1) Y_real[[2]] else as.array(double(0)) 
standata$yReal3 <- if (M > 2) Y_real[[3]] else as.array(double(0)) 

# Population level design matrices
X <- fetch(y_mod, "x", "xtemp")
standata$yX1 <- if (M > 0) X[[1]] else matrix(0,0,0)
standata$yX2 <- if (M > 1) X[[2]] else matrix(0,0,0)
standata$yX3 <- if (M > 2) X[[3]] else matrix(0,0,0)

X_bar <- fetch(y_mod, "x", "x_bar")
standata$yXbar1 <- if (M > 0) as.array(X_bar[[1]]) else as.array(double(0))
standata$yXbar2 <- if (M > 1) as.array(X_bar[[2]]) else as.array(double(0))
standata$yXbar3 <- if (M > 2) as.array(X_bar[[3]]) else as.array(double(0))

# Data for group specific terms - group factor 1
b1_varname <- cnms_nms[[1L]] # name of group factor 1
b1_nvars <- fetch_(y_mod, "z", "nvars", b1_varname, 
                   null_to_zero = TRUE, pad_length = 3)
b1_ngrps <- fetch_(y_mod, "z", "ngrps", b1_varname)
if (!n_distinct(b1_ngrps) == 1L)
  stop("The number of groups for the grouping factor '", 
       b1_varname, "' should be the same in all submodels.")

standata$bN1 <- b1_ngrps[[1L]] + 1L # add padding for _NEW_ group
standata$bK1 <- sum(b1_nvars)
standata$bK1_len <- as.array(b1_nvars)
standata$bK1_idx <- get_idx_array(b1_nvars)

Z1 <- fetch(y_mod, "z", "z", b1_varname)
Z1 <- lapply(Z1, transpose)
Z1 <- lapply(Z1, convert_null, "matrix")
standata$y1_Z1 <- if (M > 0) Z1[[1L]] else matrix(0,0,0)
standata$y2_Z1 <- if (M > 1) Z1[[2L]] else matrix(0,0,0)
standata$y3_Z1 <- if (M > 2) Z1[[3L]] else matrix(0,0,0)

Z1_id <- fetch(y_mod, "z", "group_list", b1_varname)
Z1_id <- lapply(Z1_id, groups)
Z1_id <- lapply(Z1_id, convert_null, "arrayinteger")
standata$y1_Z1_id <- if (M > 0) Z1_id[[1L]] else as.array(integer(0))
standata$y2_Z1_id <- if (M > 1) Z1_id[[2L]] else as.array(integer(0))
standata$y3_Z1_id <- if (M > 2) Z1_id[[3L]] else as.array(integer(0))

# Data for group specific terms - group factor 2
if (length(cnms) > 1L) {
  # model has a second grouping factor
  b2_varname <- cnms_nms[[2L]] # name of group factor 2
  b2_nvars <- fetch_(y_mod, "z", "nvars", b2_varname, 
                     null_to_zero = TRUE, pad_length = 3)
  b2_ngrps <- fetch_(y_mod, "z", "ngrps", b2_varname)
  if (!n_distinct(b2_ngrps) == 1L)
    stop("The number of groups for the grouping factor '", 
         b2_varname, "' should be the same in all submodels.")
  standata$bN2 <- b2_ngrps[[1L]] + 1L # add padding for _NEW_ group
  standata$bK2 <- sum(b2_nvars)
  standata$bK2_len <- as.array(b2_nvars)
  standata$bK2_idx <- get_idx_array(b2_nvars)
  
  Z2 <- fetch(y_mod, "z", "z", b2_varname)
  Z2 <- lapply(Z2, transpose)
  Z2 <- lapply(Z2, convert_null, "matrix")
  standata$y1_Z2 <- if (M > 0) Z2[[1L]] else matrix(0,0,0)
  standata$y2_Z2 <- if (M > 1) Z2[[2L]] else matrix(0,0,0)
  standata$y3_Z2 <- if (M > 2) Z2[[3L]] else matrix(0,0,0)
  
  Z2_id <- fetch(y_mod, "z", "group_list", b2_varname)
  Z2_id <- lapply(Z2_id, groups)
  Z2_id <- lapply(Z2_id, convert_null, "arrayinteger")
  standata$y1_Z2_id <- if (M > 0) Z2_id[[1L]] else as.array(integer(0))
  standata$y2_Z2_id <- if (M > 1) Z2_id[[2L]] else as.array(integer(0))
  standata$y3_Z2_id <- if (M > 2) Z2_id[[3L]] else as.array(integer(0))
  
} else {
  # no second grouping factor
  standata$bN2 <- 0L
  standata$bK2 <- 0L
  standata$bK2_len <- as.array(rep(0,3L))
  standata$bK2_idx <- get_idx_array(rep(0,3L))
  standata$y1_Z2 <- matrix(0,0,0)
  standata$y2_Z2 <- matrix(0,0,0)
  standata$y3_Z2 <- matrix(0,0,0)
  standata$y1_Z2_id <- as.array(integer(0))
  standata$y2_Z2_id <- as.array(integer(0))
  standata$y3_Z2_id <- as.array(integer(0))
}

# Priors
standata$y_prior_dist_for_intercept <- 
  fetch_array(y_prior_intercept_stuff, "prior_dist")  
standata$y_prior_mean_for_intercept <- 
  fetch_array(y_prior_intercept_stuff, "prior_mean")
standata$y_prior_scale_for_intercept <- 
  fetch_array(y_prior_intercept_stuff, "prior_scale")
standata$y_prior_df_for_intercept <- 
  fetch_array(y_prior_intercept_stuff, "prior_df")

standata$y_prior_dist_for_aux <-
  fetch_array(y_prior_aux_stuff, "prior_dist")
standata$y_prior_mean_for_aux <- 
  fetch_array(y_prior_aux_stuff, "prior_mean")
standata$y_prior_scale_for_aux <- 
  fetch_array(y_prior_aux_stuff, "prior_scale")
standata$y_prior_df_for_aux <- 
  fetch_array(y_prior_aux_stuff, "prior_df")

standata$y_prior_dist <- 
  fetch_array(y_prior_stuff, "prior_dist", pad_length = 3)

prior_mean <- fetch(y_prior_stuff, "prior_mean")
standata$y_prior_mean1 <- if (M > 0) prior_mean[[1]] else as.array(double(0))
standata$y_prior_mean2 <- if (M > 1) prior_mean[[2]] else as.array(double(0))
standata$y_prior_mean3 <- if (M > 2) prior_mean[[3]] else as.array(double(0))

prior_scale <- fetch(y_prior_stuff, "prior_scale")
standata$y_prior_scale1 <- if (M > 0) as.array(prior_scale[[1]]) else as.array(double(0))
standata$y_prior_scale2 <- if (M > 1) as.array(prior_scale[[2]]) else as.array(double(0))
standata$y_prior_scale3 <- if (M > 2) as.array(prior_scale[[3]]) else as.array(double(0))

prior_df <- fetch(y_prior_stuff, "prior_df")
standata$y_prior_df1 <- if (M > 0) prior_df[[1]] else as.array(double(0))
standata$y_prior_df2 <- if (M > 1) prior_df[[2]] else as.array(double(0))
standata$y_prior_df3 <- if (M > 2) prior_df[[3]] else as.array(double(0))

# hs priors only
standata$y_global_prior_scale <- fetch_array(y_prior_stuff, "global_prior_scale") 
standata$y_global_prior_df <- fetch_array(y_prior_stuff, "global_prior_df")
standata$y_slab_df <- fetch_array(y_prior_stuff, "slab_df")
standata$y_slab_scale <- fetch_array(y_prior_stuff, "slab_scale")

# Priors for group specific terms
standata$t <- length(cnms)
standata$p <- as.array(sapply(cnms, length))
standata$l <- as.array(
  sapply(cnms_nms, FUN = function(nm) {
    ngrps <- unique(fetch_(y_mod, "z", "ngrps", nm))
    ngrps + 1L # add padding for _NEW_ group
  }))
standata$q <- sum(standata$p * standata$l)

if (prior_covariance$dist == "decov") {
  
  # data for decov prior
  standata$prior_dist_for_cov <- b_prior_stuff$prior_dist
  standata$b_prior_shape <- b_prior_stuff$prior_shape
  standata$b_prior_scale <- b_prior_stuff$prior_scale
  standata$b_prior_concentration <- b_prior_stuff$prior_concentration
  standata$b_prior_regularization <- b_prior_stuff$prior_regularization
  standata$len_concentration <- length(standata$b_prior_concentration)
  standata$len_regularization <- length(standata$b_prior_regularization)
  standata$len_theta_L <- sum(choose(standata$p, 2), standata$p)
  
  # pass empty lkj data
  standata$b1_prior_scale <- as.array(rep(0L, standata$bK1))
  standata$b2_prior_scale <- as.array(rep(0L, standata$bK2))
  standata$b1_prior_df <- as.array(rep(0L, standata$bK1))
  standata$b2_prior_df <- as.array(rep(0L, standata$bK2))
  standata$b1_prior_regularization <- 1.0
  standata$b2_prior_regularization <- 1.0   
  
} else if (prior_covariance$dist == "lkj") {
  
  # data for lkj prior
  b1_prior_stuff <- b_prior_stuff[[b1_varname]]
  b1_prior_dist <- fetch_(b1_prior_stuff, "prior_dist")
  b1_prior_scale <- fetch_array(b1_prior_stuff, "prior_scale")
  b1_prior_df <- fetch_array(b1_prior_stuff, "prior_df")
  b1_prior_regularization <- fetch_(b1_prior_stuff, "prior_regularization")
  if (n_distinct(b1_prior_dist) > 1L)
    stop2("Bug found: covariance prior should be the same for all submodels.")
  if (n_distinct(b1_prior_regularization) > 1L) {
    stop2("Bug found: prior_regularization should be the same for all submodels.")
  }
  standata$prior_dist_for_cov <- unique(b1_prior_dist)
  standata$b1_prior_scale <- b1_prior_scale
  standata$b1_prior_df <- b1_prior_df
  standata$b1_prior_regularization <- if (length(b1_prior_regularization))
    unique(b1_prior_regularization) else 1.0
  
  if (standata$bK2 > 0) {
    # model has a second grouping factor
    b2_prior_stuff <- b_prior_stuff[[b2_varname]]
    b2_prior_scale <- fetch_array(b2_prior_stuff, "prior_scale")
    b2_prior_df    <- fetch_array(b2_prior_stuff, "prior_df")
    b2_prior_regularization <- fetch_(b2_prior_stuff, "prior_regularization")
    standata$b2_prior_scale <- b2_prior_scale
    standata$b2_prior_df    <- b2_prior_df
    standata$b2_prior_regularization <- unique(b2_prior_regularization)
  } else {
    # model does not have a second grouping factor
    standata$b2_prior_scale <- as.array(double(0))
    standata$b2_prior_df <- as.array(double(0))
    standata$b2_prior_regularization <- 1.0
  }
  
  # pass empty decov data
  standata$len_theta_L <- 0L
  standata$b_prior_shape <- as.array(rep(0L, standata$t))
  standata$b_prior_scale <- as.array(rep(0L, standata$t))
  standata$len_concentration <- 0L
  standata$len_regularization <- 0L
  standata$b_prior_concentration <- as.array(rep(0L, standata$len_concentration))
  standata$b_prior_regularization <- as.array(rep(0L, standata$len_regularization))   
}

# Names for longitudinal submodel parameters
y_intercept_nms <- uapply(1:M, function(m) {
  if (y_mod[[m]]$intercept_type$number > 0) 
    paste0(stub, m, "|(Intercept)") else NULL
})
y_beta_nms <- uapply(1:M, function(m) {
  if (!is.null(colnames(X[[m]]))) 
    paste0(stub, m, "|", colnames(X[[m]])) else NULL
})
y_aux_nms <- uapply(1:M, function(m) {
  famname_m <- family[[m]]$family
  if (is.gaussian(famname_m)) paste0(stub, m,"|sigma") else
    if (is.gamma(famname_m)) paste0(stub, m,"|shape") else
      if (is.ig(famname_m)) paste0(stub, m,"|lambda") else
        if (is.nb(famname_m)) paste0(stub, m,"|reciprocal_dispersion") else NULL
})        

# Names for group specific coefficients ("b pars")
b_nms <- uapply(seq_along(cnms), FUN = function(i) {
  nm <- cnms_nms[i]
  nms_i <- paste(cnms[[i]], nm)
  flevels[[nm]] <- c(gsub(" ", "_", flevels[[nm]]),
                     paste0("_NEW_", nm))
  if (length(nms_i) == 1) {
    paste0(nms_i, ":", flevels[[nm]])
  } else {
    c(t(sapply(nms_i, paste0, ":", flevels[[nm]])))
  }
})

# Names for Sigma matrix
Sigma_nms <- get_Sigma_nms(cnms)

#----------------
# Event submodel
#----------------

if (is_jm) { # begin jm block
} # end jm block

#---------------
# Prior summary
#---------------

prior_info <- summarize_jm_prior(
  user_priorLong = y_user_prior_stuff,
  user_priorLong_intercept = y_user_prior_intercept_stuff,
  user_priorLong_aux = y_user_prior_aux_stuff,
  if (is_jm) user_priorEvent = e_user_prior_stuff,
  if (is_jm) user_priorEvent_intercept = e_user_prior_intercept_stuff,
  if (is_jm) user_priorEvent_aux = e_user_prior_aux_stuff,
  if (is_jm) user_priorEvent_assoc = e_user_prior_assoc_stuff,
  user_prior_covariance = prior_covariance,
  b_user_prior_stuff = b_user_prior_stuff,
  b_prior_stuff = b_prior_stuff,
  y_has_intercept = fetch_(y_mod, "x", "has_intercept"),
  y_has_predictors = fetch_(y_mod, "x", "K") > 0,
  if (is_jm) e_has_intercept = standata$e_has_intercept,
  if (is_jm) e_has_predictors = standata$e_K > 0,
  if (is_jm) has_assoc = a_K > 0,
  adjusted_priorLong_scale = fetch(y_prior_stuff, "prior_scale"),
  adjusted_priorLong_intercept_scale = fetch(y_prior_intercept_stuff, "prior_scale"),
  adjusted_priorLong_aux_scale = fetch(y_prior_aux_stuff, "prior_scale"),
  if (is_jm) adjusted_priorEvent_scale = e_prior_stuff$prior_scale,
  if (is_jm) adjusted_priorEvent_intercept_scale = e_prior_intercept_stuff$prior_scale,
  if (is_jm) adjusted_priorEvent_aux_scale = e_prior_aux_stuff$prior_scale,
  if (is_jm) adjusted_priorEvent_assoc_scale = e_prior_assoc_stuff$prior_scale,
  family = family, 
  if (is_jm) basehaz = basehaz,
  stub_for_names = if (is_jm) "Long" else "y"
)  

#-----------
# Fit model
#-----------

# # call stan() to draw from posterior distribution
# stanfit <- if (is_jm) stanmodels$jm else stanmodels$mvmer
# pars <- pars_to_monitor(standata, is_jm = is_jm)
# if (M == 1L) 
#   cat("Fitting a univariate", if (is_jm) "joint" else "glmer", "model.\n\n")
# if (M  > 1L) 
#   cat("Fitting a multivariate", if (is_jm) "joint" else "glmer", "model.\n\n")
# 
# if (algorithm == "sampling") {
#   cat("Please note the warmup may be much slower than later iterations!\n")             
#   sampling_args <- set_jm_sampling_args(
#     object = stanfit,
#     cnms = cnms,
#     # user_dots = list(...), 
#     user_dots = list(), 
#     user_adapt_delta = adapt_delta,
#     user_max_treedepth = max_treedepth,
#     data = standata, 
#     pars = pars, 
#     init = init,
#     show_messages = FALSE)
#   stanfit <- do.call(sampling, sampling_args)
# } else {
#   # meanfield or fullrank vb
#   # stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
#   #                      algorithm = algorithm, ...)  
#   stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
#                        algorithm = algorithm)  
# }
# check <- check_stanfit(stanfit)
# if (!isTRUE(check)) return(standata)
# 
# # Sigma values in stanmat
# if (prior_covariance$dist == "decov" && standata$len_theta_L)
#   stanfit <- evaluate_Sigma(stanfit, cnms)
# 
# if (is_jm) { # begin jm block
#   
#   e_intercept_nms <- "Event|(Intercept)"
#   e_beta_nms <- if (e_mod$K) paste0("Event|", colnames(e_mod$Xq)) else NULL  
#   e_aux_nms <- 
#     if (basehaz$type_name == "weibull") "Event|weibull-shape" else 
#       if (basehaz$type_name == "bs") paste0("Event|b-splines-coef", seq(basehaz$df)) else
#         if (basehaz$type_name == "piecewise") paste0("Event|piecewise-coef", seq(basehaz$df)) 
#   e_assoc_nms <- character()  
#   for (m in 1:M) {
#     if (assoc["etavalue",         ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue"))
#     if (assoc["etavalue_data",    ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue:", colnames(a_mod[[m]][["X_data"]][["etavalue_data"]])))
#     if (assoc["etavalue_etavalue",][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue:Long", assoc["which_interactions",][[m]][["etavalue_etavalue"]], "|etavalue"))
#     if (assoc["etavalue_muvalue", ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue:Long", assoc["which_interactions",][[m]][["etavalue_muvalue"]],  "|muvalue"))
#     if (assoc["etaslope",         ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etaslope"))
#     if (assoc["etaslope_data",    ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etaslope:", colnames(a_mod[[m]][["X_data"]][["etaslope_data"]])))    
#     if (assoc["etaauc",           ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etaauc"))
#     if (assoc["muvalue",          ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue"))
#     if (assoc["muvalue_data",     ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue:", colnames(a_mod[[m]][["X_data"]][["muvalue_data"]])))    
#     if (assoc["muvalue_etavalue", ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_etavalue"]], "|etavalue"))
#     if (assoc["muvalue_muvalue",  ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_muvalue"]],  "|muvalue"))
#     if (assoc["muslope",          ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muslope"))
#     if (assoc["muslope_data",     ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muslope:", colnames(a_mod[[m]][["X_data"]][["muslope_data"]])))    
#     if (assoc["muauc",            ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muauc"))
#   }
#   if (sum(standata$size_which_b)) {
#     temp_g_nms <- lapply(1:M, FUN = function(m) {
#       all_nms <- paste0(paste0("Long", m, "|b["), y_mod[[m]]$z$group_cnms[[id_var]], "]")
#       all_nms[assoc["which_b_zindex",][[m]]]})
#     e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|", unlist(temp_g_nms)))
#   }
#   if (sum(standata$size_which_coef)) {
#     temp_g_nms <- lapply(1:M, FUN = function(m) {
#       all_nms <- paste0(paste0("Long", m, "|coef["), y_mod[[m]]$z$group_cnms[[id_var]], "]")
#       all_nms[assoc["which_coef_zindex",][[m]]]})
#     e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|", unlist(temp_g_nms)))
#   }
#   
# } # end jm block
# 
# new_names <- c(y_intercept_nms,
#                y_beta_nms,
#                if (is_jm) e_intercept_nms,
#                if (is_jm) e_beta_nms,
#                if (is_jm) e_assoc_nms,                   
#                if (length(standata$q)) c(paste0("b[", b_nms, "]")),
#                y_aux_nms,
#                if (is_jm) e_aux_nms,
#                paste0("Sigma[", Sigma_nms, "]"),
#                paste0(stub, 1:M, "|mean_PPD"), 
#                "log-posterior")
# stanfit@sim$fnames_oi <- new_names
# 
# stanfit_str <- nlist(.Data = stanfit, prior_info, y_mod, cnms, flevels)
# if (is_jm)
#   stanfit_str <- c(stanfit_str, nlist(e_mod, a_mod, assoc, basehaz, 
#                                       id_var, grp_stuff, scale_assoc))
# 
# do.call("structure", stanfit_str)






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OK WE HAVE STAN DATA AND THE STAN CODE NOW FOR THE MGLMM. LET'S SAVE THE STAN DATA
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
standata_mlmm <- standata


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NEXT, THE NEW MULTISTATE PORTION DATA PREPARATION
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# [1] LOAD DATA
# DONE, SEE ABOVE
head(simdata)

# [2] DEFINE DELTA INDICATORS
mylag <- function(x,k) c(rep(NA,k),head(x,-k))
df <- simdata %>% filter(outcome_number==1)
df$lag_year <- mylag(df$year,1)
df$lag_state <- mylag(df$state,1)
# now define "first" and "last"
df <- df %>% group_by(id) %>%
  mutate(first = row_number() == min( row_number() ),
         last = row_number() == max( row_number() )
  )
# next, define our indicator variables
df <- df %>% mutate(deltaI = case_when(last==FALSE & state==1 & lag_state==1 ~ 1,
                                       last==FALSE & state==2 & lag_state==1 ~ 2,
                                       last==FALSE & state==2 & lag_state==2 ~ 3),
                    deltaC = case_when(last==TRUE & state==1 & lag_state==1 ~ 1,
                                       last==TRUE & state==2 & lag_state==1 ~ 2,
                                       last==TRUE & state==2 & lag_state==2 ~ 3,
                                       last==TRUE & state==999 & lag_state==1 ~ 4,
                                       last==TRUE & state==999 & lag_state==2 ~ 5,
                                       last==TRUE & state==3 & lag_state==1 ~ 6,
                                       last==TRUE & state==3 & lag_state==2 ~ 7)
)
# check the resulting numbers:
table(df$deltaI)
table(df$deltaC)
# [3] REDUCE THE LINES OF DATA BY COMBINING INTERVALS FOR 1-> 1 (TYPE 1) AND 2->2 (TYPE 3) TRANSITIONS
# THIS FIRST STEP FORMS THE "FIRST" AND "LAST" VARIABLES
df2 <- df %>% group_by(id,deltaI) %>%
  mutate(first.deltaI = row_number() == min( row_number() ),
         last.deltaI = row_number() == max( row_number() )
  )
# TYPE 1 TRANSITIONS:
type1 <- df2 %>% filter(outcome_number==1 & (deltaI==1|deltaC==1)) # the way we define this automatically removes the "first visit" that contains no information
type1.filtered <- subset(type1,!((first.deltaI==FALSE & last.deltaI==FALSE)|(first.deltaI==TRUE & last.deltaI==TRUE)))
check <- type1.filtered %>% count(id) # make sure everyone has 2 observations, so that we can form 1 continuous interval
check %>% filter(n!=2)
type1.filtered$lag_year2 <- mylag(type1.filtered$lag_year,1)
type1.filtered2 <- type1.filtered %>% group_by(id) %>% slice_tail(n=1)
# now stack the "single-interval" people onto the above "multiple-interval" people
type1.filtered3 <- subset(type1,first.deltaI==TRUE & last.deltaI==TRUE) # notice lag_year is always 0
type1.filtered3$lag_year2 <- type1.filtered3$lag_year
type1.filtered4 <- bind_rows(type1.filtered2,type1.filtered3)
# check <- type1.filtered4 %>% count(id)
# check %>% filter(n!=1)

# TYPE 3 TRANSITIONS:
type3 <- df2 %>% filter(outcome_number==1 & (deltaI==3|deltaC==3))
type3.filtered <- subset(type3,!((first.deltaI==FALSE & last.deltaI==FALSE)|(first.deltaI==TRUE & last.deltaI==TRUE)))
check <- type3.filtered %>% count(id) # make sure everyone has 2 observations, so that we can form 1 continuous interval
check %>% filter(n!=2)
type3.filtered$lag_year2 <- mylag(type3.filtered$lag_year,1)
type3.filtered2 <- type3.filtered %>% group_by(id) %>% slice_tail(n=1)
# now stack the "single-interval" people onto the above "multiple-interval" people
type3.filtered3 <- subset(type3,first.deltaI==TRUE & last.deltaI==TRUE) # notice lag_year is always 0
type3.filtered3$lag_year2 <- type3.filtered3$lag_year
type3.filtered4 <- bind_rows(type3.filtered2,type3.filtered3)
# check <- type3.filtered4 %>% count(id)
# check %>% filter(n!=1)

#########################################################################
# [4] GENERATE THE DATASETS THAT WILL BE NEEDED FOR THE STAN LIST
#########################################################################
# GAUSS-KRONROD ELEMENTS - NEED TO GENERATE FIRST
x1=0.991455371120813	;
x2=0.949107912342759	;
x3=0.864864423359769	;
x4=0.741531185599394	;
x5=0.586087235467691	;
x6=0.405845151377397	;
x7=0.207784955007898	;
x8=0;
x9=-1*x1; x10=-1*x2; x11=-1*x3; x12=-1*x4; x13=-1*x5; x14=-1*x6; x15=-1*x7;
xk <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)

wt1=0.022935322010529	;
wt2=0.063092092629979	;
wt3=0.104790010322250	;
wt4=0.140653259715525	;
wt5=0.169004726639267	;
wt6=0.190350578064785	;
wt7=0.204432940075298	;
wt8=0.209482141084728	;
wt9=wt1; wt10=wt2; wt11=wt3; wt12=wt4; wt13=wt5; wt14=wt6; wt15=wt7;
wt <- c(wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,wt9,wt10,wt11,wt12,wt13,wt14,wt15)



# PWC knots required for some design matrices
pwc_knot1 <- 4.8
pwc_knot2 <- 17.4182



# TYPE 1: BASED ON THE "FILTERED" APPROACH
# GENERATE THE QUADRATURE POINTS + WEIGHTS
# USE DATASET type1.filtered4 GENERATED AT THE TOP OF THIS PROGRAM
type1_year <- type1.filtered4$year
type1_lag_year <- type1.filtered4$lag_year2
type1_filtered_qpoints <- list()
type1_filtered_wpoints <- list()
for (i in 1:length(type1.filtered4$year)) {
  vector1 <- data.frame("qpoints"=((xk+1)/2)*(type1.filtered4$year[i] - type1.filtered4$lag_year2[i]) + type1.filtered4$lag_year2[i])
  vector2 <- data.frame("wpoints"=((type1.filtered4$year[i] - type1.filtered4$lag_year2[i])/2)*wt)
  type1_filtered_qpoints[[i]] <- vector1
  type1_filtered_wpoints[[i]] <- vector2
}
type1_filtered_qpoints <- do.call("rbind",type1_filtered_qpoints)
type1_filtered_wpoints <- do.call("rbind",type1_filtered_wpoints)
# design matrix for baseline transition:
type1_filtered_qpoints$I1 <- ifelse(type1_filtered_qpoints$qpoints<=pwc_knot1,1,0)
type1_filtered_qpoints$I2 <- ifelse(pwc_knot1<type1_filtered_qpoints$qpoints & type1_filtered_qpoints$qpoints<=pwc_knot2,1,0)
type1_filtered_qpoints$I3 <- ifelse(pwc_knot2<type1_filtered_qpoints$qpoints,1,0)
# design matrix for covariates and longitudinal submodel defined below


# TYPE 2: LEAVE THE SAME FOR NOW, INVOLVES ITERATED INTEGRAL
type2 <- df %>% filter(outcome_number==1 & (deltaI==2|deltaC==2))
type2_N <- nrow(type2)
type2_year <- type2$year
type2_lag_year <- type2$lag_year


# TYPE 3: FILTERED APPROACH
# GENERATE THE QUAD POINTS + WEIGHTS
# USE DATASET type3.filtered4 GENERATED AT THE TOP OF THIS PROGRAM
type3_year <- type3.filtered4$year
type3_lag_year <- type3.filtered4$lag_year2
type3_filtered_qpoints <- list()
type3_filtered_wpoints <- list()
for (i in 1:length(type3.filtered4$year)) {
  vector1 <- data.frame("qpoints"=((xk+1)/2)*(type3.filtered4$year[i] - type3.filtered4$lag_year2[i]) + type3.filtered4$lag_year2[i])
  vector2 <- data.frame("wpoints"=((type3.filtered4$year[i] - type3.filtered4$lag_year2[i])/2)*wt)
  type3_filtered_qpoints[[i]] <- vector1
  type3_filtered_wpoints[[i]] <- vector2
}
type3_filtered_qpoints <- do.call("rbind",type3_filtered_qpoints)
type3_filtered_wpoints <- do.call("rbind",type3_filtered_wpoints)
# design matrix for baseline transition:
type3_filtered_qpoints$I1 <- ifelse(type3_filtered_qpoints$qpoints<=pwc_knot1,1,0)
type3_filtered_qpoints$I2 <- ifelse(pwc_knot1<type3_filtered_qpoints$qpoints & type3_filtered_qpoints$qpoints<=pwc_knot2,1,0)
type3_filtered_qpoints$I3 <- ifelse(pwc_knot2<type3_filtered_qpoints$qpoints,1,0)
# design matrix for covariates and longitudinal submodel defined below


# TYPE 4: LEAVE THE SAME FOR NOW, INVOLVES ITERATED INTEGRAL
type4 <- df %>% filter(outcome_number==1 & deltaC==4)
type4_N <- nrow(type4)
type4_year <- type4$year
type4_lag_year <- type4$lag_year


# TYPE 5: 
# GENERATE THE QUAD POINTS + WEIGHTS
type5 <- df %>% filter(outcome_number==1 & deltaC==5)
# type5_N <- nrow(type5)
type5_year <- type5$year
type5_lag_year <- type5$lag_year
type5_qpoints <- list()
type5_wpoints <- list()
for (i in 1:length(type5$year)) {
  vector1 <- data.frame("qpoints"=((xk+1)/2)*(type5$year[i] - type5$lag_year[i]) + type5$lag_year[i]) # NOT LAG_YEAR2 ANYMORE!
  vector2 <- data.frame("wpoints"=((type5$year[i] - type5$lag_year[i])/2)*wt) # NOT LAG_YEAR2 ANYMORE!
  type5_qpoints[[i]] <- vector1
  type5_wpoints[[i]] <- vector2
}
type5_qpoints <- do.call("rbind",type5_qpoints)
type5_wpoints <- do.call("rbind",type5_wpoints)
# design matrix for baseline transition:
type5_qpoints$I1 <- ifelse(type5_qpoints$qpoints<=pwc_knot1,1,0)
type5_qpoints$I2 <- ifelse(pwc_knot1<type5_qpoints$qpoints & type5_qpoints$qpoints<=pwc_knot2,1,0)
type5_qpoints$I3 <- ifelse(pwc_knot2<type5_qpoints$qpoints,1,0)
# design matrix for covariates defined below


# TYPE 6, PROBABILITIES ONLY: LEAVE THE SAME FOR NOW, INVOLVES ITERATED INTEGRAL
type6 <- df %>% filter(outcome_number==1 & deltaC==6)
type6_N <- nrow(type6)
type6_year <- type6$year
type6_lag_year <- type6$lag_year


# TYPE 7: PROBABILITIES ONLY
# GENERATE THE QUAD POINTS + WEIGHTS
type7 <- df %>% filter(outcome_number==1 & deltaC==7)
type7_year <- type7$year
type7_lag_year <- type7$lag_year
type7_qpoints <- list()
type7_wpoints <- list()
for (i in 1:length(type7$year)) {
  vector1 <- data.frame("qpoints"=((xk+1)/2)*(type7$year[i] - type7$lag_year[i]) + type7$lag_year[i]) # NOT LAG_YEAR2 ANYMORE!
  vector2 <- data.frame("wpoints"=((type7$year[i] - type7$lag_year[i])/2)*wt) # NOT LAG_YEAR2 ANYMORE!
  type7_qpoints[[i]] <- vector1
  type7_wpoints[[i]] <- vector2
}
type7_qpoints <- do.call("rbind",type7_qpoints)
type7_wpoints <- do.call("rbind",type7_wpoints)
# design matrix for baseline transition:
type7_qpoints$I1 <- ifelse(type7_qpoints$qpoints<=pwc_knot1,1,0)
type7_qpoints$I2 <- ifelse(pwc_knot1<type7_qpoints$qpoints & type7_qpoints$qpoints<=pwc_knot2,1,0)
type7_qpoints$I3 <- ifelse(pwc_knot2<type7_qpoints$qpoints,1,0)
# design matrix for covariates defined below



# # TYPE 6 AND TYPE 7: EVENTS ONLY, NO NEED FOR QUADRATURE
events <- df %>% filter(outcome_number==1 & (deltaC==6|deltaC==7))
events$I1 <- ifelse(events$year<=pwc_knot1,1,0)
events$I2 <- ifelse(pwc_knot1<events$year & events$year<=pwc_knot2,1,0)
events$I3 <- ifelse(pwc_knot2<events$year,1,0)


################################# WORK FOR THE XBAR VARIABLES
for_xbar <- df.long %>% group_by(id) %>% slice_head(n=1)






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOW BUILD THE STANDATA
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
standata_paper3 <- standata_mlmm



# FIRST, SOME PARAMETER STUFF
standata_paper3$e_K = as.array(c(2,1))
# NOW THE PWC PARAMETERS
standata_paper3$basehaz_df = 3
# x112 = ,
# x123 = ,
standata_paper3$a12_K = 3 # association parameters
standata_paper3$a23_K = 3


### TYPE 1 DATA - FULL MATRIX CALCULATION
standata_paper3$type1_N_GK = nrow(type1_filtered_qpoints)
standata_paper3$type1_id = rep(type1.filtered4$id,each=15)
standata_paper3$type1_qpoints = type1_filtered_qpoints$qpoints
standata_paper3$type1_wpoints = type1_filtered_wpoints$wpoints
standata_paper3$nrow_e_Xq_type1 = nrow(type1_filtered_qpoints) # NUMBER OF ROWS FOR COVARIATES FOR LINEAR PREDICTOR
standata_paper3$e_Xq_type1 = cbind(rep(type1.filtered4$X1cont,each=15),rep(type1.filtered4$X3cat,each=15)) # COVARIATES FOR LINEAR PREDICTOR
standata_paper3$basehaz_X_type1 = cbind(type1_filtered_qpoints$I1,type1_filtered_qpoints$I2,type1_filtered_qpoints$I3) # BASELINE TRANSITION DATA

standata_paper3$nrow_y_Xq_type1 = nrow(type1_filtered_qpoints)  # NUMBER OF ROWS TO EVALUATE THE LONGITUDINAL OUTCOMES

standata_paper3$y_Xq_type1 = cbind(#rep(1,each=nrow(type1_filtered_qpoints)), # FE DESIGN - # NO INTERCEPT
  type1_filtered_qpoints$qpoints,
  rep(type1.filtered4$X1cont,each=15),
  rep(type1.filtered4$X2cat,each=15)
)
standata_paper3$y_Zq_type1 = t(cbind(rep(1,each=nrow(type1_filtered_qpoints)), # RE DESIGN, NOTE HAS TO BE TRANSPOSED
                                     type1_filtered_qpoints$qpoints))
standata_paper3$y_Zq_id_type1 = rep(type1.filtered4$id,each=15) # ID's FOR GROUP INDEXING FOR RE DESIGN





# TYPE 3 DATA - FULL MATRIX CALCULATION
standata_paper3$type3_N_GK = nrow(type3_filtered_qpoints)
standata_paper3$type3_id = rep(type3.filtered4$id,each=15)
standata_paper3$type3_qpoints = type3_filtered_qpoints$qpoints
standata_paper3$type3_wpoints = type3_filtered_wpoints$wpoints
standata_paper3$nrow_e_Xq_type3 = nrow(type3_filtered_qpoints) # covariates
standata_paper3$e_Xq_type3 = cbind(rep(type3.filtered4$X1cont,each=15)) # covariates, only X1
standata_paper3$basehaz_X_type3 = cbind(type3_filtered_qpoints$I1,type3_filtered_qpoints$I2,type3_filtered_qpoints$I3) # baseline transition

standata_paper3$nrow_y_Xq_type3 = nrow(type3_filtered_qpoints)  # NUMBER OF ROWS TO EVALUATE THE LONGITUDINAL OUTCOME

standata_paper3$y_Xq_type3 = cbind(#rep(1,each=nrow(type3_filtered_qpoints)), # FE DESIGN # NO INTERCEPT
  type3_filtered_qpoints$qpoints,
  rep(type3.filtered4$X1cont,each=15),
  rep(type3.filtered4$X2cat,each=15)
)
standata_paper3$y_Zq_type3 = t(cbind(rep(1,each=nrow(type3_filtered_qpoints)), # RE DESIGN, NOTE HAS TO BE TRANSPOSED
                                     type3_filtered_qpoints$qpoints))
standata_paper3$y_Zq_id_type3 = rep(type3.filtered4$id,each=15) # ID's FOR GROUP INDEXING FOR RE DESIGN





# TYPE 5 DATA - FULL MATRIX CALCULATION
standata_paper3$type5_N_GK = nrow(type5_qpoints)
standata_paper3$type5_id = rep(type5$id,each=15)
standata_paper3$type5_qpoints = type5_qpoints$qpoints
standata_paper3$type5_wpoints = type5_wpoints$wpoints
standata_paper3$nrow_e_Xq_type5 = nrow(type5_qpoints) # covariates
standata_paper3$e_Xq_type5 = cbind(rep(type5$X1cont,each=15))# covariates, only X1
standata_paper3$basehaz_X_type5 = cbind(type5_qpoints$I1,type5_qpoints$I2,type5_qpoints$I3) # baseline transition

standata_paper3$nrow_y_Xq_type5 = nrow(type5_qpoints)  # NUMBER OF ROWS TO EVALUATE THE LONGITUDINAL OUTCOME

standata_paper3$y_Xq_type5 = cbind(#rep(1,each=nrow(type5_qpoints)), # FE DESIGN # NO INTERCEPT
  type5_qpoints$qpoints,
  rep(type5$X1cont,each=15),
  rep(type5$X2cat,each=15)
)
standata_paper3$y_Zq_type5 = t(cbind(rep(1,each=nrow(type5_qpoints)), # RE DESIGN, NOTE HAS TO BE TRANSPOSED
                                     type5_qpoints$qpoints))
standata_paper3$y_Zq_id_type5 = rep(type5$id,each=15) # ID's FOR GROUP INDEXING FOR RE DESIGN



# TYPE 7 DATA - FULL MATRIX CALCULATION
standata_paper3$type7_N_GK = nrow(type7_qpoints)
standata_paper3$type7_id = rep(type7$id,each=15)
standata_paper3$type7_qpoints = type7_qpoints$qpoints
standata_paper3$type7_wpoints = type7_wpoints$wpoints
standata_paper3$nrow_e_Xq_type7 = nrow(type7_qpoints) # covariates
standata_paper3$e_Xq_type7 = cbind(rep(type7$X1cont,each=15)) # covariates, only X1
standata_paper3$basehaz_X_type7 = cbind(type7_qpoints$I1,type7_qpoints$I2,type7_qpoints$I3) # baseline transition

standata_paper3$nrow_y_Xq_type7 = nrow(type7_qpoints)  # NUMBER OF ROWS TO EVALUATE THE LONGITUDINAL OUTCOME

standata_paper3$y_Xq_type7 = cbind(# rep(1,each=nrow(type7_qpoints)), # FE DESIGN # NO INTERCEPT
  type7_qpoints$qpoints,
  rep(type7$X1cont,each=15),
  rep(type7$X2cat,each=15)
)
standata_paper3$y_Zq_type7 = t(cbind(rep(1,each=nrow(type7_qpoints)), # RE DESIGN, NOTE HAS TO BE TRANSPOSED
                                     type7_qpoints$qpoints))
standata_paper3$y_Zq_id_type7 = rep(type7$id,each=15) # ID's FOR GROUP INDEXING FOR RE DESIGN


# EVENTS ALONE - FULL MATRIX CALCULATION - A LITTLE DIFFERENT FORMAT THAN ABOVE SINCE NO QUADRATURE REQUIRED
standata_paper3$events_N = nrow(events)
standata_paper3$events_year = events$year
standata_paper3$events_id = events$id
standata_paper3$nrow_e_Xq_events = nrow(events) # covariates
standata_paper3$e_Xq_events = cbind(events$X1cont) # covariates, only X1
standata_paper3$basehaz_X_events = cbind(events$I1,events$I2,events$I3) # baseline transition

standata_paper3$nrow_y_Xq_events = nrow(events)  # NUMBER OF ROWS TO EVALUATE THE LONGITUDINAL OUTCOME

standata_paper3$y_Xq_events = cbind(#rep(1,each=nrow(events)), # FE DESIGN # A LITTLE DIFFERENT FORMAT SINCE NO QUADRATURE REQUIRED # NO INTERCEPT
  events$year,
  events$X1cont,
  events$X2cat
)
standata_paper3$y_Zq_events = t(cbind(rep(1,each=nrow(events)), # RE DESIGN, NOTE HAS TO BE TRANSPOSED
                                      events$year))
standata_paper3$y_Zq_id_events = events$id # ID's FOR GROUP INDEXING FOR RE DESIGN




# TYPE 2: THE "VECTORIZED" ITERATED INTEGRAL APPROACH
standata_paper3$type2_N = type2_N
standata_paper3$type2_year = type2_year
standata_paper3$type2_lag_year = type2_lag_year
standata_paper3$type2_id = type2$id
standata_paper3$type2_X1 = type2$X1cont
standata_paper3$type2_X3 = type2$X3cat
standata_paper3$type2_X2 = type2$X2cat





# TYPE 4: THE "VECTORIZED" ITERATED INTEGRAL APPROACH
standata_paper3$type4_N = type4_N
standata_paper3$type4_year = type4_year
standata_paper3$type4_lag_year = type4_lag_year
standata_paper3$type4_id = type4$id
standata_paper3$type4_X1 = type4$X1cont
standata_paper3$type4_X3 = type4$X3cat
standata_paper3$type4_X2 = type4$X2cat



# TYPE 6: THE "VECTORIZED" ITERATED INTEGRAL APPROACH
standata_paper3$type6_N = type6_N
standata_paper3$type6_year = type6_year
standata_paper3$type6_lag_year = type6_lag_year
standata_paper3$type6_id = type6$id
standata_paper3$type6_X1 = type6$X1cont
standata_paper3$type6_X3 = type6$X3cat
standata_paper3$type6_X2 = type6$X2cat





# LOAD THIS MODEL AS REFERENCE - WE NEED SOME DUMMY VARIABLE FOR THE RSTANARM STAN FUNCTION (TO BE USED LATER)
standata_tune_jm <- readRDS('M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/AUXILLARY FILES/standata_tune_jm.rds') 
standata_paper3$y1_offset_eta <- standata_tune_jm$y1_offset_eta
standata_paper3$y2_offset_eta <- standata_tune_jm$y2_offset_eta
standata_paper3$y3_offset_eta <- standata_tune_jm$y3_offset_eta
standata_paper3$y1_z2q_eta <- matrix(,nrow=0,ncol = 0)
standata_paper3$y2_z2q_eta <- matrix(,nrow=0,ncol = 0)
standata_paper3$y3_z2q_eta <- matrix(,nrow=0,ncol = 0)
standata_paper3$y1_z2q_id_eta <- standata_tune_jm$y1_z2q_id_eta
standata_paper3$y2_z2q_id_eta <- standata_tune_jm$y2_z2q_id_eta
standata_paper3$y3_z2q_id_eta <- standata_tune_jm$y3_z2q_id_eta



# the extra hyperparameters
standata_paper3$e12_prior_mean = c(0,0)
standata_paper3$e23_prior_mean = as.array(0)
standata_paper3$a12_prior_mean = c(0,0,0)
standata_paper3$a23_prior_mean = c(0,0,0)
standata_paper3$e12_prior_mean_for_aux = c(0,0,0)
standata_paper3$e23_prior_mean_for_aux = c(0,0,0)


standata_paper3$e12_prior_scale = c(2.5,2.5)
standata_paper3$e23_prior_scale = as.array(2.5)
standata_paper3$a12_prior_scale = c(2.5,2.5,2.5)
standata_paper3$a23_prior_scale = c(2.5,2.5,2.5)
standata_paper3$e12_prior_scale_for_aux = c(20,20,20) 
standata_paper3$e23_prior_scale_for_aux = c(20,20,20)


# and lastly:
standata_paper3$xk <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)
standata_paper3$wt <- c(wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,wt9,wt10,wt11,wt12,wt13,wt14,wt15)
standata_paper3$pwc_knot1 <- 4.8
standata_paper3$pwc_knot2 <- 17.4182



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# INITIAL VALUES FOR STAN
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df.long <- simdata %>% filter(outcome_number==2)
lme.y1 <- lme(y1 ~ year + X1cont + X2cat, random =~ 1 + year | id, data = df.long,
              control = c(opt = 'optim'))
lme.y2 <- lme(y2 ~ year + X1cont + X2cat, random =~ 1 + year | id, data = df.long,
              control = c(opt = 'optim'))
lme.y3 <- lme(y3 ~ year + X1cont + X2cat, random =~ 1 + year | id, data = df.long,
              control = c(opt = 'optim'))
summary(lme.y1)
summary(lme.y2)
summary(lme.y3)
# RANDOM EFFECTS INITIAL VALUES - THESE ARE TO BE USED TO GENERATE THE "UNSTANDARDIZED" INITIAL VALUES
matrix <- t(cbind(lme.y1$coefficients$random$id[,1],lme.y1$coefficients$random$id[,2],
                  lme.y2$coefficients$random$id[,1],lme.y2$coefficients$random$id[,2],
                  lme.y3$coefficients$random$id[,1],lme.y3$coefficients$random$id[,2]))
getVarCov(lme.y1)
getVarCov(lme.y2)
getVarCov(lme.y3)

D11 = 1.4017*1.4017 
D21 = -0.04346 
D22 = 0.2248*0.2248
D31 = 0.2565
D32 = 0.01517
D33 = 0.6314*0.6314
D41 = -0.01313
D42 = 0.008311
D43 = -0.00639
D44 = 0.07677*0.07677
D51 = -1.8694
D52 = 0.7317
D53 = 1.2070
D54 = 0.2288
D55 = 11.4293*11.4293
D61 = -0.00610
D62 = -0.09324
D63 = -0.05309
D64 = -0.05461
D65 = -4.2248
D66 = 0.9635*0.9635
D <- matrix(c(D11, D21, D31, D41, D51, D61,# note the definition is symmetric in the lower diagonal
              D21, D22, D32, D42, D52, D62,
              D31, D32, D33, D43, D53, D63,
              D41, D42, D43, D44, D54, D64, 
              D51, D52, D53, D54, D55, D65,
              D61, D62, D63, D64, D65, D66), ncol=6, nrow=6) # THIS IS THE VAR-COVAR, USING THE TRUE MODEL PARAMETERS
D

cov2cor(D) # THIS IS THE REQUIRED CORRELATION
C2 <- cov2cor(D)
chol(C2)
cholC <- t(chol(C2)) # the "true" cholesky matrix needs to be lower-tri, but the default in R is upper
S <- diag(c(sqrt(D11),sqrt(D22),sqrt(D33),sqrt(D44),sqrt(D55),sqrt(D66)), nrow = 6, ncol = 6)
solve(S) # the inverse
solve(cholC) # the inverse
z_b_mat_init <- solve(cholC) %*% solve(S) %*% matrix # NOTE THAT WE COULD JUST USE THE TRUE SIMULATED RANDOM EFFECTS VALUES STORED IN B, 
# BUT WHEN WE RUN THE MODELS, WE DO NOT ATCUALLY STORE THEM


inits_paper3 <- list(
  z_yBeta1 = c(0.19/standata_paper3$y_prior_scale1[1], 0.1527/standata_paper3$y_prior_scale1[2], -2.3148/standata_paper3$y_prior_scale1[3]),     
  z_yBeta2 = c(0.02547/standata_paper3$y_prior_scale2[1], 0.05528/standata_paper3$y_prior_scale2[2], -0.3228/standata_paper3$y_prior_scale2[3]),
  z_yBeta3 = c(-1.4380/standata_paper3$y_prior_scale3[1], 0.4161/standata_paper3$y_prior_scale3[2], 1.2858/standata_paper3$y_prior_scale3[3]),
  yGamma1 = as.array(0.39),      
  yGamma2 = as.array(2.0727),
  yGamma3 = as.array(120.27),
  yAux1_unscaled = as.array(1.1362/(sd(df.long$y1)*5)), # the value divided by the scale
  yAux2_unscaled = as.array(0.6687/(sd(df.long$y2)*5)),
  yAux3_unscaled = as.array(8.3159/(sd(df.long$y3)*5)),
  e12_z_beta = c(0.07/2.5,-0.27/2.5),
  e23_z_beta = as.array(0.10/2.5),
  e12_aux_unscaled = c(-2.8/10,-3.6/10,-3.0/10),
  e23_aux_unscaled = c(-4.5/10,-5.2/10,-5.7/10),
  bSd1 = c(sqrt(D11),sqrt(D22),sqrt(D33),sqrt(D44),sqrt(D55),sqrt(D66)),
  z_bMat1 = cbind(z_b_mat_init,c(0,0)), # pad the extra value
  bCholesky1 = cholC
)

# ######################################################
# AND NOW SAVE THE OUTPUTS (OPTIONAL SINCE THEY ARE ALREADY INCLUDED IN THE FILES)
# ######################################################
setwd("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5")
saveRDS(inits_paper3,'inits_paper3.rds')
saveRDS(standata_paper3,'standata_paper3.rds')
