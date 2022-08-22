library(tidyverse)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FOR THIS EXAMPLE, WE WILL USE A SINGLE SIMULATED DATASET FOR DGP 1
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setwd("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER4")
sasdata <- readRDS('sasdata_DGP1_n1.rds')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RUN THE JMSTATEEXACT MODEL
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(JM) # JM PACKAGE REQUIRED
# TO USE THE JMSTATEMODEL FUNCTION, WE DOWNLOAD THE PACKAGE AND MANUALLY LOAD THE FUNCTIONS
# https://github.com/LoicFerrer/JMstateModel/
# LINES 123, 141, 162, AND 185 OF THE initial.surv.R PROGRAM HAD TO BE ALTERED
my_path <- c("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER4/JMstateModel-master-THESISCOPY/")
source_files <- list.files(my_path, "*.R$")  # locate all .R files
map(paste0(my_path, source_files), source)  # source all your R scripts

# nsim <- 200 # ALTER WHEN RUNNING THE LARGER SIMULATION
nsim <- 1
simlist.JMstate.coefs <- list()
simlist.JMstate.lrt <- list()
simlist.JMstate.mix <- list()
for (i in 1:nsim){
  data_long <- sasdata %>% filter(sim_number==i & outcome_number==2)
  
  data_mstate1 <- sasdata %>% filter(sim_number==i & outcome_number==1) %>%
    dplyr::group_by(id) %>%
    mutate(first_3 = state == 3 & !duplicated(state == 3))  %>%
    mutate(first_2 = state == 2 & !duplicated(state == 2))  %>%
    arrange(id)
  
  times2 <- data_mstate1 %>% dplyr::group_by(id) %>%
    filter(first_2==TRUE) %>%
    rename(time2 = year) %>%
    dplyr::select(id,time2)
  
  times3 <- data_mstate1 %>% dplyr::group_by(id) %>%
    filter(first_3==TRUE) %>%
    rename(time3 = year) %>%
    dplyr::select(id,time3)
  
  last_observed1 <- sasdata %>% filter(sim_number==i & outcome_number==1 & state==1) %>%
    dplyr::group_by(id) %>%
    slice_tail(n=1) %>%
    rename(last_observed1 = year) %>%
    dplyr::select(id,last_observed1)
  
  final_time <- data_mstate1 %>% dplyr::group_by(id) %>%
    slice_tail(n=1) %>%
    rename(totaltime = year) %>%
    dplyr::select(id,totaltime, X1cont, X3cat) # COVARIATE
  
  data_mstate2 <- left_join(final_time,times2,by="id")
  data_mstate3 <- left_join(data_mstate2,times3,by="id")
  data_mstate4 <- left_join(data_mstate3,last_observed1,by="id")
  
  ms_arrange <- function (x) {
    if (is.na(x$time2) & is.na(x$time3)) {
      x_new <- data.frame('id' = rep(x$id, 2), 'from_state' = 1:2, 'to_state' = 2:3, 
                          'transition' = 1:2, 'Tstart' = rep(0, 2), 'Tstop' = x$totaltime, 'status' = rep(0, 2),
                          'X1cont'=x$X1cont, 'X3cat'=x$X3cat)
    } else {
      if (!is.na(x$time2) & is.na(x$time3)) {
        x_new <- data.frame('id' = rep(x$id, 2), 'from_state' = 1:2, 'to_state' = 2:3, 
                            'transition' = 1:2, 'Tstart' = c(0, x$time2), 
                            'Tstop' = c(x$time2, x$totaltime), 'status' = c(1, 0),
                            'X1cont'=x$X1cont, 'X3cat'=x$X3cat)
      } else {
        if (!is.na(x$time2) & !is.na(x$time3)) {
          x_new <- data.frame('id' = rep(x$id, 2), 'from_state' = 1:2, 'to_state' = 2:3, 
                              'transition' = 1:2, 'Tstart' = c(0, x$time2), 
                              'Tstop' = c(x$time2, x$time3), 'status' = c(1, 1),
                              'X1cont'=x$X1cont, 'X3cat'=x$X3cat)
        } else {
          x_new <- data.frame('id' = rep(x$id, 2), 'from_state' = 1:2, 'to_state' = 2:3, # this is the case where only 1->3 is observed
                              'transition' = 1:2, 'Tstart' = c(0, (x$time3+x$last_observed1)/2), 
                              'Tstop' = c((x$time3+x$last_observed1)/2, x$time3), 'status' = c(1, 1),
                              'X1cont'=x$X1cont, 'X3cat'=x$X3cat)
        }
      }
    } 
    
  }
  data_mstate4_split.by.id <- split(data_mstate4, data_mstate4$id)
  data_mstate4_split.by.id.2 <- lapply(data_mstate4_split.by.id, ms_arrange)
  data_mstate.2 <- do.call(rbind, data_mstate4_split.by.id.2)
  data_mstate.2$transition <- factor(data_mstate.2$transition)
  
  data_mstate.3 <- data_mstate.2 %>% filter(!(Tstop<=Tstart))
  
  # RE-DEFINE THE TRANSITION-SPECIFIC COVARIATES
  data_mstate.4 <- data_mstate.3 %>% 
    mutate(X3cat.1 = case_when(transition=="1" ~ X3cat,
                               transition=="2" ~ 0L)) %>% # need the L for integer type
    mutate(X1cont.1 = case_when(transition=="1" ~ X1cont,
                                transition=="2" ~ 0)) %>%
    mutate(X1cont.2 = case_when(transition=="1" ~ 0,
                                transition=="2" ~ X1cont))
  
  # INITIALIZE THE TWO MODELS
  lmeFit <- lme(fixed = y ~ year + X1cont, random = ~ 1 + year | id,
                data = data_long, method = "ML", control = list(opt = "optim"))
  # coxFit <- coxph(Surv(Tstart, Tstop, status) ~ (X1cont + X3cat)*strata(transition),
  #                 data = data_mstate.3, method = "breslow", x = TRUE, model = TRUE)
  coxFit <- coxph(Surv(Tstart, Tstop, status) ~ X3cat.1 + X1cont.1 + X1cont.2 + strata(transition), 
                  data = data_mstate.4, method = "breslow", x = TRUE, model = TRUE)
  
  jointFit_JMstate <-
    JMstateModel(lmeObject = lmeFit,
                 survObject = coxFit,
                 timeVar = "year",
                 parameterization = "value",
                 method = "spline-PH-aGH",
                 interFact = list(value = ~ strata(transition) - 1,
                                  #slope = ~ strata(trans) - 1,
                                  data = data_mstate.4),
                 #derivForm = dForm,
                 Mstate = TRUE,
                 data.Mstate = data_mstate.4,
                 ID.Mstate = "id",
                 control = list(GHk = 3, equal.strata.knots = FALSE, lng.in.kn = 2))
  
  summary.JMstate <- summary(jointFit_JMstate)
  VarCov <- vcov(jointFit_JMstate)
  
  mix <- data.frame('sim_number'=i,"convergence"=summary.JMstate$conv)
  lrt <- data.frame('sim_number'=i,"logLik"=summary.JMstate$logLik, "AIC"=summary.JMstate$AIC, "BIC"=summary.JMstate$BIC)
  
  oute_rand <- data.frame('sim_number'=rep(i,3),"Parameter"=c("D11","D21","D22"),
                          "Value" = c(summary.JMstate$D[1,1],summary.JMstate$D[2,1],summary.JMstate$D[2,2]),
                          "Std.Err" = sqrt(diag(VarCov[22:24, 22:24])))
  oute_sigma <- data.frame('sim_number'=i,"Parameter"="sigma", "Value" = summary.JMstate$sigma,
                           "Std.Err" = sqrt(VarCov[4, 4]))
  oute_event <- data.frame(summary.JMstate$`CoefTable-Event`)
  oute_event <- tibble::rownames_to_column(oute_event, "Parameter")
  oute_event$sim_number <- rep(i,nrow(oute_event))
  oute_long <- data.frame(summary.JMstate$`CoefTable-Long`)
  oute_long <- tibble::rownames_to_column(oute_long, "Parameter")
  oute_long$sim_number <- rep(i,nrow(oute_long))
  
  oute_coeffs <- bind_rows(oute_event,oute_long,oute_rand,oute_sigma)
  
  simlist.JMstate.coefs[[i]] <- oute_coeffs
  simlist.JMstate.lrt[[i]] <- lrt
  simlist.JMstate.mix[[i]] <- mix
}
JMstate.coefs <- do.call("rbind", simlist.JMstate.coefs)
JMstate.lrt <- do.call("rbind", simlist.JMstate.lrt)
JMstate.mix <- do.call("rbind", simlist.JMstate.mix)


# OPTIONAL CODE TO STORE ADDITIONAL RESULTS, E.G. WHEN SETTING nsim EQUAL TO 200
# setwd("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER4")
# saveRDS(JMstate.coefs,'JMstate.coefs.DGP1n200.rds')
# saveRDS(JMstate.lrt,'JMstate.lrt.DGP1n200.rds')
# saveRDS(JMstate.mix,'JMstate.mix.DGP1n200.rds')



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RUN THE MSM MODEL
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(msm)

# nsim <- 200 # ALTER WHEN RUNNING THE LARGER SIMULATION
nsim <- 1
simlist.MSM.coefs <- list()
simlist.MSM.lrt <- list()
simlist.MSM.mix <- list()
for (i in 1:nsim) {
  data_mstate <- sasdata %>% filter(sim_number==i & outcome_number==1)
  data_long <- sasdata %>% filter(sim_number==i & outcome_number==2)
  
  # CREATE THE APPROPRIATE DATASET
  final_time <- data_mstate %>% dplyr::group_by(id) %>%
    slice_tail(n=1) %>%
    rename(totaltime = year) %>%
    dplyr::select(id,totaltime, X1cont, X3cat)
  time_range <- tmerge(final_time, final_time, id=id, tstart = 0, tstop = totaltime+0.001)
  data.with.tv <- tmerge(time_range, data_long, id=id, y.tv = tdc(year, y))
  data.with.tv2 <- tmerge(data.with.tv, data_mstate, id=id, state.tv = tdc(year, state))
  state.obs <- data_mstate %>% mutate(state.obs=year) %>% dplyr::select(id, year, state.obs)
  data.with.tv3 <- left_join(data.with.tv2, state.obs, by=c("id"="id","tstart"="year"))
  data.with.tv3$state.with.censoring <- ifelse(is.na(data.with.tv3$state.obs), 999, data.with.tv3$state.tv)
  
  # RUN THE MSM MODEL
  q12=0.5
  q23=0.5
  Q <- rbind(
    c(0,q12,0),
    c(0,0,q23),
    c(0,0,0)
  )
  Q
  Q.crude  <- crudeinits.msm( state.with.censoring ~ tstart, subject=id, data = data.with.tv3, qmatrix=Q,
                              censor = 999, censor.states = c(1,2))
  msm.tv <- msm( state.with.censoring ~ tstart, subject=id, data = data.with.tv3,
                 qmatrix = Q.crude, deathexact = 3, 
                 censor = 999, censor.states = c(1,2),
                 # covariates = ~ y.tv + X1cont + X3cat, 
                 covariates = list("1-2" = ~ y.tv + X1cont + X3cat, "2-3" = ~ y.tv + X1cont),
                 center = FALSE,
                 pci = c(4.8,17.4182),
                 opt.method="optim") # control = list(trace=1,REPORT=1), 
  
  # EXTRACT AND STORE THE MODEL PARAMETERS
  oute <- data.frame("sim_number"=i,
                     "Parameter" = c("alpha12","alpha23",
                                     "lambda121","lambda122","lambda123",
                                     "lambda231","lambda232","lambda233",
                                     "gamma121","gamma122","gamma231"),     # ,"gamma232"), # NO RANDOM EFFECTS, NO LONGITUDINAL PARAMETERS
                     "Value"=c(msm.tv$Qmatrices$y.tv[1,2],msm.tv$Qmatrices$y.tv[2,3],
                               msm.tv$Qmatrices$logbaseline[1,2],
                               msm.tv$Qmatrices$logbaseline[1,2]+msm.tv$Qmatrices$`timeperiod[4.8,17.4)`[1,2],
                               msm.tv$Qmatrices$logbaseline[1,2]+msm.tv$Qmatrices$`timeperiod[17.4,Inf)`[1,2],
                               msm.tv$Qmatrices$logbaseline[2,3],
                               msm.tv$Qmatrices$logbaseline[2,3]+msm.tv$Qmatrices$`timeperiod[4.8,17.4)`[2,3],
                               msm.tv$Qmatrices$logbaseline[2,3]+msm.tv$Qmatrices$`timeperiod[17.4,Inf)`[2,3],
                               msm.tv$Qmatrices$X1cont[1,2],msm.tv$Qmatrices$X3cat[1,2],
                               msm.tv$Qmatrices$X1cont[2,3]),        #msm.tv$Qmatrices$X3cat[2,3]),
                     "SE"=c(msm.tv$QmatricesSE$y.tv[1,2],msm.tv$QmatricesSE$y.tv[2,3],
                            msm.tv$QmatricesSE$logbaseline[1,2],
                            msm.tv$QmatricesSE$`timeperiod[4.8,17.4)`[1,2],
                            msm.tv$QmatricesSE$`timeperiod[17.4,Inf)`[1,2],
                            msm.tv$QmatricesSE$logbaseline[2,3],
                            msm.tv$QmatricesSE$`timeperiod[4.8,17.4)`[2,3],
                            msm.tv$QmatricesSE$`timeperiod[17.4,Inf)`[2,3],
                            msm.tv$QmatricesSE$X1cont[1,2],msm.tv$QmatricesSE$X3cat[1,2],
                            msm.tv$QmatricesSE$X1cont[2,3])          #,msm.tv$QmatricesSE$X3cat[2,3])
  )
  mix <- data.frame('sim_number'=i,"foundse"=msm.tv$foundse)
  lrt <- data.frame('sim_number'=i,"logLik"=msm.tv$minus2loglik/-2, "AIC"=AIC(msm.tv)) 
  
  simlist.MSM.coefs[[i]] <- oute
  simlist.MSM.lrt[[i]] <- lrt
  simlist.MSM.mix[[i]] <- mix
}
MSM.coefs <- do.call("rbind", simlist.MSM.coefs)
MSM.lrt <- do.call("rbind", simlist.MSM.lrt)
MSM.mix <- do.call("rbind", simlist.MSM.mix)


# OPTIONAL CODE TO STORE ADDITIONAL RESULTS, E.G. WHEN SETTING nsim EQUAL TO 200
# setwd("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER4")
# saveRDS(MSM.coefs,'MSM.coefs.DGP1n200.rds')
# saveRDS(MSM.lrt,'MSM.lrt.DGP1n200.rds')
# saveRDS(MSM.mix,'MSM.mix.DGP1n200.rds')