library(tidyverse)
library(haven)

#-----------------------------------------------------------------------------------------------------
# We used the rsimsum package to help analyze the results of the simulation:
# Gasparini, A. (2018). Rsimsum: Summarise results from monte carlo simulation studies. Journal of Open
# Source Software, 3, 739. https://doi.org/10.21105/joss.00739
#-----------------------------------------------------------------------------------------------------
library(rsimsum)
# AT THE END, WE OUTPUT THE RESULTS FOR LATEX:
library(kableExtra)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# (THE SIMULATION RESULTS ARE STORED, NO NEED TO RE-RUN ALL 200)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setwd("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER4")
estmlong <- readRDS('estmlong_15MAY22.rds')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DEFINE THE TRUE PARAMETER VALUES
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estmlong.2 <- estmlong %>% mutate(true_value = case_when(Parameter=="alpha12" ~ 0.12,
                                                         
                                                         Parameter=="alpha23" & DGP %in% c(1,2,4,5) ~ 0.06,
                                                         Parameter=="alpha23" & DGP==3 ~ 0.09,
                                                         
                                                         Parameter=="beta0" ~ 0.39,
                                                         Parameter=="beta1" ~ 0.19,
                                                         Parameter=="beta2" ~ 0.21,
                                                         
                                                         Parameter=="d11" ~ sqrt(3.3),#Parameter=="d11" ~ 3.3,
                                                         Parameter=="d21" ~ -0.05,
                                                         Parameter=="d22" ~ sqrt(0.05),
                                                         
                                                         Parameter=="lambda121" & DGP %in% c(1,2,4,5) ~ -3.5,
                                                         Parameter=="lambda121" & DGP==3 ~ -3.4,
                                                         
                                                         Parameter=="lambda122" & DGP %in% c(1,2,4,5) ~ -5,
                                                         Parameter=="lambda122" & DGP==3 ~ -4.9,
                                                         
                                                         Parameter=="lambda123" & DGP %in% c(1,2,4,5) ~ -3.9,
                                                         Parameter=="lambda123" & DGP==3 ~ -3.8,
                                                         
                                                         Parameter=="lambda231" & DGP %in% c(1,2,4,5) ~ -5.7,
                                                         Parameter=="lambda231" & DGP==3 ~ -5.1,
                                                         
                                                         Parameter=="lambda232" & DGP %in% c(1,2,4,5) ~ -5.9,
                                                         Parameter=="lambda232" & DGP==3 ~ -5.3,
                                                         
                                                         Parameter=="lambda233" & DGP %in% c(1,2,4,5) ~ -6.5,
                                                         Parameter=="lambda233" & DGP==3 ~ -5.9,
                                                         
                                                         Parameter=="sigma" ~ 1.0977,
                                                         
                                                         Parameter=="gamma121" ~ 0.08,
                                                         Parameter=="gamma122" ~ -0.32,
                                                         Parameter=="gamma231" ~ 0.13
))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# USING RSIMSUM, SUMMARIZE ALL OF THE PARAMETERS EXCEPT FOR THE B-SPLINE 
# PARAMETERS FOR THE JMSTATEEXACT MODEL
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estmlong.3 <- subset(estmlong.2, !(Parameter %in% c("bs1(1)","bs1(2)" ,"bs2(1)","bs2(2)",
                                                    "bs3(1)","bs3(2)","bs4(1)","bs4(2)",
                                                    "bs5(1)","bs5(2)","bs6(1)","bs6(2)")))
estmlong.3$Parameter2 <- factor(estmlong.3$Parameter)
split_estmlong <- split(estmlong.3, estmlong.3$Parameter2)
simsum.results <- list()
for (z in levels(estmlong.3$Parameter2)) {
  data <- split_estmlong[[z]]
  data$Method <- factor(data$Method)
  s <- simsum(
    data = data, estvarname = "Estimate", se = "SE", true = "true_value", ref = "NLMIXED",
    methodvar = "Method", by = c("DGP"), x = TRUE
  )
  simsum.results[[z]] <- s
  print(summary(s))
}

# USING THE STORED RESULTS, WE NOW NEED TO OUTPUT THE TABLES
# MAKE SURE THE FACTORS ARE ALL IN PROPER ORDER:
levels(estmlong.3$Parameter2)
estmlong.4 <- estmlong.3
estmlong.4$Parameter2 <- factor(estmlong.4$Parameter2, levels=c("alpha12","alpha23","lambda121", "lambda122", "lambda123",
                                                                "gamma121" , "gamma122",
                                                                "lambda231", "lambda232", "lambda233","gamma231",
                                                                "beta0" ,    "beta1"  ,   "beta2","sigma",
                                                                "d11"   ,    "d21"    ,   "d22" ))
levels(estmlong.4$Parameter2)
levels(estmlong.4$Method)
levels(estmlong.4$DGP)

# NOW HAVE TO LOOP OVER THE DIFFERENT PARAMETER, DGP, AND METHOD COMBINATIONS, AND EXTRACT EACH ROW FOR THE TABLES
rows <- list()
for (z in levels(estmlong.4$Parameter2)) {
  
  simsum1 <- simsum.results[[z]]
  
  for (k in levels(estmlong.4$DGP)) {
    
    truth <- unique(simsum1$x[simsum1$x$DGP==k,]$true_value)
    # print(truth)
    simsum2 <- simsum1$summ
    simsum2$Method <- factor(simsum2$Method)
    
    for (l in levels(simsum2$Method)) {
      
      simsum3 <- subset(simsum2, DGP==k & Method==l)
      results_vector <- data.frame("Parameter"=z,
                                   "DGP"=k,
                                   "True"=round(truth,2),
                                   "Method"=l,
                                   "Mean"=round(simsum3[simsum3$stat=="thetamean",]$est,3),
                                   "Bias"=paste0(round(simsum3[simsum3$stat=="bias",]$est,3),
                                                 " ("
                                                 ,round(simsum3[simsum3$stat=="bias",]$mcse,3),
                                                 ")"),
                                   "EmpSE"=paste0(round(simsum3[simsum3$stat=="empse",]$est,3),
                                                  " ("
                                                  ,round(simsum3[simsum3$stat=="empse",]$mcse,3),
                                                  ")"),
                                   "ModSE"=paste0(round(simsum3[simsum3$stat=="modelse",]$est,3),
                                                  " ("
                                                  ,round(simsum3[simsum3$stat=="modelse",]$mcse,3),
                                                  ")"),
                                   "MSE"=paste0(round(simsum3[simsum3$stat=="mse",]$est,3),
                                                " ("
                                                ,round(simsum3[simsum3$stat=="mse",]$mcse,3),
                                                ")"),
                                   "Cover"=paste0(round(simsum3[simsum3$stat=="cover",]$est,3)#,
                                                  # " ("
                                                  # ,round(simsum3[simsum3$stat=="cover",]$mcse,3),
                                                  # ")"
                                   ))
      rows[[z]][[k]][[l]] <- results_vector
    }
  }
}
results1 <- list()
results2 <- list()
for (i in 1:length(levels(estmlong.4$Parameter2))) {
  for (j in 1:length(levels(estmlong.4$DGP))) {
    results2[[j]] <- do.call("rbind", rows[[i]][[j]])
  }
  results1[[i]] <- do.call("rbind",results2)
}
results <- do.call("rbind",results1)

# SWITCH TO WIDE FORMAT
results.2 <- results %>%
  pivot_wider(
    names_from = Method,
    values_from = c(Mean, Bias, EmpSE, ModSE, MSE, Cover)
  )

# RE-ORDER 
names(results.2)
results.3 <- results.2 %>% 
  arrange(DGP) %>%
  relocate("Parameter","DGP","True","Mean_NLMIXED",
           "Bias_NLMIXED","EmpSE_NLMIXED","ModSE_NLMIXED","MSE_NLMIXED","Cover_NLMIXED",
           "Mean_JMStateModel","Bias_JMStateModel","EmpSE_JMStateModel","ModSE_JMStateModel","MSE_JMStateModel","Cover_JMStateModel")



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOW REPEAT FOR THE BSPLINE PARAMETERS FOR JMSTATEEXACT
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estmlong.3.bspl <- subset(estmlong.2, Parameter %in% c("bs1(1)","bs1(2)" ,"bs2(1)","bs2(2)",
                                                       "bs3(1)","bs3(2)","bs4(1)","bs4(2)",
                                                       "bs5(1)","bs5(2)","bs6(1)","bs6(2)"))
estmlong.3.bspl$Parameter2 <- factor(estmlong.3.bspl$Parameter)
split_estmlong.bspl <- split(estmlong.3.bspl, estmlong.3.bspl$Parameter2)
simsum.results.bspl <- list()
for (z in levels(estmlong.3.bspl$Parameter2)) {
  data <- split_estmlong.bspl[[z]]
  data$Method <- factor(data$Method)
  s <- simsum(
    data = data, estvarname = "Estimate", se = "SE",
    methodvar = "Method", by = c("DGP"), x = TRUE
  )
  simsum.results.bspl[[z]] <- s
  #print(summary(s))
}
# USING THE STORED RESULTS, WE NOW NEED TO OUTPUT THE TABLES
# MAKE SURE THE FACTORS ARE ALL IN PROPER ORDER:
levels(estmlong.3.bspl$Parameter2)
estmlong.4.bspl <- estmlong.3.bspl
estmlong.4.bspl$Parameter2 <- factor(estmlong.4.bspl$Parameter2, levels=c("bs1(1)","bs2(1)","bs3(1)","bs4(1)", "bs5(1)","bs6(1)",
                                                                          "bs1(2)","bs2(2)","bs3(2)","bs4(2)","bs5(2)","bs6(2)" ))
levels(estmlong.4.bspl$Parameter2)
levels(estmlong.4.bspl$Method)
levels(estmlong.4.bspl$DGP)

# NOW HAVE TO LOOP OVER THE DIFFERENT PARAMETER, DGP, AND METHOD COMBINATIONS, AND EXTRACT EACH ROW FOR THE TABLES
rows.bspl <- list()
for (z in levels(estmlong.4.bspl$Parameter2)) {
  
  simsum1 <- simsum.results.bspl[[z]]
  
  for (k in levels(estmlong.4.bspl$DGP)) {
    
    # truth <- unique(simsum1$x[simsum1$x$DGP==k,]$true_value)
    # print(truth)
    simsum2 <- simsum1$summ
    simsum2$Method <- factor(simsum2$Method)
    
    for (l in levels(simsum2$Method)) {
      
      simsum3 <- subset(simsum2, DGP==k & Method==l)
      results_vector <- data.frame("Parameter"=z,
                                   "DGP"=k,
                                   "True"=NA,
                                   "Method"=l,
                                   "Mean"=round(simsum3[simsum3$stat=="thetamean",]$est,3),
                                   "Bias"=NA,
                                   "EmpSE"=paste0(round(simsum3[simsum3$stat=="empse",]$est,3),
                                                  " ("
                                                  ,round(simsum3[simsum3$stat=="empse",]$mcse,3),
                                                  ")"),
                                   "ModSE"=paste0(round(simsum3[simsum3$stat=="modelse",]$est,3),
                                                  " ("
                                                  ,round(simsum3[simsum3$stat=="modelse",]$mcse,3),
                                                  ")"),
                                   "MSE"=NA,
                                   "Cover"=NA)
      rows.bspl[[z]][[k]][[l]] <- results_vector
    }
  }
}
results1.bspl <- list()
results2.bspl <- list()
for (i in 1:length(levels(estmlong.4.bspl$Parameter2))) {
  for (j in 1:length(levels(estmlong.4.bspl$DGP))) {
    results2.bspl[[j]] <- do.call("rbind", rows.bspl[[i]][[j]])
  }
  results1.bspl[[i]] <- do.call("rbind",results2.bspl)
}
results.bspl <- do.call("rbind",results1.bspl)

# SWITCH TO WIDE FORMAT
results.2.bspl <- results.bspl %>%
  pivot_wider(
    names_from = Method,
    values_from = c(Mean, Bias, EmpSE, ModSE, MSE, Cover)
  )



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOW COMBINE EVERYTHING AND CREATE THE TABLES FOR LATEX - EXAMPLE FOR DGP 1 SHOWN
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results.latex <- bind_rows(results.3,results.2.bspl)
table(results.latex$Parameter)
results.latex2 <- results.latex %>% arrange(DGP) %>%
  mutate(Parameter2 = case_when(Parameter=="alpha12" ~ "$\\alpha_{12}$",
                                Parameter=="alpha23" ~ "$\\alpha_{23}$",
                                Parameter=="beta0" ~ "$\\beta_{0}$",
                                Parameter=="beta1" ~ "$\\beta_{1}$",
                                Parameter=="beta2" ~ "$\\beta_{2}$",
                                Parameter=="d11" ~ "$d_{11}$",
                                Parameter=="d21" ~ "$d_{21}$",
                                Parameter=="d22" ~ "$d_{22}$",
                                Parameter=="lambda121"  ~ "$\\xi_{121}$",
                                Parameter=="lambda122"  ~ "$\\xi_{122}$",
                                Parameter=="lambda123"  ~ "$\\xi_{123}$",
                                Parameter=="lambda231"  ~ "$\\xi_{231}$",
                                Parameter=="lambda232" ~ "$\\xi_{232}$",
                                Parameter=="lambda233"  ~ "$\\xi_{233}$",
                                Parameter=="sigma" ~ "$\\sigma$",
                                Parameter=="gamma121" ~"$\\gamma_{121}$",
                                Parameter=="gamma122" ~ "$\\gamma_{122}$",
                                Parameter=="gamma231" ~ "$\\gamma_{231}$",
                                Parameter=="bs1(1)"    ~ "$\\kappa_{120}$",
                                Parameter=="bs1(2)"    ~ "$\\kappa_{230}$",
                                Parameter=="bs2(1)"    ~"$\\kappa_{121}$",
                                Parameter=="bs2(2)"    ~ "$\\kappa_{231}$",
                                Parameter=="bs3(1)"    ~ "$\\kappa_{122}$",
                                Parameter=="bs3(2)"   ~ "$\\kappa_{232}$",
                                Parameter=="bs4(1)"    ~ "$\\kappa_{123}$",
                                Parameter=="bs4(2)" ~ "$\\kappa_{233}$",
                                Parameter=="bs5(1)"   ~ "$\\kappa_{124}$",
                                Parameter=="bs5(2)"   ~ "$\\kappa_{234}$",
                                Parameter=="bs6(1)"    ~ "$\\kappa_{125}$",
                                Parameter=="bs6(2)" ~ "$\\kappa_{235}$"
                                
  )) %>%
  relocate(Parameter2)

results.latex.DGP1 <- results.latex2 %>% filter(DGP==1) %>% dplyr::select(-c(DGP,Parameter)) 

kable(results.latex.DGP1, "latex", booktabs=TRUE, caption = "Simulation results for data generating process 1", 
      label = "DGP1",
      escape = FALSE,
      align=c("l","l",rep("c",18)),
      col.names = c("Parameter","True",
                    "Mean","Bias","EmpSE","ModSE","MSE","Cover",
                    "Mean","Bias","EmpSE","ModSE","MSE","Cover",
                    "Mean","Bias","EmpSE","ModSE","MSE","Cover")) %>%
  kable_classic() %>%
  kable_styling(latex_options = c("scale_down"), position = "center") %>%
  add_header_above(c(" "=2,"NLMIXED"=6, "JMStateModel"=6,"MSM"=6)) %>%
  landscape() %>%
  pack_rows("Association parameters", 1, 2) %>%
  pack_rows("Other 1 textrightarrow 2 transition intensity parameters", 3, 7) %>%
  pack_rows("Other 2 textrightarrow 3 transition intensity parameters", 8, 11) %>%
  pack_rows("Linear mixed-effects model parameters", 12, 15) %>%
  pack_rows("Random effects parameters", 16, 18) %>%
  pack_rows("Bspline parameters (only for JMStateModel)", 19, 30) %>%
  footnote(general = "Monte Carlo standard errors shown in parantheses, when applicable. NLMIXED, the proposed joint model; JMStateModel, comparator 1; MSM, comparator (see section 4.4). EmpSE, empirical standard error; ModSE, model standard error; MSE, mean squared error.",
           footnote_as_chunk = TRUE)
