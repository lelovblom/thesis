# ----------------------------------------------------------------------------------------------------------
# THIS FILE INCLUDES 2 FUNCTIONS USED TO GENERATE SIMULATED DATA
# THEY ARE BASED ON APPROACHES OUTLINED BY FERRER ET AL., PAPGEORGIOU, AND COOK & LAWLESS

# Ferrer, L., Rondeau, V., Dignam, J., Pickles, T., Jacqmin-Gadda, H., & Proust-Lima, C. (2016). Joint
# modelling of longitudinal and multi-state processes: Application to clinical progressions in prostate
# cancer. Statistics in medicine, 35(22), 3933-3948.
# 
# Papageorgiou, G. (2022). Multi-state processes. Retrieved April 4, 2022, from https://drizopoulos.github.
# io/JMbayes2/articles/Multi State Processes.html
# 
# Cook, R. J., & Lawless, J. F. (2018). Multistate models for the analysis of life history data. CRC Press.
# SEE https://www.math.uwaterloo.ca/~rjcook/cook-lawless-multistate.html
# ----------------------------------------------------------------------------------------------------------

DGP_FUNCTION_3 <- function(N,
                           n,D11,D21,D22,true.betas,sigma.e,
                           alpha,pwc_knot1,pwc_knot2,lambda,
                           CA,RandomCensoring,prop.censored,
                           rho,X1mean, X1sd, X2min, X2max, X3p, gammas) {
  # [1] LONGITUDINAL DATA
  # vector of ids
  id <- rep(1:N, each = n)
  # minimum and maximum follow-up times  
  min.t <- 0
  max.t <- CA # ADMINISTRATIVE CENSORING TIME
  # sample time-points
  time <- replicate(N, c(0, sort(runif(n - 1, min = min.t, max = max.t))), simplify = FALSE)
  time <- do.call(c, time)
  #---------------------------------------------------------------------------------------------
  ############## sample 1 GAUSSIAN COVARIATE
  X1cont.s <- rnorm(N, mean = X1mean, sd = X1sd) # wide version
  X1cont <- rep(X1cont.s, each = n) # long
  #-------------------------------------------------------------------------------------------
  ############## Now sample the SKEWED COVARIATE THAT HAS NO EFFECT...THIS IS UNNECESSARY BUT I KEPT IT IN THE CODE....
  var1 <- runif(N,min=X2min,max=X2max)
  var2 <- runif(N,min=X2min,max=X2max)
  X2cont.s <- pmin(var1,var2) # wide version
  # hist(X2cont.s)
  X2cont <- rep(X2cont.s, each = n) # long
  #-------------------------------------------------------------------------------------------
  ############## Now sample the binary covariate (no effect in the longitudinal submodel though)
  X3cat.s <- rbinom(N, 1, X3p) # wide version
  X3cat <- rep(X3cat.s, each = n) # long
  # table(X3cat.s)
  #-------------------------------------------------------------------------------------------
  # initiate data frame to store results
  DF.long <- data.frame("id" = id, "time" = time, "X1cont" = X1cont, "X2cont" = X2cont, "X3cat" = X3cat)
  # design matrices for fixed and random effects  
  X <- model.matrix(~ 1 + time + X1cont, data = DF.long)
  Z <- model.matrix(~ 1 + time, data = DF.long)
  D <- matrix(c(D11,D21,D21,D22),2)
  b <- mvrnorm(N, mu = c(0,0), Sigma = D )
  eta.y <- as.vector(X %*% true.betas + rowSums(Z * b[id, ]))
  DF.long$y <- rnorm(N * n, eta.y, sigma.e)
  
  # [2] MULTISTATE TRANSITION DATA
  # design matrix transition intensities
  W <- cbind(X1cont[seq(1, by = n, N*n)],X3cat[seq(1, by = n, N*n)], 
             X1cont[seq(1, by = n, N*n)])
  ## linear predictor for transition: 1 -> 2
  eta.t1 <- as.vector(W[, c(1,2), drop = FALSE] %*% gammas[1:2])
  ## linear predictor for transition: 2 -> 3
  eta.t2 <- as.vector(W[, c(3), drop = FALSE] %*% gammas[3])
  
  # we simulate event times for 1->2 using inverse transform sampling
  invS12 <- function(t, u, i) { 
    h <- function(s) {
      XX <- cbind(1, s, X1cont.s[i])                                            
      ZZ <- cbind(1, s)
      f1 <- as.vector(XX %*% true.betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
      exp( lambda["lambda.121"]*(s<=pwc_knot1) + lambda["lambda.122"]*(pwc_knot1<s & s<=pwc_knot2) + lambda["lambda.123"]*(s>pwc_knot2) + eta.t1[i] + f1*alpha["alpha.12"] )
    }
    integrate(h, lower = 0, upper = t, subdivisions = 10000L)$value + log(u)
  }
  # likewise for event times for 2->3 using inverse transform sampling
  invS23 <- function (t, u, i) {
    h <- function (s) {
      XX <- cbind(1, s, X1cont.s[i])                                            
      ZZ <- cbind(1, s)
      f1 <- as.vector(XX %*% true.betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
      exp( lambda["lambda.231"]*(s<=pwc_knot1) + lambda["lambda.232"]*(pwc_knot1<s & s<=pwc_knot2) + lambda["lambda.233"]*(s>pwc_knot2) + eta.t2[i] + f1*alpha["alpha.23"] )
    }
    integrate(h, lower = 0, upper = t, subdivisions = 10000)$value + log(u)
  }
  # Probability for each transition
  u12 <- runif(N, 0, 1)
  u23 <- runif(N, 0, 1)
  # initiate vectors to save true event times
  trueT12 <- numeric(N)
  trueT23 <- numeric(N)
  # sample censoring times
  CAtimes <- rep(CA,N)
  CRtimes <- runif(N, 5, 20/prop.censored)
  # simulate time-to-event data
  for (i in 1:N) {
    Root12 <- NULL
    Root23 <- NULL
    
    Up <- 50
    tries <- 5
    # Transition 1->2
    Up <- 200
    Root12 <- try(uniroot(invS12, interval = c(1e-05, Up), u = u12[i], i = i)$root, TRUE)
    trueT12[i] <- if (!inherits(Root12, "try-error")) Root12 else 500
    
    # Transition 2->3
    if(as.numeric(trueT12[i]) < CAtimes[i]) {
      Up <- Up + 200
      Root23 <- try(uniroot(invS23, interval = c(as.numeric(trueT12[i]), Up), u = u23[i], i = i)$root, TRUE)
    } else {Root23 <- 500}
    trueT23[i] <- if (!inherits(Root23, "try-error")) Root23 else 500
  }
  # initiate multi-state dataset in wide format
  DF.trans <- data.frame('id' = 1:N, 'trueT12' = trueT12, 'trueT23' = trueT23, #'trueT12' = trueT12, 
                         'CAtimes' = CAtimes, 'CRtimes'= CRtimes,
                         'X1cont' = X1cont.s, 'X2cont' = X2cont.s, 'X3cat' = X3cat.s)
  
  # [3] MULTISTATE VISIT PROCESS
  generatedata.f <- function(N, CA, rho) { 
    outdata <- lapply(1:N, function(id, rho, CA) { #  nstates, init.pi, qmat, 
      ai <- 0
      Zai.present <- 1
      
      r <- 1;  air <- ai[1]
      while (air < CA) {
        r <- r + 1
        Delta.air <- rexp(1, rate=rho)
        air <- Delta.air + ai[r-1]
        
        if (air < CA) {
          ai <- c(ai, air)
          Zai.present <- c(Zai.present, 1)
        }
        else { break }
      }    
      
      ai <- c(ai, CA)
      Zai.present <- c(Zai.present, 999)
      
      return( data.frame(id=rep(id, length(ai)), times=ai, states.measured=Zai.present) ) 
    }, rho=rho, CA=CA)  
    outdata <- do.call("rbind", outdata)
    outdata <- outdata[order(outdata$id, outdata$times),]
    return(outdata)   
  }
  DF.visits <- generatedata.f(N, CA, rho)
  
  # [4] CLEAN AND MERGE TRANSITION DATA AND VISIT DATA TO PRODUCE THE MULTISTATE DATA
  data2 <- inner_join(DF.trans,DF.visits,by="id")
  # remove observation times occurring after censoring and/or after transition into state 3
  data3 <- subset(data2, times<=CAtimes & times<trueT23) #, select=-X) 
  # define the state observed at each time, rename time variable, and output the data
  data4 <- data3 %>% mutate(state = case_when( states.measured==1 & times < trueT12 ~ 1,
                                               states.measured==1 & times >= trueT12 & times < trueT23 ~ 2,
                                               states.measured==999 ~ 999)) %>%
    rename(year=times) %>%
    dplyr::select(id,year,state,X1cont,X2cont,X3cat)                                             
  # create a row for entry into state 3
  timesT23 <- subset(DF.trans,trueT23!=500 & trueT23<=CAtimes, select=c(id,trueT23,X1cont,X2cont,X3cat))   
  timesT23$state <- 3
  timesT23 <- timesT23 %>% rename(year=trueT23) %>%
    dplyr::select(id,year,state,X1cont,X2cont,X3cat)                                                     
  # now bind the rows, and sort
  data5 <- bind_rows(data4,timesT23)
  # data5 <- data5 %>% arrange(id,year)
  data5 <- data5[order(data5$id, data5$year),]
  
  # NOW INCORPORATE RANDOM CENSORING, IF PRESENT
  if (RandomCensoring==TRUE) {
    timesRC <- subset(DF.trans, CRtimes<trueT23 & CRtimes<CAtimes)
    timesRC$state <- 999
    timesRC <- timesRC %>% rename(year=CRtimes) %>%
      dplyr::select(id,year,state,X1cont,X2cont,X3cat)
    timesRC2 <- subset(DF.trans, CRtimes<trueT23 & CRtimes<CAtimes)
    timesRC2 <- timesRC2 %>% dplyr::select(id,CRtimes)
    data5.2 <- left_join(data5, timesRC2, by = "id")
    data5.3 <- data5.2 %>% filter(is.na(CRtimes) | (!is.na(CRtimes) & year<CRtimes)) %>% dplyr::select(-CRtimes)
    data6 <- bind_rows(data5.3,timesRC)
    data6 <- data6[order(data6$id, data6$year),]
  } else{
    data6 <- data5
  }
  
  # [5] CLEAN MULTISTATE DATA
  # REMOVE PEOPLE WITH A SINGLE MULTISTATE OBSERVATION
  data8 <- data6 %>% group_by(id) %>% filter(n()>1)
  
  # [6] CLEAN LONGITUDINAL DATA
  # REMOVE LONGITUDINAL DATA THAT OCCURS AFTER THE FINAL MULTISTATE OBSERVATION
  Tstop_Multistate <- tapply(data8$year, data8$id, max)
  Tstop_Multistate <- Tstop_Multistate[id]
  DF.long2 <- DF.long[DF.long$time <= Tstop_Multistate, ]
  
  # [7] STACK THE CLEAN MULTISTATE AND LONGITUDINAL DATA
  # ENSURE THE SAME PEOPLE ARE INCLUDED IN BOTH DATASETS
  data8.ids <- data8 %>% group_by(id) %>% filter(row_number()==1) %>% dplyr::select(id)
  DF.long2.ids <- DF.long2 %>% group_by(id) %>% filter(row_number()==1) %>% dplyr::select(id)
  inboth <- inner_join(data8.ids, DF.long2.ids, by="id")
  data9 <- inner_join(data8, inboth, by="id")
  DF.long3 <- inner_join(DF.long2, inboth, by="id")
  
  # NOW ADD SOME VARIABLES FOR SAS AND RBIND
  data9$outcome_number <- 1
  DF.long3 <- DF.long3 %>% rename(year=time)
  DF.long3$outcome_number <- 2
  # data_jm <- rbind(data9,DF.long3)
  data_jm <- bind_rows(data9,DF.long3)
  data_jm <- data_jm %>% arrange(id,outcome_number,year)
  
  # [8] AND RETURN THAT DATA SET
  return(data_jm)
}


DGP_FUNCTION_4 <- 
  function(N,n, D11, D21, D22, true.betas, sigma.e, 
           alpha, pwc_knot1, pwc_knot2, lambda,
           CA, RandomCensoring, prop.censored,rho,
           X1mean, X1sd, X2min, X2max, X3p, gammas){
    # [1] LONGITUDINAL DATA - NOT MUCH TO CHANGE HERE...
    # vector of ids
    id <- rep(1:N, each = n)
    # # minimum and maximum follow-up times  
    # min.t <- 0
    # max.t <- CA # ADMINISTRATIVE CENSORING TIME
    # # sample time-points
    # time <- replicate(N, c(0, sort(runif(n - 1, min = min.t, max = max.t))), simplify = FALSE)
    # time <- do.call(c, time)
    #---------------------------------------------------------------------------------------------
    # CODE CHANGE COMPARED TO ABOVE
    #---------------------------------------------------------------------------------------------
    # study start times
    startvar1 <- runif(N,min=0,max=10)
    startvar2 <- runif(N,min=0,max=10)
    Study.start.time.s <- pmin(startvar1,startvar2) # wide version
    # hist(Study.start.time.s)
    Study.start.time <- rep(Study.start.time.s, each = n) # long
    # minimum and maximum follow-up times  
    min.t <- Study.start.time.s
    max.t <- CA # ADMINISTRATIVE CENSORING TIME
    # sample time-points
    time <- list()
    for (i in 1:N) {
      time[[i]] <- c(min.t[i], sort(runif(n - 1, min = min.t[i], max = min.t[i] + max.t))) # this is the study time, relative to the true time 0
    }
    time <- do.call(c, time)
    #---------------------------------------------------------------------------------------------
    # RETURN TO PREVIOUS CODE
    #---------------------------------------------------------------------------------------------
    ############## sample 1 GAUSSIAN COVARIATE - THIS IS TO REPRESENT THE A1C AT RANDOMIZATION
    X1cont.s <- rnorm(N, mean = X1mean, sd = X1sd) # wide version
    X1cont <- rep(X1cont.s, each = n) # long
    # ############## NOW SAMPLE THE TRUE TIMEZERO COVARIATE: assume it's the same 
    X1cont0 <- X1cont 
    X1cont0.s <- X1cont.s 
    #-------------------------------------------------------------------------------------------
    ############## Now sample the SKEWED COVARIATE THAT HAS NO EFFECT
    var1 <- runif(N,min=X2min,max=X2max)
    var2 <- runif(N,min=X2min,max=X2max)
    X2cont.s <- pmin(var1,var2) # wide version
    # hist(X2cont.s)
    X2cont <- rep(X2cont.s, each = n) # long
    #-------------------------------------------------------------------------------------------
    ############## Now sample the binary covariate (no effect in the longitudinal submodel though)
    X3cat.s <- rbinom(N, 1, X3p) # wide version
    X3cat <- rep(X3cat.s, each = n) # long
    # table(X3cat.s)
    #-------------------------------------------------------------------------------------------
    # initiate data frame to store results
    DF.long <- data.frame("id" = id, "time" = time, "X1cont" = X1cont, "X2cont" = X2cont, "X3cat" = X3cat, "X1cont0" = X1cont0)
    # design matrices for fixed and random effects  
    X <- model.matrix(~ 1 + time + X1cont, data = DF.long)
    Z <- model.matrix(~ 1 + time, data = DF.long)
    D <- matrix(c(D11,D21,D21,D22),2)
    b <- mvrnorm(N, mu = c(0,0), Sigma = D )
    eta.y <- as.vector(X %*% true.betas + rowSums(Z * b[id, ]))
    DF.long$y <- rnorm(N * n, eta.y, sigma.e)
    # [2] MULTISTATE TRANSITION DATA - THE HARD PART IS TO CREATE A TIME-VARYING TREATMENT AND A1C VARIABLE
    # design matrix transition intensities from study start and onwards
    W <- cbind(X1cont[seq(1, by = n, N*n)],X3cat[seq(1, by = n, N*n)], 
               X1cont[seq(1, by = n, N*n)])
    # design matrix for the transition intensities from true zero to study start:
    W0 <- cbind(X1cont0[seq(1, by = n, N*n)],X1cont0[seq(1, by = n, N*n)])
    # eta.t1 <- function(s,i) {
    #   if (s < Study.start.time.s[i]) {
    #     return(as.vector(W0[i, c(1), drop = FALSE] %*% gammas[1]))
    #   } else {
    #     return(as.vector(W[i, c(1,2), drop = FALSE] %*% gammas[1:2]))
    #   }
    # }
    # we simulate event times for 1->2 using inverse transform sampling
    invS12 <- function(t, u, i) { 
      h <- function(s) {
        eta.t1 <- (s < Study.start.time.s[i])*(as.vector(W0[i, c(1), drop = FALSE] %*% gammas[1])) + 
          (s >= Study.start.time.s[i])*(as.vector(W[i, c(1,2), drop = FALSE] %*% gammas[1:2]))
        XX <- (s < Study.start.time.s[i])*(cbind(1, s, X1cont0.s[i])) + 
          (s >= Study.start.time.s[i])*cbind(1, s, X1cont.s[i]) 
        ZZ <- cbind(1, s)
        f1 <- as.vector(XX %*% true.betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ])) # the form of XX, ZZ, betas, and the random effects are always the same
        # need to fix the knots too now! We will assume that the first knot extends to the true time zero
        pwc_knot1_adj <- pwc_knot1 + Study.start.time.s[i]
        pwc_knot2_adj <- pwc_knot2 + Study.start.time.s[i]
        exp( lambda["lambda.121"]*(s<=pwc_knot1_adj) + lambda["lambda.122"]*(pwc_knot1_adj<s & s<=pwc_knot2_adj) + lambda["lambda.123"]*(s>pwc_knot2_adj) + 
               eta.t1 + f1*alpha["alpha.12"] )
      }
      integrate(h, lower = 0, upper = t, subdivisions = 10000L)$value + log(u)
    }
    # likewise for event times for 2->3 using inverse transform sampling
    invS23 <- function (t, u, i) {
      h <- function (s) {
        eta.t2 <- (s < Study.start.time.s[i])*(as.vector(W0[i, c(2), drop = FALSE] %*% gammas[3])) + 
          (s >= Study.start.time.s[i])*(as.vector(W[i, c(3), drop = FALSE] %*% gammas[3]))
        XX <- (s < Study.start.time.s[i])*(cbind(1, s, X1cont0.s[i])) + 
          (s >= Study.start.time.s[i])*cbind(1, s, X1cont.s[i]) 
        ZZ <- cbind(1, s)
        f1 <- as.vector(XX %*% true.betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
        pwc_knot1_adj <- pwc_knot1 + Study.start.time.s[i]
        pwc_knot2_adj <- pwc_knot2 + Study.start.time.s[i]
        exp( lambda["lambda.231"]*(s<=pwc_knot1_adj) + lambda["lambda.232"]*(pwc_knot1_adj<s & s<=pwc_knot2_adj) + lambda["lambda.233"]*(s>pwc_knot2_adj) + 
               eta.t2 + f1*alpha["alpha.23"] )
      }
      integrate(h, lower = 0, upper = t, subdivisions = 10000)$value + log(u)
    }
    # Probability for each transition
    u12 <- runif(N, 0, 1)
    u23 <- runif(N, 0, 1)
    # initiate vectors to save true event times
    trueT12 <- numeric(N)
    trueT23 <- numeric(N)
    #---------------------------------------------------------------------------------------------
    # CODE CHANGE
    #---------------------------------------------------------------------------------------------
    # sample censoring times
    # CAtimes <- rep(CA+10,N)
    CAtimes <- Study.start.time.s + CA # instead of 25, it's 25 after the start time, the max being 35 
    #---------------------------------------------------------------------------------------------
    # RETURN
    #---------------------------------------------------------------------------------------------
    CRtimes <- runif(N, 5, 20/prop.censored)
    # simulate time-to-event data
    for (i in 1:N) {
      Root12 <- NULL
      Root23 <- NULL
      
      Up <- 50
      tries <- 5
      # Transition 1->2
      Up <- 200
      Root12 <- try(uniroot(invS12, interval = c(1e-05, Up), u = u12[i], i = i)$root, TRUE)
      trueT12[i] <- if (!inherits(Root12, "try-error")) Root12 else 500
      
      # Transition 2->3
      if(as.numeric(trueT12[i]) < CAtimes[i]) {
        Up <- Up + 200
        Root23 <- try(uniroot(invS23, interval = c(as.numeric(trueT12[i]), Up), u = u23[i], i = i)$root, TRUE)
      } else {Root23 <- 500}
      trueT23[i] <- if (!inherits(Root23, "try-error")) Root23 else 500
    }
    #---------------------------------------------------------------------------------------------
    # initiate multi-state dataset in wide format
    DF.trans <- data.frame('id' = 1:N, 'trueT12' = trueT12, 'trueT23' = trueT23, #'trueT12' = trueT12, 
                           'CAtimes' = CAtimes, 'CRtimes'= CRtimes,
                           'X1cont' = X1cont.s, 'X2cont' = X2cont.s, 'X3cat' = X3cat.s,
                           'timezero'=Study.start.time.s,
                           'X1cont0'=X1cont0.s)
    #---------------------------------------------------------------------------------------------
    # [3] MULTISTATE VISIT PROCESS
    generatedata.f <- function(N, CA, rho) { 
      outdata <- lapply(1:N, function(id, rho, CA) { 
        ai <- Study.start.time.s[id]
        Zai.present <- 1
        
        r <- 1;  air <- ai[1]
        #---------------------------------------------------------------------------------------------
        while (air < Study.start.time.s[id] + CA) {
          #---------------------------------------------------------------------------------------------
          r <- r + 1
          Delta.air <- rexp(1, rate=rho)
          air <- Delta.air + ai[r-1]
          #---------------------------------------------------------------------------------------------
          if (air < Study.start.time.s[id] + CA) {
            #---------------------------------------------------------------------------------------------
            ai <- c(ai, air)
            Zai.present <- c(Zai.present, 1)
          }
          else { break }
        }    
        #---------------------------------------------------------------------------------------------
        ai <- c(ai, Study.start.time.s[id] + CA)
        #---------------------------------------------------------------------------------------------
        Zai.present <- c(Zai.present, 999)
        
        return( data.frame(id=rep(id, length(ai)), times=ai, states.measured=Zai.present) ) 
      }, rho=rho, CA=CA)  
      outdata <- do.call("rbind", outdata)
      outdata <- outdata[order(outdata$id, outdata$times),]
      return(outdata)   
    }
    DF.visits <- generatedata.f(N, CA, rho)
    
    # [4] CLEAN AND MERGE TRANSITION DATA AND VISIT DATA TO PRODUCE THE MULTISTATE DATA
    data2 <- inner_join(DF.trans,DF.visits,by="id")
    # remove observation times occurring after censoring and/or after transition into state 3
    data3 <- subset(data2, times<=CAtimes & times<trueT23) #, select=-X) 
    # define the state observed at each time, rename time variable, and output the data
    data4 <- data3 %>% mutate(state = case_when( states.measured==1 & times < trueT12 ~ 1,
                                                 states.measured==1 & times >= trueT12 & times < trueT23 ~ 2,
                                                 states.measured==999 ~ 999)) %>%
      #rename(year=times) %>%
      dplyr::select(id,state,X1cont,X2cont,X3cat,times,timezero,X1cont0) 
    # create a row for entry into state 3
    timesT23 <- subset(DF.trans,trueT23!=500 & trueT23<=CAtimes, select=c(id,trueT23,X1cont,X2cont,X3cat,timezero,X1cont0))   
    timesT23$state <- 3
    timesT23 <- timesT23 %>% rename(times=trueT23) %>%
      dplyr::select(id,state,X1cont,X2cont,X3cat,times,timezero,X1cont0)
    # now bind the rows, and sort
    data5 <- bind_rows(data4,timesT23)
    data5 <- data5[order(data5$id, data5$times),]
    
    # NOW INCORPORATE RANDOM CENSORING, IF PRESENT
    if (RandomCensoring==TRUE) {
      timesRC <- subset(DF.trans, CRtimes<trueT23 & CRtimes<CAtimes)
      timesRC$state <- 999
      timesRC <- timesRC %>% rename(year=CRtimes) %>%
        dplyr::select(id,year,state,X1cont,X2cont,X3cat,X1cont0)
      timesRC2 <- subset(DF.trans, CRtimes<trueT23 & CRtimes<CAtimes)
      timesRC2 <- timesRC2 %>% dplyr::select(id,CRtimes)
      data5.2 <- left_join(data5, timesRC2, by = "id")
      data5.3 <- data5.2 %>% filter(is.na(CRtimes) | (!is.na(CRtimes) & year<CRtimes)) %>% dplyr::select(-CRtimes)
      data6 <- bind_rows(data5.3,timesRC)
      data6 <- data6[order(data6$id, data6$year),]
    } else{
      data6 <- data5
    }
    
    # [5] CLEAN MULTISTATE DATA
    # REMOVE PEOPLE WITH A SINGLE MULTISTATE OBSERVATION
    data8 <- data6 %>% group_by(id) %>% filter(n()>1)   # E.G. PEOPLE WHO TRANSITION INTO STATE 3 PRIOR TO START TIME
    
    # [6] CLEAN LONGITUDINAL DATA
    # WE HAVE HAD TO EXLUCDE SOME PEOPLE FROM THE MULTISTATE PROCESS. LET'S DELETE THEM FROM THE LONGITUDINAL DATA
    data8.ids <- data8 %>% group_by(id) %>% filter(row_number()==1) %>% dplyr::select(id) 
    DF.long2 <- inner_join(DF.long,data8.ids, by = "id")
    # REMOVE LONGITUDINAL DATA THAT OCCURS AFTER THE FINAL MULTISTATE OBSERVATION
    Tstop_Multistate <- tapply(data8$times, data8$id, max)
    Tstop_Multistate2 <- data.frame("id"=data8.ids,"Tstop_Multistate"=Tstop_Multistate)
    DF.long3 <- inner_join(DF.long2,Tstop_Multistate2, by = "id")
    DF.long4 <- DF.long3 %>% filter(time<=Tstop_Multistate)
    
    # [7] STACK THE CLEAN MULTISTATE AND LONGITUDINAL DATA
    # ENSURE THE SAME PEOPLE AND COVARIATES (INCLUDING TIMEZERO) ARE INCLUDED IN BOTH DATASETS 
    data8.ids <- data8 %>% group_by(id) %>% filter(row_number()==1) %>% dplyr::select(id)
    DF.long4.ids <- DF.long4 %>% group_by(id) %>% filter(row_number()==1) %>% dplyr::select(id)
    inboth <- inner_join(data8.ids, DF.long4.ids, by="id")
    data9 <- inner_join(data8, inboth, by="id")
    DF.long5 <- inner_join(DF.long4, inboth, by="id")
    # ADD TIME ZERO TO DF.LONG5
    data8.timezero.ids <- data8 %>% group_by(id) %>% filter(row_number()==1) %>% dplyr::select(id,timezero)
    DF.long6 <- inner_join(DF.long5, data8.timezero.ids, by="id")
    
    # NOW ADD SOME VARIABLES FOR SAS AND RBIND
    data9$outcome_number <- 1
    DF.long6 <- DF.long6 %>% rename(times=time)
    DF.long6$outcome_number <- 2
    # data_jm <- rbind(data9,DF.long6)
    data_jm <- bind_rows(data9,DF.long6)
    data_jm <- data_jm %>% arrange(id,outcome_number,times)
    
    data_jm$year <- data_jm$times - data_jm$timezero
    
    # [8] AND RETURN THAT DATA SET
    return(data_jm)
  }
