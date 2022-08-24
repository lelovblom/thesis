library(rstan)

# THIS MODEL COULD TAKE SEVERAL DAYS TO RUN
# THE RESULTING MODEL OBJECT IS VERY LARGE (~500 MB), AND HAS BEEN INCLUDED IN THE CHAPTER 5 FILES
# SO RUNNING THIS MODEL AND SAVING THE OUTPUT IS OPTIONAL

setwd("M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5")
inits_paper3 <- readRDS('inits_paper3.rds')
standata_paper3 <- readRDS('standata_paper3.rds')

stan_out_1 <- stan(
  file = "M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER5/stan_DGP1.stan",
  data = standata_paper3,
  init = function() inits_paper3,
  seed = 12345,
  chains = 4, iter = 700, warmup = 300,
  cores = 4, 
  # sample_file = sample_file_path # OPTION TO STORE THE SAMPLES AS THEY RUN, FOR MONITORING
  control = list(adapt_delta=0.85)
  )

saveRDS(stan_out_1,'stan_out_1.rds')