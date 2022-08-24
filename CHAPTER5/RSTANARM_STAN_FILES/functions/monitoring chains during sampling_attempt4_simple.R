# https://github.com/stan-dev/rstan/issues/591
library(tidyverse)
library(rstan)

samples_from_chain = function(filename) {
  stanfit = rstan::read_stan_csv(filename)
  samples = data.frame(stanfit@sim$samples)
  samples %>% 
    filter(lp__ != 0) %>%  # Remove non-sampled iterations
    mutate(
      iter = 1:nrow(.),
      chain = gsub('.csv', '', last(strsplit(filename, '_')[[1]]))  # Get chain number [filename_1.csv]
    )
}

setwd('/Users/lovblom/Documents/PHD/DATASETS/Erik Data Processing/R Paper 3/RESTART 05JUN22')
setwd('M:/lovblom/PHD/DATASETS/Erik Data Processing/R Paper 3/RESTART 05JUN22')



chains = rbind(
  samples_from_chain('PAPER2_attemp4_sim_chain_data_1.csv'),
  samples_from_chain('PAPER2_attemp4_sim_chain_data_2.csv'),
  samples_from_chain('PAPER2_attemp4_sim_chain_data_3.csv')
)

chains %>%
  filter(iter>20) %>%
  # select(chain, iter, b_Intercept, sigma) %>%
  select(chain, iter, yAlpha1.1,yGamma1.1,yBeta1.1,yBeta1.2) %>%
  gather(key, value, -chain, -iter) %>%
  ggplot(aes(x=iter, y=value, color=factor(chain))) + 
  facet_wrap(~key, scales = 'free_y') +
  geom_line() 

chains %>%
  filter(iter>20) %>%
  # select(chain, iter, b_Intercept, sigma) %>%
  select(chain, iter, yAux1.1,a12_beta.1,a23_beta.1,
         bCov1.1, bCov1.2, bCov1.3) %>%
  gather(key, value, -chain, -iter) %>%
  ggplot(aes(x=iter, y=value, color=factor(chain))) + 
  facet_wrap(~key, scales = 'free_y') +
  geom_line() 

chains %>%
  filter(iter>20) %>%
  # select(chain, iter, b_Intercept, sigma) %>%
  select(chain, iter, e12_beta.1,e12_beta.2,e23_beta.1,
         e12_aux.1,e12_aux.2,e12_aux.3,
         e23_aux.1,e23_aux.2,e23_aux.3) %>%
  gather(key, value, -chain, -iter) %>%
  ggplot(aes(x=iter, y=value, color=factor(chain))) + 
  facet_wrap(~key, scales = 'free_y') +
  geom_line() 
