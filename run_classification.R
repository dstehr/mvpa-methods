
#Customize path for your own computer
setwd("/data1/2018_ActionDecoding/analysis_class/share")

source('mvpa_funcs.R')
library(dplyr)

# Classify human participant fMRI dataset ---------------------------------
load('data-roi_clean-none.RData')
#Extract data for the first subject from both somato-motor and auditory ROIs
human.data <- rbind(data[[1]]$betas.sommot,data[[1]]$betas.auditory)
#Perform leave-one-run out cross-validated classification 
#   on the type of button pressed on each trial
#   modify the last three parameters as you wish (see comments above function in 'mvpa_funcs.R')
human.result <- classify.betas(betas=human.data,labToClassify = 'resp',
                               costVal = 2^0, meanCtrRuns = T, avgTrial = '1-avg')


# Classify simulated dataset ----------------------------------------------
set.seed(123) #For reproducibility, comment out if desired ...
simulated.data <- simulate.betas(nvox=200,ntrials=24,nrun=8,
                   etvs=.8,
                   mu.alpha=1,mu.beta=0,
                   tau.alpha=1, tau.beta=.05,
                   gamma.run=1.5)
#Concatenate patters from all simulated runs
simulated.data <- do.call('rbind',simulated.data$pattern)
#Perform leave-one-run out cross-validated classification 
#   on the type of button pressed on each trial
#   modify the last three parameters as you wish (see comments above function in 'mvpa_funcs.R')
simulated.result <- classify.betas(betas=simulated.data,labToClassify = 'cond',
                                   costVal = 2^0, meanCtrRuns = T, avgTrial = '1-avg')






