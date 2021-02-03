# mvpa-methods
This is a collection of R scripts developed by Daniel A. Stehr (Daniel.A.Stehr@dartmouth.edu) that accompany Stehr et al. (2021). The purpose was to evaluate the impact of four data processing steps on multivariate pattern classification in event-related designs. In a nutshell, we found large improvements in classification accuracy when the data - individual trial least-squares separate (or LSS) beta estimates - were averaged by condition and mean centered within runs. These improvements were observed for both simulated data and human participant fMRI data in an ROI that we expected to contain signal (left somato-motor region, classifying right hand button presses while participants performed a visual discrimination task). 

The human participant data in the form of LSS beta estimates for two ROIs (left somato-motor and right auditory cortex) are stored in five separate .RData files named following the convention "data-roi_clean-X.RData" where X stands for the type of motion-related nuisance regression performed while estimating LSS beta estimates. Each file contains the pattern data from both ROIs for 24 participants. Matrices containing each participant's motion-related nuisance regressors are included when applicable. 

The script "run_classification.R" can be used to classify the button presses using all of the methods proposed in the paper. The latter half of the script will also simulate beta patterns using a multi-level modelling framework and classify the result. The bulk of this work is carried out by two functions sourced from "mvpa_functions.R" Below is a summary of what each function does. 

simulate.betas - simulates beta patterns using a multi-level model with random effects parameters set by user

classify.betas - classifies (real or simulated) beta patterns using any combination of the three methods proposed by the paper (trial averaging, run-wise mean centering, cost tuning)
