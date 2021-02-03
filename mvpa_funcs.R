# Created by Daniel A. Stehr
# University of California, Irvine
# dstehr@uci.edu
# 3 Feb 2021


# Simulate beta series using a multi-level modelling approach -------------
simulate.betas <- function(nvox=50,ntrials=24,nrun=8,
                           etvs=1,
                           mu.alpha=1,mu.beta=0,
                           tau.alpha=1, tau.beta=.1,
                           gamma.run=1.5){
# DESCRIPTION
#   Simulate multivariate pattern data using the multi-level mixed model
#   used in Stehr et al. (2021) 
#     
# USAGE
#   simulate.betas(nvox = 50, ntrials = 24, nrun = 8,
#                  etvs = 1, mu.alpha = 1, mu.beta = 0,
#                  tau.alpha = 1, tau.beta = .1,
#                  gamma.run = 1.5)
# 
# ARGUMENTS
#   nvox         The number of voxels to simulate
#   ntrials      The number of trials used in the simulation
#   nrun         The number of runs to simulate
#   etvs         The standard deviation of the Normally distributed
#                noise to add to each trial
#   mu.alpha     Mean coefficient for the intercept parameter sampled for each voxel
#   mu.beta      Mean coefficient for the slope parameter sampled for each voxel
#   tau.alpha    Standard deviation of intercepts sampled for each voxel
#   tau.beta     Standard deviation of slopes sampled for each voxel
#   gamma.run    The standard deviation of Normally distributed shift applied to each run
#   
# OUTPUT
#   A list containing the pattern data and random effects used in the simulation
#   pattern      Data frames containing the pattern data for each run of the simulation
#   tau          A two-column matrix containing the intercept and slope parameter
#                sampled for each voxel during creation of the pattern data
  
  require(MASS)
  
  #Create conditioning variable X
  xp <- rep(c(-.5,.5),each=ntrials/2)
  xp.allvox <- rep(xp,nvox)
  
  #Need this to keep track of voxels and trials
  vox.id <- rep(1:nvox,each=ntrials)
  trials.id <- rep(1:ntrials,nvox)
  
  #Preallocate list to store output
  data <- list(pattern=list(),
               tau=list())
  
  #Step 1: Generate mean shift for each run
  if(!is.null(gamma.run)){
    dc.run <- rnorm(nrun,mean=0,sd=gamma.run)
  }else{
    dc.run <- rep(0,nrun)
  }
  
  #Step 2: Generate the voxel-specific intercepts and slopes
  tau.cov <- matrix(c(tau.alpha^2,0,0,tau.beta^2),nrow=2)
  rfx <- mvrnorm(nvox, c(mu.alpha,mu.beta),tau.cov)
  rfx.alpha <- rfx[,1]    #mean activity for all conds
  rfx.beta <- rfx[,2]    #slope (or effect) of conditions
  
  for(run in 1:nrun){
    betahats <- vector()       #Activation estimates will get filled in here
    
    #Step 3: Use the voxel-specific intercepts and slopes to generate trial-specific estimates
    for(v in 1:nvox){
      betahats.loop <- (rfx.alpha[v]+dc.run[run]) + rfx.beta[v]*xp + rnorm(length(xp),sd=etvs)
      betahats <- c(betahats,betahats.loop)
    }
    
    #Organize the simulated patterns into data frame for output
    pattern.run <- data.frame(trialNum=trials.id,
                              scan=run,
                              vox=vox.id,
                              betaHat=betahats,
                              cond=xp.allvox,
                              roi='simulated-region')
    pattern.run$cond <- factor(pattern.run$cond,levels=c(-0.5,0.5),
                               labels=c('A','B'))
    
    data$pattern[[run]] <- pattern.run
    data$tau[[run]] <- rfx
    
  }
  return(data)
}


# Classify simulated or real beta series ----------------------------------
classify.betas <- function(betas=NULL,
                           labToClassify, 
                           costVal = 2^0,
                           meanCtrRuns=T,
                           avgTrials=NULL){
# DESCRIPTION
#   Classifies multivariate pattern data using SVM with leave-one-run out cross-validation
#   
# USAGE
#   classify.betas(betas, labToClassify, costVal = 2^0, meanCtrRuns = T, avgTrials = '1-avg')
#   
# ARGUMENTS
#   betas              A data frame in long format containing trial-by-trial activation estimates.
#                      It should have columns named "trialNum", "scan", "vox", "betaHat", and "roi"
#                      indicating the trial number, scan (or run) number, beta estimate, and ROI for 
#                      each observation. An additional column indicating the trial type is required 
#                      but can take any name.
#   labToClassify      The name of the column in "betas" containing the stimulus labels to classify.  
#   costVal            The value of the cost hyperparameter for the support vector classifier.
#                      If costVal is a numeric vector of values, this function will perform 
#                      cost tuning in an inner nested cross-validated fold. 
#   meanCtrRuns        Logical value specifying whether the data should be mean centered within runs
#                      by subtracting each voxel's run-level mean for both conditions from the beta 
#                      estimates of each trial
#   avgTrials          A string indicating whether (and how) trial-specific esimates should be averaged.
#                      If NULL, classification is performed on separate activation estimates for each trial.
#                      Use "1-avg" to average all trials to a single estimate per condition and run prior to
#                      classification. Use "2-avg" to randomly sample trials of each type within runs to produce
#                      2 average estimates per condition and run. 
#   
# OUTPUT
#   A list contining the results of the support vector classification
#     trueClass          The true class labels for each held-out test set
#     predClass          The predicted class lables for the held-out test set
#     fittedValues       The raw fitted (decision) values for the test set
#     svmWeight          A list of the SVM weights for each voxel in each CV fold
#     foldAccuracy       The classification accuracy within each CV fold
#     AUC                Area under the ROC curve computed across all CV folds
#     cost               The cost value that was either supplied by the user
#                        or obtained by tuning the parameter in a nested fold
#     roiSize            The size (in number of voxels) for the ROI
#     meanAccuracy       The mean cross-validated classifcation accuracy
                                            
  require(dplyr)
  require(e1071)
  require(tidyr)
  require(tools)
  require(ROCR)
  
  # set.seed(123)
  baseDir <- "/data1/2018_ActionDecoding/analysis_class"
  
  nRun <- length(unique(betas$scan))
  nTrialsPerRun <- length(unique(betas$trialNum))
  totalTrials <- nTrialsPerRun * nRun
  scanList <- as.numeric(unique(betas$scan))
  roiList <- unique(betas$roi)
  nRoi <- length(roiList)
  class1 <- levels(betas[,grep(labToClassify,names(betas))])[1]
  class2 <- levels(betas[,grep(labToClassify,names(betas))])[2]
  
  classif.svm <- list()
  
    print(paste0('Classifying ... ',labToClassify))
    
    for(roi in 1:nRoi){
      print(paste0('Working on roi ',roi,': ... ' ,roiList[roi]))
      
      betas.roi <- betas[betas$roi==roiList[roi],]
      classIdx <- grep(labToClassify,names(betas.roi))
      
      #Remove rows with missing values
      betas.roi <- na.omit(betas.roi)
      
      # Perform trial averaging, if asked
      if(avgTrials=='1-avg'){
        #Average ALL data by class by run
        betas.roi <- betas.roi %>%
          group_by(scan,vox,!!sym(labToClassify)) %>%
          summarize(betaHat=mean(betaHat),.groups='drop') %>%
          ungroup()
        
      }else if(avgTrials=='2-avg'){
        #Average HALF of the data by class
        betas.roi.sample.idx <- betas.roi %>%
          group_by(trialNum,scan,!!sym(labToClassify)) %>%
          summarize(count=n(),.groups='drop') %>%
          arrange(scan) %>%
          dplyr::select(-count) %>%
          ungroup() %>%
          group_by(scan,!!sym(labToClassify)) %>%
          sample_frac(size=0.5) %>%
          mutate(rand.samp=1)
        
        betas.roi <- left_join(betas.roi,betas.roi.sample.idx,by=c('trialNum','scan',labToClassify))
        betas.roi[is.na(betas.roi$rand.samp),'rand.samp'] <- 0
        
        betas.roi <- betas.roi %>%
          group_by(scan,vox,rand.samp,!!sym(labToClassify)) %>%
          summarize(betaHat=mean(betaHat),.groups='drop') %>%
          ungroup()
        
      }
      
      #Reshape data to make it (trial X vox) format
      betas.roi <- spread(betas.roi, vox, betaHat)
      betas.roi <- arrange(betas.roi,scan)
      
      firstVoxCol <- which(names(betas.roi)=='1')
      varsToExclude <- firstVoxCol - 1
      
      #Perform row-wise mean centering (if asked)
      if(meanCtrRuns==T){
        betas.roi.mcr <- data.frame()
        mcr.runs <- unique(betas.roi$scan)
        for(i in mcr.runs){
          mcr.rundata <- filter(betas.roi,scan==i)
          
          #Center based on BOTH classes
          mcr.rundata[,firstVoxCol:ncol(mcr.rundata)] <- mcr.rundata[,firstVoxCol:ncol(mcr.rundata)] - matrix(1,ncol=1,nrow=nrow(mcr.rundata))%*%apply(mcr.rundata[,firstVoxCol:ncol(mcr.rundata)],2,mean)
          
          betas.roi.mcr <- rbind(betas.roi.mcr,mcr.rundata)
        }
        betas.roi <- betas.roi.mcr
      }
      
      #Add class labels
      classIdx <- grep(labToClassify,names(betas.roi))
      
      classif.svm[[roiList[roi]]] <- list(trueClass = vector('list',length=nRun),
                                          predClass = vector('list',length=nRun),
                                          fittedValue = vector('list',length=nRun),
                                          svmWeight = vector('list',length=nRun),
                                          foldAccuracy = vector(),
                                          AUC = vector(),
                                          cost = vector(),
                                          roiSize = vector(),
                                          meanAccuracy = vector())
      
      roiSize <- ncol(betas.roi)-varsToExclude
      classif.svm[[roiList[roi]]]$roiSize <-  roiSize
      print(paste0('     ROI size: ',roiSize))
      
      #Leave one run out cross validation
      for(run in scanList){
        train <- filter(betas.roi,scan!=run)
        trainLabels <- dplyr::select(train,classIdx)
        names(trainLabels) <- 'trainLabels'
        train <- cbind(train[,-(1:varsToExclude)],trainLabels)
        train$trainLabels <- as.factor(train$trainLabels)
        
        test <- filter(betas.roi,scan==run)
        testLabels <- dplyr::select(test,classIdx)
        names(testLabels) <- 'testLabels'
        test <- cbind(test[,-(1:varsToExclude)],testLabels)
        test$testLabels <- as.factor(test$testLabels)
        
        #record true stim labels
        classif.svm[[roiList[roi]]]$trueClass[[run]] <- test$testLabels
        
        #TUNE the cost parameter, if asked
        if(length(costVal)==1){
          cost.best <- costVal[1]
        }else{
          
          splits.L <- scanList[!scanList==run]
          ncostval <- length(costVal)
          innerfold.counter <- 0
          cost.profile <- matrix(NA,nrow=length(splits.L),ncol=ncostval)
          
          #start inner CV fold
          for(L in splits.L){
            
            innerfold.counter <- innerfold.counter + 1
            
            #inner training split
            train.L <- betas.roi %>% 
              filter(scan!=run) %>% 
              filter(scan!=L)
            train.L.lab <- dplyr::select(train.L,classIdx)
            names(train.L.lab) <- 'classLabels'
            train.L <- cbind(train.L[,-(1:varsToExclude)],train.L.lab)
            
            #inner testing split
            test.L <- betas.roi %>% 
              filter(scan!=run) %>% 
              filter(scan==L)
            test.L.lab <- dplyr::select(test.L,classIdx)
            names(test.L.lab) <- 'classLabels'
            test.L <- cbind(test.L[,-(1:varsToExclude)],test.L.lab)
            
            #Search over all cost values in range
            for(c.par in seq_along(costVal)){
              svmfit.loop <- svm(trainLabels~.,data=train,kernel='linear',cost=costVal[c.par])
              svmpred.loop <- predict(svmfit.loop, test.L[,-(grep('classLabels',names(test.L)))])
              tree <- table(predict=svmpred.loop,truth=test.L$classLabels)
              ca.loop <- sum(diag(tree))/sum(tree)
              cost.profile[innerfold.counter,c.par] <- ca.loop
            }
          }    #END inner CV loop
          
          cost.profile.mean <- apply(cost.profile,2,mean)
          cost.best <- costVal[which.max(cost.profile.mean)]
        } #End cost tuning
        
        bestModel <- svm(trainLabels~.,data=train,kernel='linear',cost=cost.best)
        
        classif.svm[[roiList[roi]]]$cost[run] <- cost.best
        fittedValuesTrain <- bestModel$decision.values
        
        #Extract svm weights, for curiosities sake ...
        w <- abs(as.numeric(t(bestModel$coefs)%*%bestModel$SV))
        classif.svm[[roiList[roi]]]$svmWeight[[run]] <- w
        
        #Generate test predictions
        predClassTest <- predict(bestModel, test[,-(grep('testLabels',names(test)))])
        
        classif.svm[[roiList[roi]]]$predClass[[run]] <- predClassTest
        
        fittedValuesTest <- attributes(predict(bestModel, test[,-(grep('testLabels',names(test)))], 
                                               decision.values=TRUE))$decision.values
        
        classif.svm[[roiList[roi]]]$fittedValue[[run]] <- fittedValuesTest
        
        #Record classification accuracy for this fold
        tree <- table(predict=predClassTest,truth=test$testLabels)
        ca.fold <- sum(diag(tree))/sum(tree)
        classif.svm[[roiList[roi]]]$foldAccuracy[run] <- ca.fold
        rm(bestModel,train, test, trainLabels,testLabels,fittedValuesTrain)
        
      }
      
      #Record mean accuracy
      overallAcc <- sum(unlist(classif.svm[[roi]]$predClass) == unlist(classif.svm[[roi]]$trueClass))/length(unlist(classif.svm[[roi]]$predClass))
      classif.svm[[roiList[roi]]]$meanAccuracy <- overallAcc
      print(paste0('     Mean CV classification accuracy: ',round(overallAcc,2)))
      
      #Record area under ROC curve (AUC)
      predob <- prediction(unlist(classif.svm[[roi]]$fittedValue),unlist(classif.svm[[roi]]$trueClass),
                           label.ordering=rev(levels(betas.roi[[labToClassify]])))
      perf <- performance(predob,measure='auc')
      auc <- perf@y.values[[1]]
      classif.svm[[roiList[roi]]]$AUC <- auc
      # print(paste0('     Area under the ROC curve: ',round(auc,2)))
      
    }
    
  setwd(baseDir)
  return(classif.svm)
}
