# clears workspace:  
rm(list=ls()) 

setwd('C:/Users/Jeff/Documents/R/LoLComp/upload')

#load("splitsig.RData")

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)

source("BootstrapMean95.R")
source("PlotFunc_horizontal.R")

# Load data --------------------------------------------------------

#OnlineData <- read_csv("OnlineData.csv",col_names=TRUE)


## data for removing type 1 disagreement
OnlineData <- read_csv("OnlineData_psi.csv",col_names=TRUE)

nExperts = max(OnlineData$Participant)
nItems = 12

SONAData <- read_csv("SONAData.csv",col_names=TRUE)

nNovices = max(SONAData$Participant)

totalNovice = nItems*nNovices

# Prep data for raw JAGS inputs -----------------------------------------------
ConditionSize = length(which(OnlineData$SpeakerEarly==1))+length(which(SONAData$SpeakerEarly==1))
NoviceConditionSize = length(which(SONAData$SpeakerEarly==1))
ExpertConditionSize = length(which(OnlineData$SpeakerEarly==1))
ExpertStart = NoviceConditionSize + 1

NoviceIDs <- rep(1:nNovices,each=6)
IDStart = 1+nNovices
IDEnd = nNovices+nExperts
ExpertIDs <- rep(IDStart:IDEnd,each=6)

nParticipants = nExperts+nNovices
nTrials = 6
nExpResponses = nExperts*nTrials

## adjustment for removing type 1 disagreement
nBad_early <- length(which(OnlineData$Early_psi_bern>.89999))
nBad_late <- length(which(OnlineData$Late_psi_bern>.89999))
nExpResponses_early <- nExpResponses-nBad_early
nExpResponses_late <- nExpResponses-nBad_late

nNovResponses = nNovices*nTrials
nTotalResponses_early = nNovResponses+nExpResponses_early
nTotalResponses_late = nNovResponses+nExpResponses_late
nBeta = 7

XyesEarly <- matrix(data=rep(0,nTotalResponses_early*(nBeta+1)),nrow = nTotalResponses_early,ncol = (nBeta+1))
YyesEarly <- rep(NaN,nTotalResponses_early)
zyesEarly <- matrix(data=rep(0,nTotalResponses_early*2),nrow = nTotalResponses_early,ncol = 2)

XyesLate <- matrix(data=rep(0,nTotalResponses_late*(nBeta+1)),nrow = nTotalResponses_late,ncol = (nBeta+1))
YyesLate <- rep(NaN,nTotalResponses_late)
zyesLate <- matrix(data=rep(0,nTotalResponses_late*2),nrow = nTotalResponses_late,ncol = 2)

XyesEarly[(nExpResponses_early+1):(nTotalResponses_early),1] <- 1
zyesEarly[(nExpResponses_early+1):(nTotalResponses_early),2] <- 1
zyesEarly[1:nExpResponses_early,1] <- 1

XyesLate[(nExpResponses_late+1):(nTotalResponses_late),1] <- 1
zyesLate[(nExpResponses_late+1):(nTotalResponses_late),2] <- 1
zyesLate[1:nExpResponses_late,1] <- 1

tickeryesEarly = 1
tickeryesLate = 1

for(i in 1:nrow(OnlineData)){
  if(OnlineData$Early_psi_bern[i]<.9){
    if(OnlineData$SpeakerEarly[i]==1){
      
      YyesEarly[tickeryesEarly] <- OnlineData$EarlyResponse[i]
      
      XyesEarly[tickeryesEarly,(OnlineData$CompType[i]+1)] <- 1
      XyesEarly[tickeryesEarly,(nBeta+1)] <- OnlineData$Participant[i]
      
      tickeryesEarly = tickeryesEarly+1
    }
  }
  if(OnlineData$Late_psi_bern[i]<.9){  
    if(OnlineData$SpeakerLate[i]==1){
      
      YyesLate[tickeryesLate] <- OnlineData$LateResponse[i]
      
      XyesLate[tickeryesLate,(OnlineData$CompType[i]+1)] <- 1
      XyesLate[tickeryesLate,(nBeta+1)] <- OnlineData$Participant[i]
      
      tickeryesLate = tickeryesLate+1
    }
  }
}
  
for(i in 1:nrow(SONAData)){
    
    if(SONAData$SpeakerEarly[i]==1){
      
      YyesEarly[tickeryesEarly] <- SONAData$EarlyResponse[i]
      XyesEarly[tickeryesEarly,(nBeta+1)] <- SONAData$Participant[i]+49
      
      #XyesEarly[tickeryesEarly,(SONAData$CompType[i]+1)] <- 1
      
      tickeryesEarly = tickeryesEarly+1
    }
    
    if(SONAData$SpeakerLate[i]==1){
      
      YyesLate[tickeryesLate] <- SONAData$LateResponse[i]
      XyesLate[tickeryesLate,(nBeta+1)] <- SONAData$Participant[i]+49
      
      #XyesLate[tickeryesLate,(SONAData$CompType[i]+1)] <- 1
      
      tickeryesLate = tickeryesLate+1
    }
}

# MVN Precision Matrix (for raw JAGS betas) ----------------------------------------
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Run raw linear regression in JAGS -------------------------------------------
nSamples = 4000

rawBetayesEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
rawDIC <- rep(NaN,2)

# X <- XyesEarly
# Y<- YyesEarly
# z <- zyesEarly
# 
# data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta") # to be passed on to JAGS
# # parameters to be monitored:	
# parameters <- c("BETA","tau","sigma","alpha","Ypred")
# 
# #Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# myinits <- list(list("BETA" = runif(7,-.5,.5)))
# 
# # The following command calls JAGS with specific options.
# #This is for running 1 chain (use code below for faster multiple chains)
# samples <- jags(data, #inits=myinits, 
#                 parameters.to.save=parameters,
#                 model.file="LoLComp_LinearReg_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

## model for removing type 1 disagreement
X <- XyesEarly[,1:nBeta]
Y<- YyesEarly
Participant_vector <- XyesEarly[,(nBeta+1)]
nData <- nrow(XyesEarly)
nParticipants <- max(Participant_vector)

data <- list("nData","nParticipants","X", "Y","bPrec","nBeta","Participant_vector") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg__latentcheck_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

rawBetayesEarlysamples <- samples$BUGSoutput$sims.list$BETA

rawsigma_yesEarlysamples <- samples$BUGSoutput$sims.list$sigma

rawtau_yesEarlysamples <- samples$BUGSoutput$sims.list$tau

rawsummaryyesEarly <- samples$BUGSoutput$summary

rawDIC[1] <- samples$BUGSoutput$DIC

##Yes Late Regression##

rawBetayesLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

# X <- XyesLate
# Y<- YyesLate
# z <- zyesLate
# 
# data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta") # to be passed on to JAGS
# # parameters to be monitored:	
# parameters <- c("BETA","tau","sigma","alpha","Ypred")
# 
# #Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# myinits <- list(list("BETA" = runif(7,-.5,.5)))
# 
# # The following command calls JAGS with specific options.
# #This is for running 1 chain (use code below for faster multiple chains)
# samples <- jags(data, #inits=myinits, 
#                 parameters.to.save=parameters,
#                 model.file="LoLComp_LinearReg_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

## model for removing type 1 disagreement
X <- XyesLate[,1:nBeta]
Y<- YyesLate
Participant_vector <- XyesLate[,(nBeta+1)]
nData <- nrow(XyesLate)
nParticipants <- max(Participant_vector)

data <- list("nData","nParticipants","X", "Y","bPrec","nBeta","Participant_vector") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg__latentcheck_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

rawBetayesLatesamples <- samples$BUGSoutput$sims.list$BETA

rawsigma_yesLatesamples <- samples$BUGSoutput$sims.list$sigma

rawtau_yesLatesamples <- samples$BUGSoutput$sims.list$tau

rawsummaryyesLate <- samples$BUGSoutput$summary

rawDIC[2] <- samples$BUGSoutput$DIC

# Run split sigma linear regression in JAGS -------------------------------------------
nSamples = 4000

splitsig_BetayesEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
splitsig_DIC <- rep(NaN,2)

X <- XyesEarly[,1:nBeta]
Y<- YyesEarly
Participant_vector <- XyesEarly[,(nBeta+1)]
nData <- nrow(XyesEarly)
nParticipants <- max(Participant_vector)
z <- zyesEarly

data <- list("nData","nParticipants","X", "Y","bPrec","nBeta","Participant_vector","z") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_nov","sigma_exp","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_latentcheck_jags_splitsig.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

splitsig_BetayesEarlysamples <- samples$BUGSoutput$sims.list$BETA

sigmanov_yesEarlysamples <- samples$BUGSoutput$sims.list$sigma_nov

sigmaexp_yesEarlysamples <- samples$BUGSoutput$sims.list$sigma_exp

tau_yesEarlysamples <- samples$BUGSoutput$sims.list$tau

splitsig_summaryyesEarly <- samples$BUGSoutput$summary

splitsig_DIC[1] <- samples$BUGSoutput$DIC

##Yes Late Regression##

splitsig_BetayesLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XyesLate[,1:nBeta]
Y<- YyesLate
Participant_vector <- XyesLate[,(nBeta+1)]
nData <- nrow(XyesLate)
nParticipants <- max(Participant_vector)
z <- zyesLate

data <- list("nData","nParticipants","X", "Y","bPrec","nBeta","Participant_vector","z") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_nov","sigma_exp","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_latentcheck_jags_splitsig.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

splitsig_BetayesLatesamples <- samples$BUGSoutput$sims.list$BETA

sigmanov_yesLatesamples <- samples$BUGSoutput$sims.list$sigma_nov

sigmaexp_yesLatesamples <- samples$BUGSoutput$sims.list$sigma_exp

tau_yesLatesamples <- samples$BUGSoutput$sims.list$tau

splitsig_summaryyesLate <- samples$BUGSoutput$summary

splitsig_DIC[2] <- samples$BUGSoutput$DIC

# Analyze sigmas and taus from standard model -----------------------------
RhoTest_yesEarly <- rep(NaN,length(rawsigma_yesEarlysamples))
RhoTest_yesLate <- matrix(NaN,length(rawsigma_yesLatesamples))

for(i in 1:length(rawsigma_yesEarlysamples)){
  RhoTest_yesEarly[i] <- rawtau_yesEarlysamples[i]^2/
    (rawtau_yesEarlysamples[i]^2 + rawsigma_yesEarlysamples[i]^2)
  
  RhoTest_yesLate[i] <- rawtau_yesLatesamples[i]^2/
    (rawtau_yesLatesamples[i]^2 + rawsigma_yesLatesamples[i]^2)
}

ind975 = length(rawsigma_yesEarlysamples)*.975
ind075 = length(rawsigma_yesEarlysamples)*.075

sortRhoTest_yesEarly <- sort(RhoTest_yesEarly)
sortRhoTest_yesLate <- sort(RhoTest_yesLate)

mean(RhoTest_yesEarly)
sortRhoTest_yesEarly[ind975]
sortRhoTest_yesEarly[ind075]

mean(RhoTest_yesLate)
sortRhoTest_yesLate[ind975]
sortRhoTest_yesLate[ind075]

# Compare expert and novice sigmas ----------------------------------------

SigmaTest_yesEarly <- matrix(NaN,nrow = length(sigmanov_yesEarlysamples),ncol = 2)
SigmaTest_yesLate <- matrix(NaN,nrow = length(sigmanov_yesLatesamples),ncol = 2)

SigmaDiff_yesEarly <- rep(NaN,length(sigmanov_yesEarlysamples))
SigmaDiff_yesLate <- rep(NaN,length(sigmanov_yesEarlysamples))

for(i in 1:length(sigmanov_yesLatesamples)){
  SigmaTest_yesEarly[i,1] <- tau_yesEarlysamples[i]^2/
    (tau_yesEarlysamples[i]^2 + sigmaexp_yesEarlysamples[i]^2)
  SigmaTest_yesEarly[i,2] <- tau_yesEarlysamples[i]^2/
    (tau_yesEarlysamples[i]^2 + sigmanov_yesEarlysamples[i]^2)
  
  SigmaDiff_yesEarly[i] <- sigmaexp_yesEarlysamples[i]-sigmanov_yesEarlysamples[i]
  SigmaDiff_yesLate[i] <- sigmaexp_yesLatesamples[i]-sigmanov_yesLatesamples[i]
  
  SigmaTest_yesLate[i,1] <- tau_yesLatesamples[i]^2/
    (tau_yesLatesamples[i]^2 + sigmaexp_yesLatesamples[i]^2)
  SigmaTest_yesLate[i,2] <- tau_yesEarlysamples[i]^2/
    (tau_yesLatesamples[i]^2 + sigmanov_yesLatesamples[i]^2)
}

sortSigmaDiff_yesEarly <- sort(SigmaDiff_yesEarly)
sortSigmaDiff_yesLate <- sort(SigmaDiff_yesLate)

ind975 = length(sortSigmaDiff_yesEarly)*.975
ind075 = length(sortSigmaDiff_yesEarly)*.075

mean(SigmaDiff_yesEarly)
sortSigmaDiff_yesEarly[ind975]
sortSigmaDiff_yesEarly[ind075]

mean(SigmaDiff_yesLate)
sortSigmaDiff_yesLate[ind975]
sortSigmaDiff_yesLate[ind075]

sortSigmaTest_yesEarly <- matrix(NaN,nrow = length(sigmanov_yesEarlysamples),ncol = 2)
sortSigmaTest_yesLate <- matrix(NaN,nrow = length(sigmanov_yesLatesamples),ncol = 2)

sortSigmaTest_yesEarly[,1] <- sort(SigmaTest_yesEarly[,1])
sortSigmaTest_yesEarly[,2] <- sort(SigmaTest_yesEarly[,2])
sortSigmaTest_yesLate[,1] <- sort(SigmaTest_yesLate[,1])
sortSigmaTest_yesLate[,2] <- sort(SigmaTest_yesLate[,2])

mean(SigmaTest_yesEarly[,1])
sortSigmaTest_yesEarly[ind975,1]
sortSigmaTest_yesEarly[ind075,1]

mean(SigmaTest_yesEarly[,2])
sortSigmaTest_yesEarly[ind975,2]
sortSigmaTest_yesEarly[ind075,2]

mean(SigmaTest_yesLate[,1])
sortSigmaTest_yesLate[ind975,1]
sortSigmaTest_yesLate[ind075,1]

mean(SigmaTest_yesLate[,2])
sortSigmaTest_yesLate[ind975,2]
sortSigmaTest_yesLate[ind075,2]

# Savage Dickey for sigma difference --------------------------------------

SigmaSimPulls1 <- rgamma(1000000, shape = 1.5, rate = 2)
SigmaSimPulls2 <- rgamma(1000000, shape = 1.5, rate = 2)

SigmaSimDiff <- SigmaSimPulls1-SigmaSimPulls2

hist(SigmaSimDiff)

SigmaDiffPrior_mean = mean(SigmaSimDiff)
SigmaDiffPrior_sd = sd(SigmaSimDiff)

prior_density_SigDiff = dnorm(0,0,.867)

posterior_density_sigDiff_yesEarly <- dlogspline(0,logspline(SigmaDiff_yesEarly))
posterior_density_sigDiff_yesLate <- dlogspline(0,logspline(SigmaDiff_yesLate))

SavDick_sig <- c(posterior_density_sigDiff_yesEarly/prior_density_SigDiff,
                 posterior_density_sigDiff_yesLate/prior_density_SigDiff)

# Plot beta inferences from raw -------------------------------------------
rawPlotValues_EarlyYes <- matrix(NaN, nrow = 7,ncol = 3)
rawPlotValues_LateYes <- matrix(NaN, nrow = 7,ncol = 3)

loconf = nrow(rawBetayesEarlysamples) * .025
hiconf = nrow(rawBetayesEarlysamples) * .975

tempsort <- sort(rawBetayesEarlysamples[,1])
rawPlotValues_EarlyYes[1,1] <- mean(tempsort[1500:1001])
rawPlotValues_EarlyYes[1,2] <- tempsort[loconf]
rawPlotValues_EarlyYes[1,3] <- tempsort[hiconf]

tempsort <- sort(rawBetayesLatesamples[,1])
rawPlotValues_LateYes[1,1] <- mean(tempsort[1500:1001])
rawPlotValues_LateYes[1,2] <- tempsort[loconf]
rawPlotValues_LateYes[1,3] <- tempsort[hiconf]

for(i in 2:7){
  
  tempsort <- sort(rawBetayesEarlysamples[,i])
  rawPlotValues_EarlyYes[i,1] <- mean(tempsort[1500:1001])
  rawPlotValues_EarlyYes[i,2] <- tempsort[loconf]
  rawPlotValues_EarlyYes[i,3] <- tempsort[hiconf]
  
  
  tempsort <- sort(rawBetayesLatesamples[,i])
  rawPlotValues_LateYes[i,1] <- mean(tempsort[1500:1001])
  rawPlotValues_LateYes[i,2] <- tempsort[loconf]
  rawPlotValues_LateYes[i,3] <- tempsort[hiconf]
}

plot(1:6,rawPlotValues_EarlyYes[2:7,1],ylim = c(.15,.9),main = "Yes Early",xlab = 'beta',ylab = 'Interpretation')
abline(h = rawPlotValues_EarlyYes[1,1],col = '#DC143C')
abline(h= rawPlotValues_EarlyYes[1,2],col = '#DC143C',lty = 2)
abline(h= rawPlotValues_EarlyYes[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = rawPlotValues_EarlyYes[i+1,2],
         x1 = i, y1 = rawPlotValues_EarlyYes[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}
legend(3.5,.45,legend = c('Novice Interpretation','Expert Interpretation'),pch = c(NA,1),lty = c(1,NA),
       col = c('#DC143C','black'),cex = .8)

plot(1:6,rawPlotValues_LateYes[2:7,1],ylim = c(.15,.9),main = "Yes Late",xlab = 'beta',ylab = 'Interpretation')
abline(h = rawPlotValues_LateYes[1,1],col = '#DC143C')
abline(h= rawPlotValues_LateYes[1,2],col = '#DC143C',lty = 2)
abline(h= rawPlotValues_LateYes[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = rawPlotValues_LateYes[i+1,2],
         x1 = i, y1 = rawPlotValues_LateYes[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}

# Test for expected order in "Yes" conditions -----------------------------

OrderTestEarly <- matrix(NaN,nrow=nrow(rawBetayesEarlysamples),ncol = 7)

OrderTestEarly[,1:6] <- rawBetayesEarlysamples[,2:7]

##Early hypothesis: 1 > 2 > 3, 6 > 5 > 4
for(i in 1:nrow(rawBetayesEarlysamples)){
  firstcheck <- OrderTestEarly[i,1] > OrderTestEarly[i,2] && 
    OrderTestEarly[i,2] > OrderTestEarly[i,3]
  secondcheck <- OrderTestEarly[i,6] > OrderTestEarly[i,5] && 
    OrderTestEarly[i,5] > OrderTestEarly[i,4]
  OrderTestEarly[i,7] <- firstcheck && secondcheck
}

OrderTestLate <- matrix(NaN,nrow=nrow(rawBetayesLatesamples),ncol = 7)

OrderTestLate[,1:6] <- rawBetayesLatesamples[,2:7]

##Late hypothesis: 3 > 2 > 1, 4 > 5 > 6
for(i in 1:nrow(rawBetayesLatesamples)){
  firstcheck <- OrderTestLate[i,3] > OrderTestLate[i,2] && 
    OrderTestLate[i,2] > OrderTestLate[i,1]
  secondcheck <- OrderTestLate[i,4] > OrderTestLate[i,5] && 
    OrderTestLate[i,5] > OrderTestLate[i,6]
  OrderTestLate[i,7] <- firstcheck && secondcheck
}

EarlyOrderProb <- mean(OrderTestEarly[,7])
LateOrderProb <- mean(OrderTestLate[,7])

EarlyOrderOdds <- EarlyOrderProb/(1-EarlyOrderProb)
LateOrderOdds <- LateOrderProb/(1-LateOrderProb)

PriorOdds <- (1/36)/(1-(1/36))

BFEarlyOrder <- EarlyOrderOdds/PriorOdds
BFLateOrder <- LateOrderOdds/PriorOdds

# Prep data for overfull JAGS inputs -----------------------------------------------
nBeta = 12

XyesEarly_overfull <- matrix(data=rep(0,nTotalResponses_early*(nBeta+1)),nrow = nTotalResponses_early,ncol = (nBeta+1))
YyesEarly <- rep(NaN,nTotalResponses_early)

XyesLate_overfull <- matrix(data=rep(0,nTotalResponses_late*(nBeta+1)),nrow = nTotalResponses_late,ncol = (nBeta+1))
YyesLate <- rep(NaN,nTotalResponses_late)

tickeryesEarly = 1
tickeryesLate = 1

for(i in 1:nrow(OnlineData)){
  if(OnlineData$Early_psi_bern[i]<.9){
    if(OnlineData$SpeakerEarly[i]==1){
      
      YyesEarly[tickeryesEarly] <- OnlineData$EarlyResponse[i]
      
      XyesEarly_overfull[tickeryesEarly,(OnlineData$CompType[i]+1)] <- 1
      XyesEarly_overfull[tickeryesEarly,(nBeta+1)] <- OnlineData$Participant[i]
      
      tickeryesEarly = tickeryesEarly+1
    }
  }
  if(OnlineData$Late_psi_bern[i]<.9){  
    if(OnlineData$SpeakerLate[i]==1){
      
      YyesLate[tickeryesLate] <- OnlineData$LateResponse[i]
      
      XyesLate_overfull[tickeryesLate,(OnlineData$CompType[i]+1)] <- 1
      XyesLate_overfull[tickeryesLate,(nBeta+1)] <- OnlineData$Participant[i]
      
      tickeryesLate = tickeryesLate+1
    }
  }
}

for(i in 1:nrow(SONAData)){
  
  
  if(SONAData$SpeakerEarly[i]==1){
    
    YyesEarly[tickeryesEarly] <- SONAData$EarlyResponse[i]
    
    XyesEarly_overfull[tickeryesEarly,(SONAData$CompType[i]+6)] <- 1
    XyesEarly_overfull[tickeryesEarly,(nBeta+1)] <- SONAData$Participant[i]+49
    
    tickeryesEarly = tickeryesEarly+1
  }
  
  if(SONAData$SpeakerLate[i]==1){
    
    YyesLate[tickeryesLate] <- SONAData$LateResponse[i]
    
    XyesLate_overfull[tickeryesLate,(SONAData$CompType[i]+6)] <- 1
    XyesLate_overfull[tickeryesLate,(nBeta+1)] <- SONAData$Participant[i]+49
    
    tickeryesLate = tickeryesLate+1
  }
}

# MVN Precision Matrix (for overfull JAGS betas) ----------------------------------------
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Run overfull linear regression in JAGS -------------------------------------------
nSamples = 4000

overfullBetayesEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
overfullDIC <- rep(NaN,2)

X <- XyesEarly_overfull[,1:nBeta]
Y<- YyesEarly
Participant_vector <- XyesEarly_overfull[,(nBeta+1)]
nData <- nrow(XyesEarly_overfull)
nParticipants <- max(Participant_vector)

data <- list("nData","nParticipants","X", "Y","bPrec","nBeta","Participant_vector") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg__latentcheck_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

overfullBetayesEarlysamples <- samples$BUGSoutput$sims.list$BETA

overfullsummaryyesEarly <- samples$BUGSoutput$summary

overfullDIC[1] <- samples$BUGSoutput$DIC

##Yes Late Regression##

overfullBetayesLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XyesLate_overfull[,1:nBeta]
Y<- YyesLate
Participant_vector <- XyesLate_overfull[,(nBeta+1)]
nData <- nrow(XyesLate_overfull)
nParticipants <- max(Participant_vector)

data <- list("nData","nParticipants","X", "Y","bPrec","nBeta","Participant_vector") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg__latentcheck_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

overfullBetayesLatesamples <- samples$BUGSoutput$sims.list$BETA

overfullsummaryyesLate <- samples$BUGSoutput$summary

overfullDIC[2] <- samples$BUGSoutput$DIC

# Compare DIC's (between main and overfull) -----------------------------------------------------------

par(mfrow=c(1,1))

plot(1:2,rawDIC,ylim = c(-900,-300),xlab = 'Condition',ylab = 'DIC', main = "Main Model and Overfull DIC Comparison")

points(1:2,overfullDIC,col = '#DC143C',pch=3)

legend(1, -400, legend=c("Main Model", "Overfull"),
       col=c("black", "red"), pch=c(1,3), cex=0.8)

text_rawDIC <- matrix(NaN,nrow = 2, ncol =2)
text_rawDIC[,1] <- 1:2
text_rawDIC[,2] <- round(rawDIC,0)
text_rawDIC <- data.frame(text_rawDIC)
names(text_rawDIC) <- c("Condition","DIC")

text(DIC~Condition, labels=DIC,data=text_rawDIC, cex=0.9, font=1,pos=1)

text_overfullDIC <- matrix(NaN,nrow = 2, ncol =2)
text_overfullDIC[,1] <- 1:2
text_overfullDIC[,2] <- round(overfullDIC,0)
text_overfullDIC <- data.frame(text_overfullDIC)
names(text_overfullDIC) <- c("Condition","DIC")

text(DIC~Condition, labels=DIC,data=text_overfullDIC, cex=0.9, font=1,pos=3,col='#DC143C')

# Expected order test for experts in "Yes" condition, overfull model ----------------------

OrderTestEarly_expertoverfull <- matrix(NaN,nrow=nrow(overfullBetayesEarlysamples),ncol = 7)

OrderTestEarly_expertoverfull[,1:6] <- overfullBetayesEarlysamples[,1:6]

##Early hypothesis: 1 > 2 > 3, 6 > 5 > 4
for(i in 1:nrow(overfullBetayesEarlysamples)){
  firstcheck <- OrderTestEarly_expertoverfull[i,1] > OrderTestEarly_expertoverfull[i,2] && 
    OrderTestEarly_expertoverfull[i,2] > OrderTestEarly_expertoverfull[i,3]
  secondcheck <- OrderTestEarly_expertoverfull[i,6] > OrderTestEarly_expertoverfull[i,5] && 
    OrderTestEarly_expertoverfull[i,5] > OrderTestEarly_expertoverfull[i,4]
  OrderTestEarly_expertoverfull[i,7] <- firstcheck && secondcheck
}

OrderTestLate_expertoverfull <- matrix(NaN,nrow=nrow(overfullBetayesLatesamples),ncol = 7)

OrderTestLate_expertoverfull[,1:6] <- overfullBetayesLatesamples[,1:6]

##Late hypothesis: 3 > 2 > 1, 4 > 5 > 6
for(i in 1:nrow(overfullBetayesLatesamples)){
  firstcheck <- OrderTestLate_expertoverfull[i,3] > OrderTestLate_expertoverfull[i,2] && 
    OrderTestLate_expertoverfull[i,2] > OrderTestLate_expertoverfull[i,1]
  secondcheck <- OrderTestLate_expertoverfull[i,4] > OrderTestLate_expertoverfull[i,5] && 
    OrderTestLate_expertoverfull[i,5] > OrderTestLate_expertoverfull[i,6]
  OrderTestLate_expertoverfull[i,7] <- firstcheck && secondcheck
}

EarlyOrderProb_expertoverfull <- mean(OrderTestEarly_expertoverfull[,7])
LateOrderProb_expertoverfull <- mean(OrderTestLate_expertoverfull[,7])

EarlyOrderOdds_expertoverfull <- EarlyOrderProb_expertoverfull/(1-EarlyOrderProb_expertoverfull)
LateOrderOdds_expertoverfull <- LateOrderProb_expertoverfull/(1-LateOrderProb_expertoverfull)

PriorOdds_expertoverfull <- (1/36)/(1-(1/36))

BFEarlyOrder_expertoverfull <- EarlyOrderOdds_expertoverfull/PriorOdds_expertoverfull
BFLateOrder_expertoverfull <- LateOrderOdds_expertoverfull/PriorOdds_expertoverfull

# Expected order test for novices in "Yes" condition, overfull model ----------------------

OrderTestEarly_noviceoverfull <- matrix(NaN,nrow=nrow(overfullBetayesEarlysamples),ncol = 7)

OrderTestEarly_noviceoverfull[,1:6] <- overfullBetayesEarlysamples[,7:12]

##Early hypothesis: 1 > 2 > 3, 6 > 5 > 4
for(i in 1:nrow(overfullBetayesEarlysamples)){
  firstcheck <- OrderTestEarly_noviceoverfull[i,1] > OrderTestEarly_noviceoverfull[i,2] && 
    OrderTestEarly_noviceoverfull[i,2] > OrderTestEarly_noviceoverfull[i,3]
  secondcheck <- OrderTestEarly_noviceoverfull[i,6] > OrderTestEarly_noviceoverfull[i,5] && 
    OrderTestEarly_noviceoverfull[i,5] > OrderTestEarly_noviceoverfull[i,4]
  OrderTestEarly_noviceoverfull[i,7] <- firstcheck && secondcheck
}

OrderTestLate_noviceoverfull <- matrix(NaN,nrow=nrow(overfullBetayesLatesamples),ncol = 7)

OrderTestLate_noviceoverfull[,1:6] <- overfullBetayesLatesamples[,7:12]

##Late hypothesis: 3 > 2 > 1, 4 > 5 > 6
for(i in 1:nrow(overfullBetayesLatesamples)){
  firstcheck <- OrderTestLate_noviceoverfull[i,3] > OrderTestLate_noviceoverfull[i,2] && 
    OrderTestLate_noviceoverfull[i,2] > OrderTestLate_noviceoverfull[i,1]
  secondcheck <- OrderTestLate_noviceoverfull[i,4] > OrderTestLate_noviceoverfull[i,5] && 
    OrderTestLate_noviceoverfull[i,5] > OrderTestLate_noviceoverfull[i,6]
  OrderTestLate_noviceoverfull[i,7] <- firstcheck && secondcheck
}

EarlyOrderProb_noviceoverfull <- mean(OrderTestEarly_noviceoverfull[,7])
LateOrderProb_noviceoverfull <- mean(OrderTestLate_noviceoverfull[,7])

EarlyOrderOdds_noviceoverfull <- EarlyOrderProb_noviceoverfull/(1-EarlyOrderProb_noviceoverfull)
LateOrderOdds_noviceoverfull <- LateOrderProb_noviceoverfull/(1-LateOrderProb_noviceoverfull)

PriorOdds_noviceoverfull <- (1/36)/(1-(1/36))

BFEarlyOrder_noviceoverfull <- EarlyOrderOdds_noviceoverfull/PriorOdds_noviceoverfull
BFLateOrder_noviceoverfull <- LateOrderOdds_noviceoverfull/PriorOdds_noviceoverfull

# Posterior predictive accuracy by comp type ------------------------------
yesEarly_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
noEarly_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
yesLate_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
noLate_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)



# Plot beta inferences from raw w/ actual data -------------------------------------------
par(mfrow=c(1,1))
par(mar = c(3.5, 3, 1, 2))

yesEarly_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
noEarly_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
yesLate_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
noLate_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
testTotal = nExpResponses+nNovResponses

for(i in 1:6){
  tempindex <- which(XyesEarly[1:nExpResponses,(i+1)]==1)
  yesEarly_bycomp_expert[1,i] = mean(YyesEarly[tempindex])
  yesEarly_bycomp_expert[2,i] = sd(YyesEarly[tempindex])
  yesEarly_bycomp_expert[3,i] = mean(YyesEarly[(nExpResponses+1):testTotal])
  yesEarly_bycomp_expert[4,i] = sd(YyesEarly[(nExpResponses+1):testTotal])
  tempindex_novice <- SONAData[which(SONAData$CompType==i),]
  yesEarly_bycomp_expert[5,i] = 
    mean(tempindex_novice$EarlyResponse[which(tempindex_novice$SpeakerEarly==1)])
}

for(i in 1:6){
  tempindex = which(XyesLate[1:nExpResponses,i+1]==1)
  yesLate_bycomp_expert[1,i] = mean(YyesLate[tempindex])
  yesLate_bycomp_expert[2,i] = sd(YyesLate[tempindex])
  yesLate_bycomp_expert[3,i] = mean(YyesLate[(nExpResponses+1):testTotal])
  yesLate_bycomp_expert[4,i] = sd(YyesLate[(nExpResponses+1):testTotal])
  tempindex_novice <- SONAData[which(SONAData$CompType==i),]
  yesLate_bycomp_expert[5,i] = 
    mean(tempindex_novice$LateResponse[which(tempindex_novice$SpeakerLate==1)])
}

ConfInt <- function(Dev,n){
  ConfInt = 1.96*(Dev/sqrt(n))
  return(ConfInt)
}

plot1 = MainPlot(rawPlotValues_EarlyYes[1,1:3],yesEarly_bycomp_expert[3,1],yesEarly_bycomp_expert[5,1:6],
                 rawPlotValues_EarlyYes[2:7,1:3],yesEarly_bycomp_expert[1,1:6],'',
                 if_legend = 1,leftmost = 1)
par(mar = c(3.5, 3, 1, 2))
plot2 = MainPlot(rawPlotValues_LateYes[1,1:3],yesLate_bycomp_expert[3,1],yesLate_bycomp_expert[5,1:6],
                 rawPlotValues_LateYes[2:7,1:3],yesLate_bycomp_expert[1,1:6],'',
                 if_legend = 0, leftmost = 0)

# Hypothesis test using Savage-Dickey -------------------------------------

prior_density = dnorm(0,0,.3537)
posterior_density_yesEarly <- rep(NaN,6)
posterior_density_yesLate <- rep(NaN,6)

posterior_12B_density_yesEarly <- rep(NaN,6)
posterior_12B_density_yesLate <- rep(NaN,6)

for(i in 1:6){
  tempdif <- rawBetayesEarlysamples[,1]-rawBetayesEarlysamples[,i+1]
  posterior_density_yesEarly[i] <- dlogspline(0,logspline(tempdif))
  tempdif <- overfullBetayesEarlysamples[,i+6]-overfullBetayesEarlysamples[,i]
  posterior_12B_density_yesEarly[i] <- dlogspline(0,logspline(tempdif))
  
  tempdif <- rawBetayesLatesamples[,1]-rawBetayesLatesamples[,i+1]
  posterior_density_yesLate[i] <- dlogspline(0,logspline(tempdif))
  tempdif <- overfullBetayesLatesamples[,i+6]-overfullBetayesLatesamples[,i]
  posterior_12B_density_yesLate[i] <- dlogspline(0,logspline(tempdif))
}

SavDick <- data.frame(matrix(NaN,nrow=6,ncol=4))
names(SavDick) <- c("Yes Early","12B Yes Early","Yes Late","12B Yes Late")

SavDick[,1] <- posterior_density_yesEarly/prior_density
SavDick[,2] <- posterior_12B_density_yesEarly/prior_density
SavDick[,3] <- posterior_density_yesLate/prior_density
SavDick[,4] <- posterior_12B_density_yesLate/prior_density

# Expert TVJ by comp type -------------------------------------

par(mfrow=c(1,2))
EarlyGame_Ill = c(.70,.63,.55,.55,.63,.70)
LateGame_Ill = c(.55,.63,.70,.70,.63,.55)
par(mar = c(3, 4, 2, 2))

TVJ_Early_bycomp <- c(mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==1)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==2)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==3)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==4)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==5)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==6)]))

TVJ_Late_bycomp <- c(mean(OnlineData$LateTVJ[which(OnlineData$CompType==1)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==2)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==3)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==4)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==5)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==6)]))

TVJ_Early_bycomp_boot = matrix(NaN, 6, 2)
TVJ_Late_bycomp_boot = matrix(NaN, 6, 2)

for(i in 1:6){
  nIterations = 1000
  nSamples = 100
  dataVector_EarlyTVJ = OnlineData$EarlyTVJ[which(OnlineData$CompType==i)]
  bootRange_EarlyTVJ = BootstrapMean95(dataVector_EarlyTVJ,nIterations = nIterations,nSamples = nSamples)
  TVJ_Early_bycomp_boot[i,1] <- bootRange_EarlyTVJ[1]
  TVJ_Early_bycomp_boot[i,2] <- bootRange_EarlyTVJ[2]
  
  dataVector_LateTVJ = OnlineData$LateTVJ[which(OnlineData$CompType==i)]
  bootRange_LateTVJ = BootstrapMean95(dataVector_LateTVJ,nIterations = nIterations,nSamples = nSamples)
  TVJ_Late_bycomp_boot[i,1] <- bootRange_LateTVJ[1]
  TVJ_Late_bycomp_boot[i,2] <- bootRange_LateTVJ[2]
}  

plot(1:6,TVJ_Early_bycomp,ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'black',bty = 'n',xlim = c(.75,6),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 23)

for(i in 1:6){
  arrows(x0 = i, y0 = TVJ_Early_bycomp_boot[i,1],
         x1 = i, y1 = TVJ_Early_bycomp_boot[i,2],code = 3,
         col = 'black',angle = 90, length = .1,lty = 1)
}

xlabels = c('E+','E0','E-','L+','L0','L-')
xat = c(1:6)

title(main="Early Game ", line=.5, cex.main=1.8)
title(ylab="Endorsement Rate", line=2.1, cex.lab=1.8)
title(xlab="Composition Type", line=1.6, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,6.25),y=c(0,0),lwd =2)
abline(v=0)

plot(1:6,TVJ_Late_bycomp,ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'black',bty = 'n',xlim = c(.75,6),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 23)

for(i in 1:6){
  arrows(x0 = i, y0 = TVJ_Late_bycomp_boot[i,1],
         x1 = i, y1 = TVJ_Late_bycomp_boot[i,2],code = 3,
         col = 'black',angle = 90, length = .1,lty = 1)
}

title(main="Late Game", line=.5, cex.main=1.8)
title(xlab="Composition Type", line=1.6, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,at = xat,labels = xlabels)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,6.25),y=c(0,0),lwd =2)

# Agreement Test ----------------------------------------------------------

Early_yes <- OnlineData[which(OnlineData$EarlyTVJ==1),]
Early_no <- OnlineData[which(OnlineData$EarlyTVJ==0),]

Late_yes <- OnlineData[which(OnlineData$LateTVJ==1),]
Late_no <- OnlineData[which(OnlineData$LateTVJ==0),]

Agree_Early_yesyes <- Early_yes[which(Early_yes$SpeakerEarly==1),]
#Agree_Early_nono <- Early_no[which(OnlineData$SpeakerEarly==0),]

Agree_Late_yesyes <- Late_yes[which(Late_yes$SpeakerLate==1),]
#Agree_Late_nono <- Late_no[which(OnlineData$SpeakerLate==0),]

#Disagree_Early_yesno <- Early_yes[which(OnlineData$SpeakerEarly==0),]
Disagree_Early_noyes <- Early_no[which(Early_no$SpeakerEarly==1),]

#Disagree_Late_yesno <- Late_yes[which(OnlineData$SpeakerLate==0),]
Disagree_Late_noyes <- Late_no[which(Late_no$SpeakerLate==1),]

Diff_Agree_Early_yesyes <- Agree_Early_yesyes$EarlyResponse - Agree_Early_yesyes$EarlyPrior
#Diff_Agree_Early_nono <- Agree_Early_nono$EarlyResponse - Agree_Early_nono$EarlyPrior

Diff_Agree_Late_yesyes <- Agree_Late_yesyes$LateResponse - Agree_Late_yesyes$LatePrior
#Diff_Agree_Late_nono <- Agree_Late_nono$LateResponse - Agree_Late_nono$LatePrior

#Diff_Disagree_Early_yesno <- Disagree_Early_yesno$EarlyResponse - Disagree_Early_yesno$EarlyPrior
Diff_Disagree_Early_noyes <- Disagree_Early_noyes$EarlyResponse - Disagree_Early_noyes$EarlyPrior

#Diff_Disagree_Late_yesno <- Disagree_Late_yesno$LateResponse - Disagree_Late_yesno$LatePrior
Diff_Disagree_Late_noyes <- Disagree_Late_noyes$LateResponse - Disagree_Late_noyes$LatePrior

Diff_Agree_Early <- c(Diff_Agree_Early_yesyes#,Diff_Agree_Early_nono
                      )
Diff_Agree_Late <- c(Diff_Agree_Late_yesyes#,#Diff_Agree_Late_nono
                     )

Diff_Disagree_Early <- c(Diff_Disagree_Early_noyes#,Diff_Disagree_Early_yesno
                         )
Diff_Disagree_Late <- c(Diff_Disagree_Late_noyes#,Diff_Disagree_Late_yesno
                        )

Diff_Agree <- c(Diff_Agree_Early,Diff_Agree_Late)
Diff_Disagree <- c(Diff_Disagree_Early,Diff_Disagree_Late)

c1 <- rgb(173,216,230,max = 255, alpha = 130, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 100, names = "lt.pink")

b <- min(c(Diff_Agree,Diff_Disagree)) - 0.001 # Set the minimum for the breakpoints
e <- max(c(Diff_Agree,Diff_Disagree)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 12) # Make a neat vector for the breakpoints
hist1 <- hist(Diff_Agree)
hist2 <- hist(Diff_Disagree)
par(mfrow=c(1,1))
par(mar = c(5, 4, 1, 2))
plot(hist1,col = c1,ylim = c(0,200),xlim = c(-1,1),xlab = 'Interpretation - Prior Estimate',main = "Interpretation Shift",
     axes = FALSE)
plot(hist2,col = c2,add=TRUE)
abline(v=0,lty=3)
axis(side = 2, pos= -1)
axis(side = 1, pos=0)
legend(-.9,100,legend = c('Agreement','Disagreement'),pch = c(15,15),
       col = c(c1,c2),cex = 1.3)


# Sanity Check ------------------------------------------------------------
disagree_plot <- matrix(NaN, nrow = 251, ncol = 3)

ticker=1
for(i in 1:nrow(OnlineData)){
  if(OnlineData$SpeakerEarly[i]==1){
    if(OnlineData$EarlyTVJ[i]==0){
      disagree_plot[ticker,1] <- OnlineData$EarlyPrior[i]
      disagree_plot[ticker,2] <- OnlineData$EarlyResponse[i]
      disagree_plot[ticker,3] <- OnlineData$Early_psi_bern[i]
      ticker = ticker+1
    }
  }
  
  if(OnlineData$SpeakerLate[i]==1){
    if(OnlineData$LateTVJ[i]==0){
      disagree_plot[ticker,1] <- OnlineData$LatePrior[i]
      disagree_plot[ticker,2] <- OnlineData$LateResponse[i]
      disagree_plot[ticker,3] <- OnlineData$Late_psi_bern[i]
      ticker = ticker+1
    }
  }
}

plot(disagree_plot[,1],disagree_plot[,2],col=rgb(disagree_plot[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation")

abline(a=0,b=1,lty = 3)
