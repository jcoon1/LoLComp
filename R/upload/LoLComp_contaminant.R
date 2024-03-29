# clears workspace:  
rm(list=ls()) 

setwd('C:/Users/Jeff/Documents/R/LoLComp/upload')

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)

source("BootstrapMean95.R")
source("PlotFunc_horizontal.R")

# Load data --------------------------------------------------------

OnlineData <- read_csv("OnlineData.csv",col_names=TRUE)

nExperts = max(OnlineData$Participant)
nItems = 12

OnlineData$AgreeEarly <- rep(99,nrow(OnlineData))
OnlineData$AgreeLate <- rep(99,nrow(OnlineData))

for(i in 1:nrow(OnlineData)){
  if(OnlineData$SpeakerEarly[i]==1){
    OnlineData$AgreeEarly[i] <- OnlineData$EarlyTVJ[i] == 1
  }
  if(OnlineData$SpeakerLate[i]==1){
    OnlineData$AgreeLate[i] <- OnlineData$LateTVJ[i] == 1
  }
}

SONAData <- read_csv("SONAData.csv",col_names=TRUE)

nNovices = max(SONAData$Participant)

totalNovice = nItems*nNovices

# Prep Data for Early Condition ---------------------------------------------------------------
nBeta = 7
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

nAgree_early <- rep(NaN,nExperts)
nDisagree_early <- rep(NaN,nExperts)

for(i in 1:nExperts){
  temp <- OnlineData[which(OnlineData$Participant==i),]
  nAgree_early[i] <- length(which(temp$AgreeEarly==1))
  nDisagree_early[i] <- length(which(temp$AgreeEarly==0))
}

nExperts_agree_early <- length(which(nAgree_early>0))
nExperts_disagree_early <- length(which(nDisagree_early>0))

nAgree_early <- rep(NaN,nExperts_agree_early)
nDisagree_early <- rep(NaN,nExperts_disagree_early)

ticker = 1
for(i in 1:nExperts){
  temp <- OnlineData[which(OnlineData$Participant==i),]
  if(length(which(temp$AgreeEarly==1))>0){
    nAgree_early[ticker] <- length(which(temp$AgreeEarly==1))
    ticker = ticker+1
  }
}

ticker = 1
for(i in 1:nExperts){
  temp <- OnlineData[which(OnlineData$Participant==i),]
  if(length(which(temp$AgreeEarly==0))>0){
    nDisagree_early[ticker] <- length(which(temp$AgreeEarly==0))
    ticker = ticker+1
  }
}

X_agree_early <- matrix(0, nrow = sum(nAgree_early),ncol = nBeta)
Y_agree_early <- rep(NaN,sum(nAgree_early))
Participant_list_agree_early <- rep(NaN,sum(nAgree_early))

X_disagree_early <- matrix(0, nrow = sum(nDisagree_early),ncol = nBeta)
Y_disagree_early <- rep(NaN,sum(nDisagree_early))
Participant_list_disagree_early <- rep(NaN,sum(nDisagree_early))
Prior_early <- rep(NaN,sum(nDisagree_early))

tickerAgree_early = 1
tickerDisagree_early = 1

for(i in 1:nrow(OnlineData)){
  if(OnlineData$AgreeEarly[i]==1){
    Y_agree_early[tickerAgree_early] <- OnlineData$EarlyResponse[i]
      
    X_agree_early[tickerAgree_early,(OnlineData$CompType[i]+1)] <- 1
    Participant_list_agree_early[tickerAgree_early] <- OnlineData$Participant[i]
      
    tickerAgree_early = tickerAgree_early+1
  }
  if(OnlineData$AgreeEarly[i]==0){
    Y_disagree_early[tickerDisagree_early] <- OnlineData$EarlyResponse[i]
    Prior_early[tickerDisagree_early] <- OnlineData$EarlyPrior[i]
    
    X_disagree_early[tickerDisagree_early,(OnlineData$CompType[i]+1)] <- 1
    Participant_list_disagree_early[tickerDisagree_early] <- OnlineData$Participant[i]
    
    tickerDisagree_early = tickerDisagree_early+1
  }
}

index_agree_early <- matrix(NaN, nrow = nExperts_agree_early, ncol = (nItems/2))
index_disagree_early <- matrix(NaN, nrow = nExperts_disagree_early, ncol = (nItems/2))

for(i in 1:nExperts_agree_early){
  for(j in 1:nAgree_early[i]){
    index_agree_early[i,j] <- j + sum(nAgree_early[1:i]) - nAgree_early[i]
  }
}

for(i in 1:nExperts_agree_early){
  for(j in 1:nDisagree_early[i]){
    index_disagree_early[i,j] <- j + sum(nDisagree_early[1:i]) - nDisagree_early[i]
  }
}

# Run early condition in JAGS -------------------------------------------
nSamples = 4000

DIC <- rep(NaN,2)

X_agree <- X_agree_early
Y_agree <- Y_agree_early
Participant_list_agree <- Participant_list_agree_early
nAgree <- nAgree_early

X_disagree <- X_disagree_early
Y_disagree <- Y_disagree_early
Participant_list_disagree <- Participant_list_disagree_early
Prior <- Prior_early
nDisagree <- nDisagree_early

bPrec <- bPrec
nBeta <- nBeta
nExperts <- nExperts
nExperts_agree <- nExperts_agree_early
nExperts_disagree <- nExperts_disagree_early
index_agree <- index_agree_early
index_disagree <- index_disagree_early

data <- list("X_agree","Y_agree","Participant_list_agree","nAgree", 
             "X_disagree","Y_disagree","Participant_list_disagree","Prior","nDisagree",
             "bPrec","nBeta","index_agree","index_disagree",
             "nExperts","nExperts_agree","nExperts_disagree") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_exp","sigma_con","sigma_nov","alpha",
                "phi_disagree","psi_disagree","psi_disagree_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_contaminant.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

Beta_Early_samples <- samples$BUGSoutput$sims.list$BETA

SigmaExp_Early_samples <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Early_samples <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Early_samples <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Early_samples <- samples$BUGSoutput$sims.list$tau

summary_Early <- samples$BUGSoutput$summary

DIC[1] <- samples$BUGSoutput$DIC

test_vector_agree <- rep(NaN, sum(nAgree_early))

ticker = 1
for(i in 1:nExperts_agree){
  for(j in 1:nAgree[i]){
    test_vector_agree[ticker] <- index_agree[i,j]
    ticker = ticker+1
  }
}

test_vector_disagree <- rep(NaN, sum(nDisagree_early))

ticker = 1
for(i in 1:nExperts_disagree){
  for(j in 1:nDisagree[i]){
    test_vector_disagree[ticker] <- index_disagree[i,j]
    ticker = ticker+1
  }
}

# Prep Data for Late Condition ---------------------------------------------------------------
nBeta = 7
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

nAgree_late <- rep(NaN,nExperts)
nDisagree_late <- rep(NaN,nExperts)

for(i in 1:nExperts){
  temp <- OnlineData[which(OnlineData$Participant==i),]
  nAgree_late[i] <- length(which(temp$AgreeLate==1))
  nDisagree_late[i] <- length(which(temp$AgreeLate==0))
}

nExperts_agree_late <- length(which(nAgree_late>0))
nExperts_disagree_late <- length(which(nDisagree_late>0))

nAgree_late <- rep(NaN,nExperts_agree_late)
nDisagree_late <- rep(NaN,nExperts_disagree_late)

ticker = 1
for(i in 1:nExperts){
  temp <- OnlineData[which(OnlineData$Participant==i),]
  if(length(which(temp$AgreeLate==1))>0){
    nAgree_late[ticker] <- length(which(temp$AgreeLate==1))
    ticker = ticker+1
  }
}

ticker = 1
for(i in 1:nExperts){
  temp <- OnlineData[which(OnlineData$Participant==i),]
  if(length(which(temp$AgreeLate==0))>0){
    nDisagree_late[ticker] <- length(which(temp$AgreeLate==0))
    ticker = ticker+1
  }
}

X_agree_late <- matrix(0, nrow = sum(nAgree_late),ncol = nBeta)
Y_agree_late <- rep(NaN,sum(nAgree_late))
Participant_list_agree_late <- rep(NaN,sum(nAgree_late))

X_disagree_late <- matrix(0, nrow = sum(nDisagree_late),ncol = nBeta)
Y_disagree_late <- rep(NaN,sum(nDisagree_late))
Participant_list_disagree_late <- rep(NaN,sum(nDisagree_late))
Prior_late <- rep(NaN,sum(nDisagree_late))

tickerAgree_late = 1
tickerDisagree_late = 1

for(i in 1:nrow(OnlineData)){
  if(OnlineData$AgreeLate[i]==1){
    Y_agree_late[tickerAgree_late] <- OnlineData$LateResponse[i]
    
    X_agree_late[tickerAgree_late,(OnlineData$CompType[i]+1)] <- 1
    Participant_list_agree_late[tickerAgree_late] <- OnlineData$Participant[i]
    
    tickerAgree_late = tickerAgree_late+1
  }
  if(OnlineData$AgreeLate[i]==0){
    Y_disagree_late[tickerDisagree_late] <- OnlineData$LateResponse[i]
    Prior_late[tickerDisagree_late] <- OnlineData$LatePrior[i]
    
    X_disagree_late[tickerDisagree_late,(OnlineData$CompType[i]+1)] <- 1
    Participant_list_disagree_late[tickerDisagree_late] <- OnlineData$Participant[i]
    
    tickerDisagree_late = tickerDisagree_late+1
  }
}

index_agree_late <- matrix(NaN, nrow = nExperts_agree_late, ncol = (nItems/2))
index_disagree_late <- matrix(NaN, nrow = nExperts_disagree_late, ncol = (nItems/2))

for(i in 1:nExperts_agree_late){
  for(j in 1:nAgree_late[i]){
    index_agree_late[i,j] <- j + sum(nAgree_late[1:i]) - nAgree_late[i]
  }
}

for(i in 1:nExperts_agree_late){
  for(j in 1:nDisagree_late[i]){
    index_disagree_late[i,j] <- j + sum(nDisagree_late[1:i]) - nDisagree_late[i]
  }
}

# Run late condition in JAGS ----------------------------------------------

nSamples = 4000

DIC <- rep(NaN,2)

X_agree <- X_agree_late
Y_agree <- Y_agree_late
Participant_list_agree <- Participant_list_agree_late
nAgree <- nAgree_late

X_disagree <- X_disagree_late
Y_disagree <- Y_disagree_late
Participant_list_disagree <- Participant_list_disagree_late
Prior <- Prior_late
nDisagree <- nDisagree_late

bPrec <- bPrec
nBeta <- nBeta
nExperts <- nExperts
nExperts_agree <- nExperts_agree_late
nExperts_disagree <- nExperts_disagree_late
index_agree <- index_agree_late
index_disagree <- index_disagree_late

data <- list("X_agree","Y_agree","Participant_list_agree","nAgree", 
             "X_disagree","Y_disagree","Participant_list_disagree","Prior","nDisagree",
             "bPrec","nBeta","index_agree","index_disagree",
             "nExperts","nExperts_agree","nExperts_disagree") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_exp","sigma_con","sigma_nov","alpha",
                "phi_disagree","psi_disagree","psi_disagree_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_contaminant.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

Beta_Late_samples <- samples$BUGSoutput$sims.list$BETA

SigmaExp_Late_samples <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Late_samples <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Late_samples <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Late_samples <- samples$BUGSoutput$sims.list$tau

summary_Late <- samples$BUGSoutput$summary

DIC[1] <- samples$BUGSoutput$DIC

test_vector_agree <- rep(NaN, sum(nAgree_late))

ticker = 1
for(i in 1:nExperts_agree){
  for(j in 1:nAgree[i]){
    test_vector_agree[ticker] <- index_agree[i,j]
    ticker = ticker+1
  }
}

test_vector_disagree <- rep(NaN, sum(nDisagree_late))

ticker = 1
for(i in 1:nExperts_disagree){
  for(j in 1:nDisagree[i]){
    test_vector_disagree[ticker] <- index_disagree[i,j]
    ticker = ticker+1
  }
}

# Compare expert and novice sigmas ----------------------------------------

SigmaTest_yesEarly <- matrix(NaN,nrow = length(SigmaExp_Late_samples),ncol = 2)
SigmaTest_yesLate <- matrix(NaN,nrow = length(SigmaExp_Late_samples),ncol = 2)

SigmaDiff_yesEarly <- rep(NaN,length(SigmaExp_Late_samples))
SigmaDiff_yesLate <- rep(NaN,length(SigmaExp_Late_samples))

for(i in 1:length(SigmaExp_Late_samples)){
  SigmaTest_yesEarly[i,1] <- Tau_Early_samples[i]^2/
    (Tau_Early_samples[i]^2 + SigmaExp_Early_samples[i]^2)
  SigmaTest_yesEarly[i,2] <- Tau_Early_samples[i]^2/
    (Tau_Early_samples[i]^2 + SigmaNov_Early_samples[i]^2)
  
  SigmaDiff_yesEarly[i] <- SigmaExp_Early_samples[i]-SigmaNov_Early_samples[i]
  SigmaDiff_yesLate[i] <- SigmaExp_Late_samples[i]-SigmaNov_Late_samples[i]
  
  SigmaTest_yesLate[i,1] <- Tau_Late_samples[i]^2/
    (Tau_Late_samples[i]^2 + SigmaExp_Late_samples[i]^2)
  SigmaTest_yesLate[i,2] <- Tau_Late_samples[i]^2/
    (Tau_Late_samples[i]^2 + SigmaNov_Late_samples[i]^2)
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

sortSigmaTest_yesEarly <- matrix(NaN,nrow = length(SigmaExp_Late_samples),ncol = 2)
sortSigmaTest_yesLate <- matrix(NaN,nrow = length(SigmaExp_Late_samples),ncol = 2)

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
