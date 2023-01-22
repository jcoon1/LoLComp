# clears workspace:  
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)

source("BootstrapMean95.R")
source("PlotFunc_horizontal.R")

# Load data --------------------------------------------------------

## Each row of data is a participant responding to a new composition. They each saw
## 12 compositions and interpreted two generalizations (early & late) for each
## composition. The hypothetical speaker said they would endorse half of each 
## but wouldn't endorse the other half. We're only interested in the statements
## that speakers WOULD endorse (SpeakerEarly or SpeakerLate = 1). Thus there are
## 6 statements of interest for each type of statement (early & late), evenly
## distrubited so that there's one for each typ of composition.

## Comp is the specific composition because there are 2 exemplars
## for each type of composition. CompType goes 1:6, meaning (E+,E0,E-,L+,LO,L-).
## CompBorderline collapses this further into 1:3 meaning (+,0,-).

## Response is participants' interpretation, TVJ is whether they would endorse
## the generalization themselves, and Prior is their own expectations for how
## often the generalization would apply (Early & Late for each)
OnlineData <- read_csv("OnlineData.csv",col_names=TRUE)

## How many LoL people participated
nExperts = max(OnlineData$Participant)
## Total number of trials of interest for each participant
nItems = 12

## Identify whether this is an example of agreement or disagreement. 0 = disagree,
## 1 = agree, 99 = not a trial of interest

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

## Same structure for naive data; no need to worry about agreement/disagreement
SONAData <- read_csv("SONAData.csv",col_names=TRUE)
nNovices = max(SONAData$Participant)

# Prep Data for Early Condition ---------------------------------------------------------------
## Treating generlizations about early game and late game as separate statements
## (and thus essentially replications of each other)

## First beta is naive participant interpretation, the next 6 are experienced 
## by condition
nBeta = 7
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}


## Set up a count of how many times each expert participant agreed and disagreed
## with the speaker
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

## preallocate stuff
X_agree_early <- matrix(0, nrow = sum(nAgree_early),ncol = nBeta)
Y_agree_early <- rep(NaN,sum(nAgree_early))
Participant_list_agree_early <- rep(NaN,sum(nAgree_early))

X_disagree_early <- matrix(0, nrow = sum(nDisagree_early),ncol = nBeta)
Y_disagree_early <- rep(NaN,sum(nDisagree_early))
Participant_list_disagree_early <- rep(NaN,sum(nDisagree_early))
Prior_early <- rep(NaN,sum(nDisagree_early))

tickerAgree_early = 1
tickerDisagree_early = 1

## Make data, containing only examples of expert agreement. Y is their interpretation.
## "Prior" is the X variable (input) for the latent mixture part of the model and is thus
## only needed for disagreements, since that's the only part of the model with the mixture.
## X is the indicator matrix for turning the betas on and off as appropriate in the main
## model. Since these are only ever expert participants, the first column is always 0.
## For each row, there will always be a single column with a 1 corresponding to the condition.

## Participant_list is a vector identifying which participant is responsible for
## each data point.
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

## Fancy indexing matrix that keeps track of trial number for looping purposes.
## The agree version just needs to move through Participant_list, X, and Y. The
## disagree version needs to also move through Prior. In both cases, row is
## participant, and column is how many agrees/disagrees that participant has
## done so far. There are only 48 rows because there is a participant who
## always disagreed and a participant who always agreed. It should still be
## fine because even though the model only cycles through 48 participants,
## it identifies participant number by using this index matrix to find where
## it should be in the Participant_list vectors, which skip a participant ID
## where appropriate.

index_agree_early <- matrix(NaN, nrow = nExperts_agree_early, ncol = (nItems/2))
index_disagree_early <- matrix(NaN, nrow = nExperts_disagree_early, ncol = (nItems/2))

for(i in 1:nExperts_agree_early){
  for(j in 1:nAgree_early[i]){
    index_agree_early[i,j] <- j + sum(nAgree_early[1:i]) - nAgree_early[i]
  }
}

for(i in 1:nExperts_disagree_early){
  for(j in 1:nDisagree_early[i]){
    index_disagree_early[i,j] <- j + sum(nDisagree_early[1:i]) - nDisagree_early[i]
  }
}

Novice_early <- length(which(SONAData$SpeakerEarly==1))

## X and Y for novices. X is always just a 1 in the first column because these
## are all novice participants

X_nov_early <- matrix(0, nrow = Novice_early, ncol = nBeta)
X_nov_early[,1] <- rep(1,Novice_early)
ticker = 1
Y_nov_early <- rep(NaN, Novice_early)

for(i in 1:nrow(SONAData)){
  if(SONAData$SpeakerEarly[i]==1){
    Y_nov_early[ticker] <- SONAData$EarlyResponse[i]
    ticker = ticker + 1
  }
}

# Run early condition in JAGS -------------------------------------------
nSamples = 4000

DIC <- rep(NaN,2)

X_agree <- X_agree_early
Y_agree <- Y_agree_early
## Which participant produced each data point
Participant_list_agree <- Participant_list_agree_early
## How many generalizations did each participant agree with
nAgree <- nAgree_early

X_disagree <- X_disagree_early
Y_disagree <- Y_disagree_early
## Which participant produced each data point
Participant_list_disagree <- Participant_list_disagree_early
## Prior estimates for each disagreement
Prior <- Prior_early
## How many generalizations did each participant disagree with
nDisagree <- nDisagree_early

X_nov <- X_nov_early
Y_nov <- Y_nov_early
nNovices <- nNovices
nTrials <- length(which(SONAData$SpeakerEarly[1:12]==1))

bPrec <- bPrec
nBeta <- nBeta
## Total experts
nExperts <- nExperts
## Number of experts who agreed and disagreed with any generalizations
nExperts_agree <- nExperts_agree_early
nExperts_disagree <- nExperts_disagree_early
index_agree <- index_agree_early
index_disagree <- index_disagree_early

data <- list("X_agree","Y_agree","Participant_list_agree","nAgree", 
             "X_disagree","Y_disagree","Participant_list_disagree","Prior","nDisagree",
             "bPrec","nBeta","index_agree","index_disagree",
             "nExperts","nExperts_agree","nExperts_disagree",
             "nNovices","X_nov","Y_nov","nTrials") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_exp","sigma_con","sigma_nov","alpha",
                "phi_disagree",
                "psi_disagree","psi_disagree_bern","sigma_psi")

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

alpha_Early_samples <- samples$BUGSoutput$sims.list$alpha

phi_Early_samples <- samples$BUGSoutput$sims.list$psi_disagree

psi_Early_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

Phi_Early_samples <- samples$BUGSoutput$sims.list$phi_disagree

SigmaExp_Early_samples <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Early_samples <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Early_samples <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Early_samples <- samples$BUGSoutput$sims.list$tau

summary_Early <- samples$BUGSoutput$summary

DIC[1] <- samples$BUGSoutput$DIC


## Test the fancy indexing matrix to make sure it's built properly
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
##  Repeat everything for the generalizations about the late game
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

for(i in 1:nExperts_disagree_late){
  for(j in 1:nDisagree_late[i]){
    index_disagree_late[i,j] <- j + sum(nDisagree_late[1:i]) - nDisagree_late[i]
  }
}

Novice_late <- length(which(SONAData$SpeakerLate==1))

X_nov_late <- matrix(0, nrow = Novice_late, ncol = nBeta)
X_nov_late[,1] <- rep(1,Novice_late)

ticker = 1
Y_nov_late <- rep(NaN, Novice_late)

for(i in 1:nrow(SONAData)){
  if(SONAData$SpeakerLate[i]==1){
    Y_nov_late[ticker] <- SONAData$LateResponse[i]
    ticker = ticker + 1
    
  }
}

# Run late condition in JAGS ----------------------------------------------

nSamples = 4000

X_agree <- X_agree_late
Y_agree <- Y_agree_late
Participant_list_agree <- Participant_list_agree_late
nAgree <- nAgree_late

X_disagree <- X_disagree_late
Y_disagree <- Y_disagree_late
Participant_list_disagree <- Participant_list_disagree_late
Prior <- Prior_late
nDisagree <- nDisagree_late

X_nov <- X_nov_late
Y_nov <- Y_nov_late
nNovices <- nNovices
nTrials <- length(which(SONAData$SpeakerLate[1:12]==1))

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
             "nExperts","nExperts_agree","nExperts_disagree",
             "nNovices","X_nov","Y_nov","nTrials") # to be passed on to JAGS
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

alpha_Late_samples <- samples$BUGSoutput$sims.list$alpha

phi_Late_samples <- samples$BUGSoutput$sims.list$psi_disagree

Phi_Late_samples <- samples$BUGSoutput$sims.list$phi_disagree

psi_Late_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

SigmaExp_Late_samples <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Late_samples <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Late_samples <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Late_samples <- samples$BUGSoutput$sims.list$tau

summary_Late <- samples$BUGSoutput$summary

DIC[2] <- samples$BUGSoutput$DIC

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


# Plot beta inferences from raw -------------------------------------------
plotValues_early <- matrix(NaN, nrow = 7,ncol = 3)
plotValues_late <- matrix(NaN, nrow = 7,ncol = 3)

loconf = nrow(Beta_Early_samples) * .025
hiconf = nrow(Beta_Early_samples) * .975

tempsort <- sort(Beta_Early_samples[,1])
plotValues_early[1,1] <- mean(tempsort[1500:1001])
plotValues_early[1,2] <- tempsort[loconf]
plotValues_early[1,3] <- tempsort[hiconf]

tempsort <- sort(Beta_Late_samples[,1])
plotValues_late[1,1] <- mean(tempsort[1500:1001])
plotValues_late[1,2] <- tempsort[loconf]
plotValues_late[1,3] <- tempsort[hiconf]

for(i in 2:7){
  
  tempsort <- sort(Beta_Early_samples[,i])
  plotValues_early[i,1] <- mean(tempsort[1500:1001])
  plotValues_early[i,2] <- tempsort[loconf]
  plotValues_early[i,3] <- tempsort[hiconf]
  
  
  tempsort <- sort(Beta_Late_samples[,i])
  plotValues_late[i,1] <- mean(tempsort[1500:1001])
  plotValues_late[i,2] <- tempsort[loconf]
  plotValues_late[i,3] <- tempsort[hiconf]
}

plot(1:6,plotValues_early[2:7,1],ylim = c(.15,.9),main = "Yes Early",xlab = 'beta',ylab = 'Interpretation')
abline(h = plotValues_early[1,1],col = '#DC143C')
abline(h= plotValues_early[1,2],col = '#DC143C',lty = 2)
abline(h= plotValues_early[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = plotValues_early[i+1,2],
         x1 = i, y1 = plotValues_early[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}
legend(3.5,.45,legend = c('Novice Interpretation','Expert Interpretation'),pch = c(NA,1),lty = c(1,NA),
       col = c('#DC143C','black'),cex = .8)

plot(1:6,plotValues_late[2:7,1],ylim = c(.15,.9),main = "Yes Late",xlab = 'beta',ylab = 'Interpretation')
abline(h = plotValues_late[1,1],col = '#DC143C')
abline(h= plotValues_late[1,2],col = '#DC143C',lty = 2)
abline(h= plotValues_late[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = plotValues_late[i+1,2],
         x1 = i, y1 = plotValues_late[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}


nNovResponses = nNovices*nTrials
SONA_empirical_Early <- matrix(NaN,nrow = nNovResponses, ncol = 2)
SONA_empirical_Late <- matrix(NaN,nrow = nNovResponses, ncol = 2)

SONA_empirical_Early[,1] <- SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)]
SONA_empirical_Early[,2] <- SONAData$CompType[which(SONAData$SpeakerEarly==1)]

SONA_empirical_Late[,1] <- SONAData$LateResponse[which(SONAData$SpeakerLate==1)]
SONA_empirical_Late[,2] <- SONAData$CompType[which(SONAData$SpeakerLate==1)]

SONA_empirical_Early_mean <- mean(SONA_empirical_Early[,1])
SONA_empirical_Late_mean <- mean(SONA_empirical_Late[,1])

SONA_empirical_bycomp_Early <- rep(NaN,6)
SONA_empirical_bycomp_Late <- rep(NaN,6)

for(i in 1:max(SONAData$CompType)){
  SONA_empirical_bycomp_Early[i] <- mean(SONA_empirical_Early[which(SONA_empirical_Early[,2]==i),1])
  SONA_empirical_bycomp_Late[i] <- mean(SONA_empirical_Late[which(SONA_empirical_Late[,2]==i),1])
}

nExpResponses = nExperts*nTrials
Online_empirical_Early <- matrix(NaN,nrow = nExpResponses, ncol = 2)
Online_empirical_Late <- matrix(NaN,nrow = nExpResponses, ncol = 2)

Online_empirical_Early[,1] <- OnlineData$EarlyResponse[which(OnlineData$SpeakerEarly==1)]
Online_empirical_Early[,2] <- OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]

Online_empirical_Late[,1] <- OnlineData$LateResponse[which(OnlineData$SpeakerLate==1)]
Online_empirical_Late[,2] <- OnlineData$CompType[which(OnlineData$SpeakerLate==1)]

Online_empirical_Early_mean <- mean(Online_empirical_Early[,1])
Online_empirical_Late_mean <- mean(Online_empirical_Late[,1])

Online_empirical_bycomp_Early <- rep(NaN,6)
Online_empirical_bycomp_Late <- rep(NaN,6)

for(i in 1:max(OnlineData$CompType)){
  Online_empirical_bycomp_Early[i] <- mean(Online_empirical_Early[which(Online_empirical_Early[,2]==i),1])
  Online_empirical_bycomp_Late[i] <- mean(Online_empirical_Late[which(Online_empirical_Late[,2]==i),1])
}


par(mfrow=c(1,1))
par(mar = c(3.5, 3, 1, 2))

plot1 = MainPlot(plotValues_early[1,1:3],SONA_empirical_Early_mean,SONA_empirical_bycomp_Early,
                 plotValues_early[2:7,1:3],Online_empirical_bycomp_Early,'',
                 if_legend = 1,leftmost = 1)
par(mar = c(3.5, 3, 1, 2))
plot2 = MainPlot(plotValues_late[1,1:3],SONA_empirical_Late_mean,SONA_empirical_bycomp_Late,
                 plotValues_late[2:7,1:3],Online_empirical_bycomp_Late,'',
                 if_legend = 0, leftmost = 0)

# Make a fancy plot where color indicates psi_bern ---------------------------
par(mar = c(4.5, 4, 1, 1))

X_prior_graph_latentMix_psi <- cbind(Prior_early,
                                     summary_Early[140:(139+sum(nDisagree_early)),1])

disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi),ncol = 5)

disagree_plot[,1] = X_prior_graph_latentMix_psi[,1]
disagree_plot[,2] = Y_disagree_early
disagree_plot[,3] = X_prior_graph_latentMix_psi[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi)){
  disagree_plot[i,4] = which(X_disagree_early[i,]==1)-1
  disagree_plot[i,5] = which(X_disagree_early[i,]==1)-1
  if(disagree_plot[i,4]>3){
    disagree_plot[i,5] <- 7-disagree_plot[i,4]
  }
}


plot(disagree_plot[,1],disagree_plot[,2],col=rgb(disagree_plot[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation", main= "early 7b")

abline(a=0,b=1,lty = 3)

text(disagree_plot[,1],disagree_plot[,2],label = disagree_plot[,4],col=rgb(disagree_plot[,3],0,0))




X_prior_graph_latentMix_psi_late <- cbind(Prior_late,
                                     summary_Late[141:(140+sum(nDisagree_late)),1])

disagree_plot_late = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi_late),ncol = 5)

disagree_plot_late[,1] = X_prior_graph_latentMix_psi_late[,1]
disagree_plot_late[,2] = Y_disagree_late
disagree_plot_late[,3] = X_prior_graph_latentMix_psi_late[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi_late)){
  disagree_plot_late[i,4] = which(X_disagree_late[i,]==1)-1
  disagree_plot_late[i,5] = which(X_disagree_late[i,]==1)-1
  if(disagree_plot_late[i,4]>3){
    disagree_plot_late[i,5] <- 7-disagree_plot_late[i,4]
  }
}


plot(disagree_plot_late[,1],disagree_plot_late[,2],col=rgb(disagree_plot_late[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation",type = 'n', main = "late 7b")

text(disagree_plot_late[,1],disagree_plot_late[,2],label = disagree_plot_late[,4],col=rgb(disagree_plot_late[,3],0,0))

abline(a=0,b=1,lty = 3)


# By participant (watch out) -----------------------------------------------
par(mfrow=c(6,4))

for(i in 1:nExperts){
  if(length(which(Participant_list_disagree_early==i))>1){
    temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
    plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
    abline(a=0,b=1,lty = 3)
  }
  if(length(which(Participant_list_disagree_early==i))==1){
    temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
    plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
    abline(a=0,b=1,lty = 3)
  }
}

for(i in 1:nExperts){
  if(length(which(Participant_list_disagree_late==i))>1){
    temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
    plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
    abline(a=0,b=1,lty = 3)
  }
  if(length(which(Participant_list_disagree_late==i))==1){
    temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
    plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
    abline(a=0,b=1,lty = 3)
  }
}


# Run Early and Late Together (shared participant-level info) -------------

nSamples = 4000

nBeta <- 14

X_agree_early14 <- matrix(0,nrow = sum(nAgree_early), ncol = nBeta)
X_agree_late14 <- matrix(0,nrow = sum(nAgree_late), ncol = nBeta)

X_agree_early14[1:nrow(X_agree_early14),1:7] <- X_agree_early
X_agree_late14[1:nrow(X_agree_late14),8:14] <- X_agree_late
Y_agree_early <- Y_agree_early
Y_agree_late <- Y_agree_late
Participant_list_agree_early <- Participant_list_agree_early
Participant_list_agree_late <- Participant_list_agree_late
nAgree_early <- nAgree_early
nAgree_late <- nAgree_late

X_disagree_early14 <- matrix(0,nrow = nrow(X_disagree_early), ncol = nBeta)
X_disagree_late14 <- matrix(0,nrow = nrow(X_disagree_late), ncol = nBeta)

X_disagree_early14[1:nrow(X_disagree_early),1:7] <- X_disagree_early
X_disagree_late14[1:nrow(X_disagree_late),8:14] <- X_disagree_late
Y_disagree_early <- Y_disagree_early
Y_disagree_late <- Y_disagree_late
Participant_list_disagree_early <- Participant_list_disagree_early
Participant_list_disagree_late <- Participant_list_disagree_late
nDisagree_early <- nDisagree_early
nDisagree_late <- nDisagree_late
Prior_early <- Prior_early
Prior_late <- Prior_late

X_nov_early14 <- matrix(0,nrow = nrow(X_nov_early), ncol = nBeta)
X_nov_early14[1:nrow(X_nov_early),1:7] <- X_nov_early

X_nov_late14 <- matrix(0,nrow = nrow(X_nov_late), ncol = nBeta)
X_nov_late14[1:nrow(X_nov_late),8:14] <- X_nov_late
nNovices <- nNovices
nTrials <- length(which(SONAData$SpeakerLate[1:12]==1))

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

bPrec <- bPrec
nBeta <- nBeta
nExperts <- nExperts
nExperts_agree_early <- nExperts_agree_early
nExperts_disagree_early <- nExperts_disagree_early
nExperts_agree_late <- nExperts_agree_late
nExperts_disagree_late <- nExperts_disagree_late
index_agree_early <- index_agree_early
index_disagree_early <- index_disagree_early
index_agree_late <- index_agree_late
index_disagree_late <- index_disagree_late

data <- list("X_agree_early14","X_agree_late14","Y_agree_early","Y_agree_late",
             "Participant_list_agree_early", "Participant_list_agree_late",
             "X_disagree_early14","X_disagree_late14","Y_disagree_early","Y_disagree_late",
             "Participant_list_disagree_early","Participant_list_disagree_late",
             "Prior_early","Prior_late","bPrec","nBeta",
             "nAgree_early","nAgree_late","nDisagree_early","nDisagree_late",
             "index_agree_early","index_disagree_early","index_agree_late","index_disagree_late",
             "nExperts","nExperts_agree_early","nExperts_disagree_early","nExperts_agree_late","nExperts_disagree_late",
             "nNovices","X_nov_early14","X_nov_late14","Y_nov_early","Y_nov_late","nTrials") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_exp_early","sigma_con_early","sigma_nov_early",
                "sigma_exp_late","sigma_con_late","sigma_nov_late",
                "phi_disagree","psi_disagree","psi_disagree_bern_early","psi_disagree_bern_late",
                "alpha","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_contaminant_14.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

Beta_All_samples <- samples$BUGSoutput$sims.list$BETA

alpha_All_samples <- samples$BUGSoutput$sims.list$alpha

Phi_All_samples <- samples$BUGSoutput$sims.list$phi_disagree

phi_All_samples <- samples$BUGSoutput$sims.list$psi_disagree

psi_All_samples_early <- samples$BUGSoutput$sims.list$psi_disagree_bern_early
psi_All_samples_late <- samples$BUGSoutput$sims.list$psi_disagree_bern_late


SigmaExp_All_early_samples <- samples$BUGSoutput$sims.list$sigma_exp_early
SigmaCon_All_early_samples <- samples$BUGSoutput$sims.list$sigma_con_early
SigmaNov_All_early_samples <- samples$BUGSoutput$sims.list$sigma_nov_early

SigmaExp_All_late_samples <- samples$BUGSoutput$sims.list$sigma_exp_late
SigmaCon_All_late_samples <- samples$BUGSoutput$sims.list$sigma_con_late
SigmaNov_All_late_samples <- samples$BUGSoutput$sims.list$sigma_nov_late

Tau_All_samples <- samples$BUGSoutput$sims.list$tau

summary_All <- samples$BUGSoutput$summary


# Main plot ---------------------------------------------------------------

plotValues_early <- matrix(NaN, nrow = 7,ncol = 3)
plotValues_late <- matrix(NaN, nrow = 7,ncol = 3)

loconf = nrow(Beta_All_samples) * .025
hiconf = nrow(Beta_All_samples) * .975

tempsort <- sort(Beta_All_samples[,1])
plotValues_early[1,1] <- mean(tempsort[1500:1001])
plotValues_early[1,2] <- tempsort[loconf]
plotValues_early[1,3] <- tempsort[hiconf]

tempsort <- sort(Beta_All_samples[,8])
plotValues_late[1,1] <- mean(tempsort[1500:1001])
plotValues_late[1,2] <- tempsort[loconf]
plotValues_late[1,3] <- tempsort[hiconf]

for(i in 2:7){
  
  tempsort <- sort(Beta_All_samples[,i])
  plotValues_early[i,1] <- mean(tempsort[1500:1001])
  plotValues_early[i,2] <- tempsort[loconf]
  plotValues_early[i,3] <- tempsort[hiconf]
  
  late_ticker = i + 7
  tempsort <- sort(Beta_All_samples[,late_ticker])
  plotValues_late[i,1] <- mean(tempsort[1500:1001])
  plotValues_late[i,2] <- tempsort[loconf]
  plotValues_late[i,3] <- tempsort[hiconf]
}

plot(1:6,plotValues_early[2:7,1],ylim = c(.15,.9),main = "Yes Early",xlab = 'beta',ylab = 'Interpretation')
abline(h = plotValues_early[1,1],col = '#DC143C')
abline(h= plotValues_early[1,2],col = '#DC143C',lty = 2)
abline(h= plotValues_early[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = plotValues_early[i+1,2],
         x1 = i, y1 = plotValues_early[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}
legend(3.5,.45,legend = c('Novice Interpretation','Expert Interpretation'),pch = c(NA,1),lty = c(1,NA),
       col = c('#DC143C','black'),cex = .8)

plot(1:6,plotValues_late[2:7,1],ylim = c(.15,.9),main = "Yes Late",xlab = 'beta',ylab = 'Interpretation')
abline(h = plotValues_late[1,1],col = '#DC143C')
abline(h= plotValues_late[1,2],col = '#DC143C',lty = 2)
abline(h= plotValues_late[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = plotValues_late[i+1,2],
         x1 = i, y1 = plotValues_late[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}


par(mfrow=c(1,1))
par(mar = c(3.5, 3, 1, 2))

plot1 = MainPlot(plotValues_early[1,1:3],SONA_empirical_Early_mean,SONA_empirical_bycomp_Early,
                 plotValues_early[2:7,1:3],Online_empirical_bycomp_Early,'',
                 if_legend = 1,leftmost = 1)
par(mar = c(3.5, 3, 1, 2))
plot2 = MainPlot(plotValues_late[1,1:3],SONA_empirical_Late_mean,SONA_empirical_bycomp_Late,
                 plotValues_late[2:7,1:3],Online_empirical_bycomp_Late,'',
                 if_legend = 0, leftmost = 0)

# Make a fancy plot where color indicates psi_bern ---------------------------
par(mar = c(4.5, 4, 1, 1))

X_prior_graph_latentMix_psi <- cbind(Prior_early,
                                     summary_All[148:(147+sum(nDisagree_early)),1])

disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi),ncol = 5)

disagree_plot[,1] = X_prior_graph_latentMix_psi[,1]
disagree_plot[,2] = Y_disagree_early
disagree_plot[,3] = X_prior_graph_latentMix_psi[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi)){
  disagree_plot[i,4] = which(X_disagree_early[i,]==1)-1
  disagree_plot[i,5] = which(X_disagree_early[i,]==1)-1
  if(disagree_plot[i,4]>3){
    disagree_plot[i,5] <- 7-disagree_plot[i,4]
  }
}


plot(disagree_plot[,1],disagree_plot[,2],col=rgb(disagree_plot[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation", main= "early 14b")

abline(a=0,b=1,lty = 3)

text(disagree_plot[,1],disagree_plot[,2],label = disagree_plot[,4],col=rgb(disagree_plot[,3],0,0))




X_prior_graph_latentMix_psi_late <- cbind(Prior_late,
                                          summary_All[272:(271+sum(nDisagree_late)),1])

disagree_plot_late = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi_late),ncol = 5)

disagree_plot_late[,1] = X_prior_graph_latentMix_psi_late[,1]
disagree_plot_late[,2] = Y_disagree_late
disagree_plot_late[,3] = X_prior_graph_latentMix_psi_late[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi_late)){
  disagree_plot_late[i,4] = which(X_disagree_late[i,]==1)-1
  disagree_plot_late[i,5] = which(X_disagree_late[i,]==1)-1
  if(disagree_plot_late[i,4]>3){
    disagree_plot_late[i,5] <- 7-disagree_plot_late[i,4]
  }
}


plot(disagree_plot_late[,1],disagree_plot_late[,2],col=rgb(disagree_plot_late[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation",type = 'n', main = "late 14b")

text(disagree_plot_late[,1],disagree_plot_late[,2],label = disagree_plot_late[,4],col=rgb(disagree_plot_late[,3],0,0))

abline(a=0,b=1,lty = 3)

# By participant (watch out) -----------------------------------------------
par(mfrow=c(6,4))

for(i in 1:nExperts){
  if(length(which(Participant_list_disagree_early==i))>1){
    temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
    plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
    abline(a=0,b=1,lty = 3)
  }
  if(length(which(Participant_list_disagree_early==i))==1){
    temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
    plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
    abline(a=0,b=1,lty = 3)
  }
}

for(i in 1:nExperts){
  if(length(which(Participant_list_disagree_late==i))>1){
    temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
    plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
    abline(a=0,b=1,lty = 3)
  }
  if(length(which(Participant_list_disagree_late==i))==1){
    temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
    plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
    abline(a=0,b=1,lty = 3)
  }
}

# Run Early and Late Together (correlated Alphas) -------------

nSamples = 4000

nBeta <- 14

X_agree_early14 <- matrix(0,nrow = sum(nAgree_early), ncol = nBeta)
X_agree_late14 <- matrix(0,nrow = sum(nAgree_late), ncol = nBeta)

X_agree_early14[1:nrow(X_agree_early14),1:7] <- X_agree_early
X_agree_late14[1:nrow(X_agree_late14),8:14] <- X_agree_late
Y_agree_early <- Y_agree_early
Y_agree_late <- Y_agree_late
Participant_list_agree_early <- Participant_list_agree_early
Participant_list_agree_late <- Participant_list_agree_late
nAgree_early <- nAgree_early
nAgree_late <- nAgree_late

X_disagree_early14 <- matrix(0,nrow = nrow(X_disagree_early), ncol = nBeta)
X_disagree_late14 <- matrix(0,nrow = nrow(X_disagree_late), ncol = nBeta)

X_disagree_early14[1:nrow(X_disagree_early),1:7] <- X_disagree_early
X_disagree_late14[1:nrow(X_disagree_late),8:14] <- X_disagree_late
Y_disagree_early <- Y_disagree_early
Y_disagree_late <- Y_disagree_late
Participant_list_disagree_early <- Participant_list_disagree_early
Participant_list_disagree_late <- Participant_list_disagree_late
nDisagree_early <- nDisagree_early
nDisagree_late <- nDisagree_late
Prior_early <- Prior_early
Prior_late <- Prior_late

X_nov_early14 <- matrix(0,nrow = nrow(X_nov_early), ncol = nBeta)
X_nov_early14[1:nrow(X_nov_early),1:7] <- X_nov_early

X_nov_late14 <- matrix(0,nrow = nrow(X_nov_late), ncol = nBeta)
X_nov_late14[1:nrow(X_nov_late),8:14] <- X_nov_late
nNovices <- nNovices
nTrials <- length(which(SONAData$SpeakerLate[1:12]==1))

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

bPrec <- bPrec
nBeta <- nBeta
nExperts <- nExperts
nExperts_agree_early <- nExperts_agree_early
nExperts_disagree_early <- nExperts_disagree_early
nExperts_agree_late <- nExperts_agree_late
nExperts_disagree_late <- nExperts_disagree_late
index_agree_early <- index_agree_early
index_disagree_early <- index_disagree_early
index_agree_late <- index_agree_late
index_disagree_late <- index_disagree_late

I <- diag(2) # identity matrix for Wishart distribution

data <- list("X_agree_early14","X_agree_late14","Y_agree_early","Y_agree_late",
             "Participant_list_agree_early", "Participant_list_agree_late",
             "X_disagree_early14","X_disagree_late14","Y_disagree_early","Y_disagree_late",
             "Participant_list_disagree_early","Participant_list_disagree_late",
             "Prior_early","Prior_late","bPrec","nBeta","I",
             "nAgree_early","nAgree_late","nDisagree_early","nDisagree_late",
             "index_agree_early","index_disagree_early","index_agree_late","index_disagree_late",
             "nExperts","nExperts_agree_early","nExperts_disagree_early","nExperts_agree_late","nExperts_disagree_late",
             "nNovices","X_nov_early14","X_nov_late14","Y_nov_early","Y_nov_late","nTrials") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_exp_early","sigma_con_early","sigma_nov_early",
                "sigma_exp_late","sigma_con_late","sigma_nov_late",
                "phi_disagree","psi_disagree","psi_disagree_bern_early","psi_disagree_bern_late",
                "alpha_early","alpha_late","sigma_psi","rho","r")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_contaminant_14_corrAlpha.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

Beta_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$BETA

alpha_early_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$alpha_early
alpha_late_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$alpha_late

rho_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$rho
r_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$r[,1,2]

Phi_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$phi_disagree

phi_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$psi_disagree

psi_All_corrAlpha_samples_early <- samples$BUGSoutput$sims.list$psi_disagree_bern_early
psi_All_corrAlpha_samples_late <- samples$BUGSoutput$sims.list$psi_disagree_bern_late

SigmaExp_All_corrAlpha_early_samples <- samples$BUGSoutput$sims.list$sigma_exp_early
SigmaCon_All_corrAlpha_early_samples <- samples$BUGSoutput$sims.list$sigma_con_early
SigmaNov_All_corrAlpha_early_samples <- samples$BUGSoutput$sims.list$sigma_nov_early

SigmaExp_All_corrAlpha_late_samples <- samples$BUGSoutput$sims.list$sigma_exp_late
SigmaCon_All_corrAlpha_late_samples <- samples$BUGSoutput$sims.list$sigma_con_late
SigmaNov_All_corrAlpha_late_samples <- samples$BUGSoutput$sims.list$sigma_nov_late

Tau_All_corrAlpha_samples <- samples$BUGSoutput$sims.list$tau

summary_All_corrAlpha <- samples$BUGSoutput$summary


# Main plot ---------------------------------------------------------------

plotValues_early <- matrix(NaN, nrow = 7,ncol = 3)
plotValues_late <- matrix(NaN, nrow = 7,ncol = 3)

loconf = nrow(Beta_All_corrAlpha_samples) * .025
hiconf = nrow(Beta_All_corrAlpha_samples) * .975

tempsort <- sort(Beta_All_corrAlpha_samples[,1])
plotValues_early[1,1] <- mean(tempsort[1500:1001])
plotValues_early[1,2] <- tempsort[loconf]
plotValues_early[1,3] <- tempsort[hiconf]

tempsort <- sort(Beta_All_corrAlpha_samples[,8])
plotValues_late[1,1] <- mean(tempsort[1500:1001])
plotValues_late[1,2] <- tempsort[loconf]
plotValues_late[1,3] <- tempsort[hiconf]

for(i in 2:7){
  
  tempsort <- sort(Beta_All_corrAlpha_samples[,i])
  plotValues_early[i,1] <- mean(tempsort[1500:1001])
  plotValues_early[i,2] <- tempsort[loconf]
  plotValues_early[i,3] <- tempsort[hiconf]
  
  late_ticker = i + 7
  tempsort <- sort(Beta_All_corrAlpha_samples[,late_ticker])
  plotValues_late[i,1] <- mean(tempsort[1500:1001])
  plotValues_late[i,2] <- tempsort[loconf]
  plotValues_late[i,3] <- tempsort[hiconf]
}

plot(1:6,plotValues_early[2:7,1],ylim = c(.15,.9),main = "Yes Early",xlab = 'beta',ylab = 'Interpretation')
abline(h = plotValues_early[1,1],col = '#DC143C')
abline(h= plotValues_early[1,2],col = '#DC143C',lty = 2)
abline(h= plotValues_early[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = plotValues_early[i+1,2],
         x1 = i, y1 = plotValues_early[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}
legend(3.5,.45,legend = c('Novice Interpretation','Expert Interpretation'),pch = c(NA,1),lty = c(1,NA),
       col = c('#DC143C','black'),cex = .8)

plot(1:6,plotValues_late[2:7,1],ylim = c(.15,.9),main = "Yes Late",xlab = 'beta',ylab = 'Interpretation')
abline(h = plotValues_late[1,1],col = '#DC143C')
abline(h= plotValues_late[1,2],col = '#DC143C',lty = 2)
abline(h= plotValues_late[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = plotValues_late[i+1,2],
         x1 = i, y1 = plotValues_late[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}


par(mfrow=c(1,1))
par(mar = c(3.5, 3, 1, 2))

plot1 = MainPlot(plotValues_early[1,1:3],SONA_empirical_Early_mean,SONA_empirical_bycomp_Early,
                 plotValues_early[2:7,1:3],Online_empirical_bycomp_Early,'',
                 if_legend = 1,leftmost = 1)
par(mar = c(3.5, 3, 1, 2))
plot2 = MainPlot(plotValues_late[1,1:3],SONA_empirical_Late_mean,SONA_empirical_bycomp_Late,
                 plotValues_late[2:7,1:3],Online_empirical_bycomp_Late,'',
                 if_legend = 0, leftmost = 0)

# Make a fancy plot where color indicates psi_bern ---------------------------
par(mar = c(4.5, 4, 1, 1))

X_prior_graph_latentMix_psi <- cbind(Prior_early,
                                     summary_All_splitAlpha[230:(229+sum(nDisagree_early)),1])

disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi),ncol = 5)

disagree_plot[,1] = X_prior_graph_latentMix_psi[,1]
disagree_plot[,2] = Y_disagree_early
disagree_plot[,3] = X_prior_graph_latentMix_psi[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi)){
  disagree_plot[i,4] = which(X_disagree_early[i,]==1)-1
  disagree_plot[i,5] = which(X_disagree_early[i,]==1)-1
  if(disagree_plot[i,4]>3){
    disagree_plot[i,5] <- 7-disagree_plot[i,4]
  }
}


plot(disagree_plot[,1],disagree_plot[,2],col=rgb(disagree_plot[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation", main= "early 14b")

abline(a=0,b=1,lty = 3)

text(disagree_plot[,1],disagree_plot[,2],label = disagree_plot[,4],col=rgb(disagree_plot[,3],0,0))




X_prior_graph_latentMix_psi_late <- cbind(Prior_late,
                                          summary_All_splitAlpha[354:(353+sum(nDisagree_late)),1])

disagree_plot_late = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi_late),ncol = 5)

disagree_plot_late[,1] = X_prior_graph_latentMix_psi_late[,1]
disagree_plot_late[,2] = Y_disagree_late
disagree_plot_late[,3] = X_prior_graph_latentMix_psi_late[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi_late)){
  disagree_plot_late[i,4] = which(X_disagree_late[i,]==1)-1
  disagree_plot_late[i,5] = which(X_disagree_late[i,]==1)-1
  if(disagree_plot_late[i,4]>3){
    disagree_plot_late[i,5] <- 7-disagree_plot_late[i,4]
  }
}


plot(disagree_plot_late[,1],disagree_plot_late[,2],col=rgb(disagree_plot_late[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation",type = 'n', main = "late 14b")

text(disagree_plot_late[,1],disagree_plot_late[,2],label = disagree_plot_late[,4],col=rgb(disagree_plot_late[,3],0,0))

abline(a=0,b=1,lty = 3)

# By participant (watch out) -----------------------------------------------
par(mfrow=c(6,4))

for(i in 1:nExperts){
  if(length(which(Participant_list_disagree_early==i))>1){
    temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
    plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
    abline(a=0,b=1,lty = 3)
  }
  if(length(which(Participant_list_disagree_early==i))==1){
    temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
    plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
    abline(a=0,b=1,lty = 3)
  }
}

for(i in 1:nExperts){
  if(length(which(Participant_list_disagree_late==i))>1){
    temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
    plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
    abline(a=0,b=1,lty = 3)
  }
  if(length(which(Participant_list_disagree_late==i))==1){
    temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
    plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
         xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
    
    text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
    abline(a=0,b=1,lty = 3)
  }
}


c1 <- rgb(173,216,230,max = 255, alpha = 130, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 100, names = "lt.pink")

b <- min(c(alpha_diff_exp,alpha_diff_nov)) - 0.001 # Set the minimum for the breakpoints
e <- max(c(alpha_diff_exp,alpha_diff_nov)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 12) # Make a neat vector for the breakpoints
hist1 <- hist(alpha_diff_exp)
hist2 <- hist(alpha_diff_nov)
par(mfrow=c(1,1))
par(mar = c(5, 4, 1, 2))
plot(hist1,col = c1,ylim = c(0,30),xlim = c(-.5,.5),xlab = 'Alpha_late - Alpha_early',main = "Shift",
     axes = FALSE)
plot(hist2,col = c2,add=TRUE)
abline(v=0,lty=3)
axis(side = 2, pos= -.5)
axis(side = 1, pos=0)
legend(-.4,15,legend = c('expert','novice'),pch = c(15,15),
       col = c(c1,c2),cex = 1.3)


# Investigating Results ---------------------------------------------------

latent_results <- data.frame(matrix(NaN,nrow = (length(Y_disagree_early)+length(Y_disagree_late)), ncol = 17))

col_names <- c("Prior","Interpretation","Participant","Composition","Early"
               ,"psi_Sep","phi_Sep","alpha_Sep","beta_Sep",
               "psi_Joint","phi_Joint","alpha_Joint","beta_Joint",
               "psi_alphaSep","phi_alphaSep","alpha_alphaSep","beta_alphaSep")

colnames(latent_results) <- col_names

latent_results$Prior[1:sum(nDisagree_early)] <- Prior_early
latent_results$Prior[(sum(nDisagree_early)+1):nrow(latent_results)] <- Prior_late

latent_results$Interpretation[1:sum(nDisagree_early)] <- Y_disagree_early
latent_results$Interpretation[(sum(nDisagree_early)+1):nrow(latent_results)] <- Y_disagree_late

latent_results$Participant[1:sum(nDisagree_early)] <- Participant_list_disagree_early
latent_results$Participant[(sum(nDisagree_early)+1):nrow(latent_results)] <- Participant_list_disagree_late

latent_results$Composition[1:sum(nDisagree_early)] <- disagree_plot[,4]
latent_results$Composition[(sum(nDisagree_early)+1):nrow(latent_results)] <- disagree_plot_late[,4]

latent_results$psi_Sep[1:sum(nDisagree_early)] <- colMeans(psi_Early_samples)
latent_results$psi_Sep[(sum(nDisagree_early)+1):nrow(latent_results)] <- colMeans(psi_Late_samples)

latent_results$psi_Joint[1:sum(nDisagree_early)] <- colMeans(psi_All_samples_early)
latent_results$psi_Joint[(sum(nDisagree_early)+1):nrow(latent_results)] <- colMeans(psi_All_samples_late)


latent_results$psi_alphaSep[1:sum(nDisagree_early)] <- colMeans(psi_All_splitAlpha_samples_early)
latent_results$psi_alphaSep[(sum(nDisagree_early)+1):nrow(latent_results)] <- colMeans(psi_All_splitAlpha_samples_late)

individually_Betas_early <- colMeans(Beta_Early_samples[,2:7])
individually_Betas_late <- colMeans(Beta_Late_samples[,2:7])
individually_alphas_early <- colMeans(alpha_Early_samples[,1:nExperts])
individually_alphas_late <-colMeans(alpha_Late_samples[,1:nExperts])
individually_phis_early <- colMeans(phi_Early_samples)
individually_phis_late <- colMeans(phi_Late_samples)

joint_phis <- colMeans(phi_All_samples)
joint_alphas <- colMeans(alpha_All_samples[,1:nExperts])
joint_Betas_early <- colMeans(Beta_All_samples[,2:7])
joint_Betas_late <- colMeans(Beta_All_samples[,9:14])

splitAlpha_phis <- colMeans(phi_All_splitAlpha_samples)
splitAlpha_Betas_early <- colMeans(Beta_All_splitAlpha_samples[,2:7])
splitAlpha_Betas_late <- colMeans(Beta_All_splitAlpha_samples[,9:14])
splitAlpha_alpha_early <- colMeans(alpha_early_All_splitAlpha_samples[,1:nExperts])
splitAlpha_alpha_late <- colMeans(alpha_late_All_splitAlpha_samples[,1:nExperts])


for(i in 1:sum(nDisagree_early)){
  latent_results$beta_Sep[i] <- individually_Betas_early[latent_results$Composition[i]]
  latent_results$alpha_Sep[i] <- individually_alphas_early[latent_results$Participant[i]]
  latent_results$phi_Sep[i] <- individually_phis_early[latent_results$Participant[i]]
  
  latent_results$beta_Joint[i] <- joint_Betas_early[latent_results$Composition[i]]
  latent_results$alpha_Joint[i] <- joint_alphas[latent_results$Participant[i]]
  latent_results$phi_Joint[i] <- joint_phis[latent_results$Participant[i]]
  
  latent_results$beta_alphaSep[i] <- splitAlpha_Betas_early[latent_results$Composition[i]]
  latent_results$alpha_alphaSep[i] <- splitAlpha_alpha_early[latent_results$Participant[i]]
  latent_results$phi_alphaSep[i] <- splitAlpha_phis[latent_results$Participant[i]]
}

for(i in 1:sum(nDisagree_late)){
  ticker = i+sum(nDisagree_early)
  latent_results$beta_Sep[ticker] <- individually_Betas_late[latent_results$Composition[ticker]]
  latent_results$alpha_Sep[ticker] <- individually_alphas_late[latent_results$Participant[ticker]]
  latent_results$phi_Sep[ticker] <- individually_phis_late[latent_results$Participant[ticker]]
  
  latent_results$beta_Joint[ticker] <- joint_Betas_late[latent_results$Composition[ticker]]
  latent_results$alpha_Joint[ticker] <- joint_alphas[latent_results$Participant[ticker]]
  latent_results$phi_Joint[ticker] <- joint_phis[latent_results$Participant[ticker]]
  
  latent_results$beta_alphaSep[ticker] <- splitAlpha_Betas_late[latent_results$Composition[ticker]]
  latent_results$alpha_alphaSep[ticker] <- splitAlpha_alpha_late[latent_results$Participant[ticker]]
  latent_results$phi_alphaSep[ticker] <- splitAlpha_phis[latent_results$Participant[ticker]]
}

latent_results$Early[1:sum(nDisagree_early)] <- 1

latent_results$Early[(sum(nDisagree_early)+1):nrow(latent_results)] <- 0

write_csv(latent_results,"latent_result.csv")

# Pick a model ------------------------------------------------------------

joint_Betas_early - individually_Betas_early
joint_Betas_late - individually_Betas_late

plot(latent_results$phi_alphaSep,latent_results$phi_Joint)
abline(a=0,b=1)