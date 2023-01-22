# clears workspace:  
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)

source("BootstrapMean95.R")
source("PlotFunc.R")

# Load data --------------------------------------------------------

## Each row of data is a participant responding to a new composition. They each saw
## 12 compositions and interpreted two generalizations (early & late) for each
## composition. The hypothetical speaker said they would endorse half of each 
## but wouldn't endorse the other half. We're only interested in the statements
## that speakers WOULD endorse (SpeakerEarly or SpeakerLate = 1). Thus there are
## 6 statements of interest for each type of statement (early & late), evenly
## distributed so that there's one for each typ of composition.

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

# Run early condition in JAGS (add phi) -------------------------------------------
nSamples = 4000

DIC <- rep(NaN,6)

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

Beta_Early_add_phi_samples <- samples$BUGSoutput$sims.list$BETA

alpha_Early_add_phi_samples <- samples$BUGSoutput$sims.list$alpha

phi_Early_add_phi_samples <- samples$BUGSoutput$sims.list$psi_disagree

psi_Early_add_phi_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

Phi_Early_add_phi_samples <- samples$BUGSoutput$sims.list$phi_disagree

SigmaExp_Early_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Early_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Early_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Early_add_phi_samples <- samples$BUGSoutput$sims.list$tau

summary_add_phi_Early <- samples$BUGSoutput$summary

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

# Run late condition in JAGS (add phi) ----------------------------------------------

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

Beta_Late_add_phi_samples <- samples$BUGSoutput$sims.list$BETA

alpha_Late_add_phi_samples <- samples$BUGSoutput$sims.list$alpha

phi_Late_add_phi_samples <- samples$BUGSoutput$sims.list$psi_disagree

Phi_Late_add_phi_samples <- samples$BUGSoutput$sims.list$phi_disagree

psi_Late_add_phi_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

SigmaExp_Late_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Late_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Late_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Late_add_phi_samples <- samples$BUGSoutput$sims.list$tau

summary_add_phi_Late <- samples$BUGSoutput$summary

DIC[4] <- samples$BUGSoutput$DIC

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

# Make a fancy plot where color indicates psi_bern (add phi) ---------------------------
par(mar = c(4.5, 4, 1, 1))

X_prior_graph_latentMix_psi <- cbind(Prior_early,
                                     colMeans(psi_Early_add_phi_samples))

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
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation", main= "early add phi")

abline(a=0,b=1,lty = 3)

text(disagree_plot[,1],disagree_plot[,2],label = disagree_plot[,4],col=rgb(disagree_plot[,3],0,0))


X_prior_graph_latentMix_psi_late <- cbind(Prior_late,
                                          colMeans(psi_Late_add_phi_samples))

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
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation",type = 'n', main = "late add phi")

text(disagree_plot_late[,1],disagree_plot_late[,2],label = disagree_plot_late[,4],col=rgb(disagree_plot_late[,3],0,0))

abline(a=0,b=1,lty = 3)

# Compare expert and novice sigmas ----------------------------------------

SigmaTest_yesEarly <- matrix(NaN,nrow = length(SigmaExp_Early_add_phi_samples),ncol = 2)
SigmaTest_yesLate <- matrix(NaN,nrow = length(SigmaExp_Late_add_phi_samples),ncol = 2)

SigmaDiff_yesEarly <- rep(NaN,length(SigmaExp_Early_add_phi_samples))
SigmaDiff_yesLate <- rep(NaN,length(SigmaExp_Late_add_phi_samples))

for(i in 1:length(SigmaExp_Early_add_phi_samples)){
  SigmaTest_yesEarly[i,1] <- Tau_Early_add_phi_samples[i]^2/
    (Tau_Early_add_phi_samples[i]^2 + SigmaExp_Early_add_phi_samples[i]^2)
  SigmaTest_yesEarly[i,2] <- Tau_Early_add_phi_samples[i]^2/
    (Tau_Early_add_phi_samples[i]^2 + SigmaNov_Early_add_phi_samples[i]^2)
  
  SigmaDiff_yesEarly[i] <- SigmaExp_Early_add_phi_samples[i]-SigmaNov_Early_add_phi_samples[i]
  SigmaDiff_yesLate[i] <- SigmaExp_Late_add_phi_samples[i]-SigmaNov_Late_add_phi_samples[i]
  
  SigmaTest_yesLate[i,1] <- Tau_Late_add_phi_samples[i]^2/
    (Tau_Late_add_phi_samples[i]^2 + SigmaExp_Late_add_phi_samples[i]^2)
  SigmaTest_yesLate[i,2] <- Tau_Late_add_phi_samples[i]^2/
    (Tau_Late_add_phi_samples[i]^2 + SigmaNov_Late_add_phi_samples[i]^2)
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

sortSigmaTest_yesEarly <- matrix(NaN,nrow = length(SigmaExp_Early_add_phi_samples),ncol = 2)
sortSigmaTest_yesLate <- matrix(NaN,nrow = length(SigmaExp_Late_add_phi_samples),ncol = 2)

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


# Plot beta inferences (add phi) -------------------------------------------
plotValues_early <- matrix(NaN, nrow = 7,ncol = 3)
plotValues_late <- matrix(NaN, nrow = 7,ncol = 3)

loconf = nrow(Beta_Early_add_phi_samples) * .025
hiconf = nrow(Beta_Early_add_phi_samples) * .975

tempsort <- sort(Beta_Early_add_phi_samples[,1])
plotValues_early[1,1] <- mean(tempsort[1500:1001])
plotValues_early[1,2] <- tempsort[loconf]
plotValues_early[1,3] <- tempsort[hiconf]

tempsort <- sort(Beta_Late_add_phi_samples[,1])
plotValues_late[1,1] <- mean(tempsort[1500:1001])
plotValues_late[1,2] <- tempsort[loconf]
plotValues_late[1,3] <- tempsort[hiconf]

for(i in 2:7){
  
  tempsort <- sort(Beta_Early_add_phi_samples[,i])
  plotValues_early[i,1] <- mean(tempsort[1500:1001])
  plotValues_early[i,2] <- tempsort[loconf]
  plotValues_early[i,3] <- tempsort[hiconf]
  
  
  tempsort <- sort(Beta_Late_add_phi_samples[,i])
  plotValues_late[i,1] <- mean(tempsort[1500:1001])
  plotValues_late[i,2] <- tempsort[loconf]
  plotValues_late[i,3] <- tempsort[hiconf]
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
                 if_legend = 1,leftmost = 1,bottommost = 0)
par(mar = c(3.5, 3, 1, 2))
plot2 = MainPlot(plotValues_late[1,1:3],SONA_empirical_Late_mean,SONA_empirical_bycomp_Late,
                 plotValues_late[2:7,1:3],Online_empirical_bycomp_Late,'',
                 if_legend = 0, leftmost = 1,bottommost = 1)

# Fancy latent mixture plot for publication -------------------------------

X_prior_graph_latentMix_psi_early <- cbind(Prior_early,
                                     colMeans(psi_Early_add_phi_samples))

disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi_early),ncol = 5)

disagree_plot[,1] = X_prior_graph_latentMix_psi_early[,1]
disagree_plot[,2] = Y_disagree_early
disagree_plot[,3] = X_prior_graph_latentMix_psi_early[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi_early)){
  disagree_plot[i,4] = which(X_disagree_early[i,]==1)-1
  disagree_plot[i,5] = which(X_disagree_early[i,]==1)-1
  if(disagree_plot[i,4]>3){
    disagree_plot[i,5] <- 7-disagree_plot[i,4]
  }
}

par(mfrow=c(2,3))
par(mar = c(3.3, 3.3, 1, 1))

early_titles <- c("Early (E+)","Early (E0)","Early (E-)","Early (L+)","Early (L0)","Early (L-)")
late_titles <- c("Late (E+)","Late (E0)","Late (E-)","Late (L+)","Late (L0)","Late (L-)")

for(i in 1:nTrials){
  id_comp <- i+1
  plot_index <- which(disagree_plot[,4]==i)
  
  plot(disagree_plot[plot_index,1],disagree_plot[plot_index,2],col = rgb(0,0,0),pch = 21,
       bg=rgb(disagree_plot[plot_index,3],disagree_plot[plot_index,3],disagree_plot[plot_index,3]),
       xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ", main= " ",
       axes = FALSE,cex = 2,type="n")
  
  lines(x=c(0,1),y=c(0,1),lty = 3,lwd =2)
  
  temp_lines_early <- Beta_Early_add_phi_samples[,id_comp]
  
  text(x = .7, y = 1,labels = early_titles[i], cex = 2)
  
  for(j in 1:length(temp_lines_early)){
    lines(x = c(0,1), y = c(temp_lines_early[j],temp_lines_early[j]),col = rgb(0, 0, 0, max = 255, alpha = 1.2))
  }
  
  points(disagree_plot[plot_index,1],disagree_plot[plot_index,2],col = rgb(0,0,0),pch = 21,
         bg=rgb(disagree_plot[plot_index,3],disagree_plot[plot_index,3],disagree_plot[plot_index,3]),cex = 2)
  
  axis(side=1,pos=0,cex.axis = 1.6,lwd = 2)
  axis(side = 2, pos=0,cex.axis = 1.6,lwd = 2)
  
  if(i == 1){
    title(ylab="Interpretation ( /100 games)", line=1.7, cex.lab=1.8)
  }
  
  if(i == 4){
    title(ylab="Interpretation ( /100 games)", line=1.7, cex.lab=1.8)
  }
  
  if(i > 3){
    title(xlab="Prior ( /100 games)", line=1.7, cex.lab=1.8)
  }
  
  if(i == 6){
    legend(x = .6, y = .5,col = c("black","black","darkgray"), cex = 1.7,lty = c(NaN,NaN,1),
           lwd = c(NaN,NaN,3), pch = c(21,21,NaN), pt.bg = c("black","white","white"), bty = 'n',
           legend = c(expression(paste(psi, " = ", 0)),
                      expression(paste(psi, " = ", 1)),
                      expression(beta[j])))
  }
  
}

X_prior_graph_latentMix_psi_late <- cbind(Prior_late,
                                           colMeans(psi_Late_add_phi_samples))

disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi_late),ncol = 5)

disagree_plot[,1] = X_prior_graph_latentMix_psi_late[,1]
disagree_plot[,2] = Y_disagree_late
disagree_plot[,3] = X_prior_graph_latentMix_psi_late[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi_late)){
  disagree_plot[i,4] = which(X_disagree_late[i,]==1)-1
  disagree_plot[i,5] = which(X_disagree_late[i,]==1)-1
  if(disagree_plot[i,4]>3){
    disagree_plot[i,5] <- 7-disagree_plot[i,4]
  }
}

par(mfrow=c(2,3))
par(mar = c(3.3, 3.3, 1, 1))

for(i in 1:nTrials){
  id_comp <- i+1
  plot_index <- which(disagree_plot[,4]==i)
  
  plot(disagree_plot[plot_index,1],disagree_plot[plot_index,2],col = rgb(0,0,0),pch = 21,
       bg=rgb(disagree_plot[plot_index,3],disagree_plot[plot_index,3],disagree_plot[plot_index,3]),
       xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ", main= " ",
       axes = FALSE,cex = 2,type = "n")
  
  lines(x=c(0,1),y=c(0,1),lty = 3,lwd =2)
  
  text(x = .7, y = 1,labels = late_titles[i], cex = 2)

  temp_lines_late <- Beta_Late_add_phi_samples[,id_comp]
  
  for(j in 1:length(temp_lines_late)){
    lines(x = c(0,1), y = c(temp_lines_late[j],temp_lines_late[j]),col = rgb(0, 0, 0, max = 255, alpha = 1.2))
  }
  
  points(disagree_plot[plot_index,1],disagree_plot[plot_index,2],col = rgb(0,0,0),pch = 21,
         bg=rgb(disagree_plot[plot_index,3],disagree_plot[plot_index,3],disagree_plot[plot_index,3]), cex = 2)
  
  axis(side=1,pos=0,cex.axis = 1.6,lwd = 2)
  axis(side = 2, pos=0,cex.axis = 1.6,lwd = 2)
  
  if(i == 1){
    title(ylab="Interpretation ( /100 games)", line=1.7, cex.lab=1.8)
  }
  
  if(i == 4){
    title(ylab="Interpretation ( /100 games)", line=1.7, cex.lab=1.8)
  }
  
  if(i > 3){
    title(xlab="Prior ( /100 games)", line=1.7, cex.lab=1.8)
  }
  
  if(i == 6){
    legend(x = .6, y = .5,col = c("black","black","darkgrey"), cex = 1.7,lty = c(NaN,NaN,1),
           lwd = c(NaN,NaN,3), pch = c(21,21,NaN), pt.bg = c("black","white","white"),bty = 'n',
           legend = c(expression(paste(psi, " = ", 0)),
                      expression(paste(psi, " = ", 1)),
                      expression(beta[j])))
  }
  
}

# Prep overfull version ---------------------------------------------------

## First 6 betas are naive participant interpretations by condition, 
## the next 6 are experienced by condition
nBeta = 12
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
X_agree_early_overfull <- matrix(0, nrow = sum(nAgree_early),ncol = nBeta)
Y_agree_early <- rep(NaN,sum(nAgree_early))
Participant_list_agree_early <- rep(NaN,sum(nAgree_early))

X_disagree_early_overfull <- matrix(0, nrow = sum(nDisagree_early),ncol = nBeta)
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
    
    X_agree_early_overfull[tickerAgree_early,(OnlineData$CompType[i]+6)] <- 1
    Participant_list_agree_early[tickerAgree_early] <- OnlineData$Participant[i]
    
    tickerAgree_early = tickerAgree_early+1
  }
  if(OnlineData$AgreeEarly[i]==0){
    Y_disagree_early[tickerDisagree_early] <- OnlineData$EarlyResponse[i]
    Prior_early[tickerDisagree_early] <- OnlineData$EarlyPrior[i]
    
    X_disagree_early_overfull[tickerDisagree_early,(OnlineData$CompType[i]+6)] <- 1
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

X_nov_early_overfull <- matrix(0, nrow = Novice_early, ncol = nBeta)
ticker = 1
Y_nov_early <- rep(NaN, Novice_early)

for(i in 1:nrow(SONAData)){
  if(SONAData$SpeakerEarly[i]==1){
    Y_nov_early[ticker] <- SONAData$EarlyResponse[i]
    X_nov_early_overfull[ticker,SONAData$CompType[i]] <- 1
    ticker = ticker + 1
  }
}

# Run JAGS model for early condition (overfull version) -------------------------------------
nSamples = 4000

X_agree <- X_agree_early_overfull
Y_agree <- Y_agree_early
## Which participant produced each data point
Participant_list_agree <- Participant_list_agree_early
## How many generalizations did each participant agree with
nAgree <- nAgree_early

X_disagree <- X_disagree_early_overfull
Y_disagree <- Y_disagree_early
## Which participant produced each data point
Participant_list_disagree <- Participant_list_disagree_early
## Prior estimates for each disagreement
Prior <- Prior_early
## How many generalizations did each participant disagree with
nDisagree <- nDisagree_early

X_nov <- X_nov_early_overfull
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

SigmaExp_Early_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Early_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Early_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Early_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$tau

summary_add_phi_Early_overfull <- samples$BUGSoutput$summary

DIC[2] <- samples$BUGSoutput$DIC

# Prep overfull version for Late condition---------------------------------------

## First 6 betas are naive participant interpretations by condition, 
## the next 6 are experienced by condition
nBeta = 12
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}


## Set up a count of how many times each expert participant agreed and disagreed
## with the speaker
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

## preallocate stuff
X_agree_late_overfull <- matrix(0, nrow = sum(nAgree_late),ncol = nBeta)
Y_agree_late <- rep(NaN,sum(nAgree_late))
Participant_list_agree_late <- rep(NaN,sum(nAgree_late))

X_disagree_late_overfull <- matrix(0, nrow = sum(nDisagree_late),ncol = nBeta)
Y_disagree_late <- rep(NaN,sum(nDisagree_late))
Participant_list_disagree_late <- rep(NaN,sum(nDisagree_late))
Prior_late <- rep(NaN,sum(nDisagree_late))

tickerAgree_late = 1
tickerDisagree_late = 1

## Make data, containing only examples of expert agreement. Y is their interpretation.
## "Prior" is the X variable (input) for the latent mixture part of the model and is thus
## only needed for disagreements, since that's the only part of the model with the mixture.
## X is the indicator matrix for turning the betas on and off as appropriate in the main
## model. Since these are only ever expert participants, the first column is always 0.
## For each row, there will always be a single column with a 1 corresponding to the condition.

## Participant_list is a vector identifying which participant is responsible for
## each data point.
for(i in 1:nrow(OnlineData)){
  if(OnlineData$AgreeLate[i]==1){
    Y_agree_late[tickerAgree_late] <- OnlineData$LateResponse[i]
    
    X_agree_late_overfull[tickerAgree_late,(OnlineData$CompType[i]+6)] <- 1
    Participant_list_agree_late[tickerAgree_late] <- OnlineData$Participant[i]
    
    tickerAgree_late = tickerAgree_late+1
  }
  if(OnlineData$AgreeLate[i]==0){
    Y_disagree_late[tickerDisagree_late] <- OnlineData$LateResponse[i]
    Prior_late[tickerDisagree_late] <- OnlineData$LatePrior[i]
    
    X_disagree_late_overfull[tickerDisagree_late,(OnlineData$CompType[i]+6)] <- 1
    Participant_list_disagree_late[tickerDisagree_late] <- OnlineData$Participant[i]
    
    tickerDisagree_late = tickerDisagree_late+1
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

## X and Y for novices. X is always just a 1 in the first column because these
## are all novice participants

X_nov_late_overfull <- matrix(0, nrow = Novice_late, ncol = nBeta)
ticker = 1
Y_nov_late <- rep(NaN, Novice_late)

for(i in 1:nrow(SONAData)){
  if(SONAData$SpeakerLate[i]==1){
    Y_nov_late[ticker] <- SONAData$LateResponse[i]
    X_nov_late_overfull[ticker,SONAData$CompType[i]] <- 1
    ticker = ticker + 1
  }
}

# Run JAGS model for late condition (overfull version) -------------------------------------
nSamples = 4000

X_agree <- X_agree_late_overfull
Y_agree <- Y_agree_late
## Which participant produced each data point
Participant_list_agree <- Participant_list_agree_late
## How many generalizations did each participant agree with
nAgree <- nAgree_late

X_disagree <- X_disagree_late_overfull
Y_disagree <- Y_disagree_late
## Which participant produced each data point
Participant_list_disagree <- Participant_list_disagree_late
## Prior estimates for each disagreement
Prior <- Prior_late
## How many generalizations did each participant disagree with
nDisagree <- nDisagree_late

X_nov <- X_nov_late_overfull
Y_nov <- Y_nov_late
nNovices <- nNovices
nTrials <- length(which(SONAData$SpeakerLate[1:12]==1))

bPrec <- bPrec
nBeta <- nBeta
## Total experts
nExperts <- nExperts
## Number of experts who agreed and disagreed with any generalizations
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
                "phi_disagree",
                "psi_disagree","psi_disagree_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(112,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_contaminant.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

SigmaExp_Late_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Late_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Late_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Late_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$tau

summary_add_phi_Late_overfull <- samples$BUGSoutput$summary

DIC[5] <- samples$BUGSoutput$DIC

# Compare expert and novice sigmas ----------------------------------------

SigmaTest_yesEarly <- matrix(NaN,nrow = length(SigmaExp_Early_add_phi_samples_overfull),ncol = 2)
SigmaTest_yesLate <- matrix(NaN,nrow = length(SigmaExp_Late_add_phi_samples_overfull),ncol = 2)

SigmaDiff_yesEarly <- rep(NaN,length(SigmaExp_Early_add_phi_samples_overfull))
SigmaDiff_yesLate <- rep(NaN,length(SigmaExp_Late_add_phi_samples_overfull))

for(i in 1:length(SigmaExp_Early_add_phi_samples_overfull)){
  SigmaTest_yesEarly[i,1] <- Tau_Early_add_phi_samples_overfull[i]^2/
    (Tau_Early_add_phi_samples_overfull[i]^2 + SigmaExp_Early_add_phi_samples_overfull[i]^2)
  SigmaTest_yesEarly[i,2] <- Tau_Early_add_phi_samples_overfull[i]^2/
    (Tau_Early_add_phi_samples_overfull[i]^2 + SigmaNov_Early_add_phi_samples_overfull[i]^2)
  
  SigmaDiff_yesEarly[i] <- SigmaExp_Early_add_phi_samples_overfull[i]-SigmaNov_Early_add_phi_samples_overfull[i]
  SigmaDiff_yesLate[i] <- SigmaExp_Late_add_phi_samples_overfull[i]-SigmaNov_Late_add_phi_samples_overfull[i]
  
  SigmaTest_yesLate[i,1] <- Tau_Late_add_phi_samples_overfull[i]^2/
    (Tau_Late_add_phi_samples_overfull[i]^2 + SigmaExp_Late_add_phi_samples_overfull[i]^2)
  SigmaTest_yesLate[i,2] <- Tau_Late_add_phi_samples_overfull[i]^2/
    (Tau_Late_add_phi_samples_overfull[i]^2 + SigmaNov_Late_add_phi_samples_overfull[i]^2)
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

sortSigmaTest_yesEarly <- matrix(NaN,nrow = length(SigmaExp_Early_add_phi_samples_overfull),ncol = 2)
sortSigmaTest_yesLate <- matrix(NaN,nrow = length(SigmaExp_Late_add_phi_samples_overfull),ncol = 2)

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

# Savage Dickey for sigma difference (overfull) --------------------------------------

SigmaSimPulls1 <- rgamma(1000000, shape = 1.5, rate = 2)
SigmaSimPulls2 <- rgamma(1000000, shape = 1.5, rate = 2)

SigmaSimDiff <- SigmaSimPulls1-SigmaSimPulls2

hist(SigmaSimDiff)

SigmaDiffPrior_mean = mean(SigmaSimDiff)
SigmaDiffPrior_sd = sd(SigmaSimDiff)

prior_density_SigDiff = dnorm(0,0,.867)

posterior_density_sigDiff_yesEarly <- dlogspline(0,logspline(SigmaDiff_yesEarly))
posterior_density_sigDiff_yesLate <- dlogspline(0,logspline(SigmaDiff_yesLate))

SavDick_sig_overfull <- c(posterior_density_sigDiff_yesEarly/prior_density_SigDiff,
                          posterior_density_sigDiff_yesLate/prior_density_SigDiff)

# Prepare early condition data for underfull model ----------------------------

## Treating generlizations about early game and late game as separate statements
## (and thus essentially replications of each other)

## First beta is naive participant interpretation, the next 6 are experienced 
## by condition
nBeta = 2
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
X_agree_early_underfull <- matrix(0, nrow = sum(nAgree_early),ncol = nBeta)
Y_agree_early <- rep(NaN,sum(nAgree_early))
Participant_list_agree_early <- rep(NaN,sum(nAgree_early))

X_disagree_early_underfull <- matrix(0, nrow = sum(nDisagree_early),ncol = nBeta)
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
    
    X_agree_early_underfull[tickerAgree_early,2] <- 1
    Participant_list_agree_early[tickerAgree_early] <- OnlineData$Participant[i]
    
    tickerAgree_early = tickerAgree_early+1
  }
  if(OnlineData$AgreeEarly[i]==0){
    Y_disagree_early[tickerDisagree_early] <- OnlineData$EarlyResponse[i]
    Prior_early[tickerDisagree_early] <- OnlineData$EarlyPrior[i]
    
    X_disagree_early_underfull[tickerDisagree_early,2] <- 1
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

X_nov_early_underfull <- matrix(0, nrow = Novice_early, ncol = nBeta)
X_nov_early_underfull[,1] <- rep(1,Novice_early)
ticker = 1
Y_nov_early <- rep(NaN, Novice_early)

for(i in 1:nrow(SONAData)){
  if(SONAData$SpeakerEarly[i]==1){
    Y_nov_early[ticker] <- SONAData$EarlyResponse[i]
    ticker = ticker + 1
  }
}

# Run early condition in JAGS (underfull) -------------------------------------------
nSamples = 4000

X_agree <- X_agree_early_underfull
Y_agree <- Y_agree_early
## Which participant produced each data point
Participant_list_agree <- Participant_list_agree_early
## How many generalizations did each participant agree with
nAgree <- nAgree_early

X_disagree <- X_disagree_early_underfull
Y_disagree <- Y_disagree_early
## Which participant produced each data point
Participant_list_disagree <- Participant_list_disagree_early
## Prior estimates for each disagreement
Prior <- Prior_early
## How many generalizations did each participant disagree with
nDisagree <- nDisagree_early

X_nov <- X_nov_early_underfull
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
myinits <- list(list("BETA" = runif(2,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_contaminant.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

SigmaExp_Early_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Early_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Early_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Early_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$tau

summary_add_phi_Early_underfull <- samples$BUGSoutput$summary

DIC[3] <- samples$BUGSoutput$DIC

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

# Prepare late condition data for underfull model ----------------------------

## Treating generlizations about late game and late game as separate statements
## (and thus essentially replications of each other)

## First beta is naive participant interpretation, the next 6 are experienced 
## by condition
nBeta = 2
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}


## Set up a count of how many times each expert participant agreed and disagreed
## with the speaker
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

## preallocate stuff
X_agree_late_underfull <- matrix(0, nrow = sum(nAgree_late),ncol = nBeta)
Y_agree_late <- rep(NaN,sum(nAgree_late))
Participant_list_agree_late <- rep(NaN,sum(nAgree_late))

X_disagree_late_underfull <- matrix(0, nrow = sum(nDisagree_late),ncol = nBeta)
Y_disagree_late <- rep(NaN,sum(nDisagree_late))
Participant_list_disagree_late <- rep(NaN,sum(nDisagree_late))
Prior_late <- rep(NaN,sum(nDisagree_late))

tickerAgree_late = 1
tickerDisagree_late = 1

## Make data, containing only examples of expert agreement. Y is their interpretation.
## "Prior" is the X variable (input) for the latent mixture part of the model and is thus
## only needed for disagreements, since that's the only part of the model with the mixture.
## X is the indicator matrix for turning the betas on and off as appropriate in the main
## model. Since these are only ever expert participants, the first column is always 0.
## For each row, there will always be a single column with a 1 corresponding to the condition.

## Participant_list is a vector identifying which participant is responsible for
## each data point.
for(i in 1:nrow(OnlineData)){
  if(OnlineData$AgreeLate[i]==1){
    Y_agree_late[tickerAgree_late] <- OnlineData$LateResponse[i]
    
    X_agree_late_underfull[tickerAgree_late,2] <- 1
    Participant_list_agree_late[tickerAgree_late] <- OnlineData$Participant[i]
    
    tickerAgree_late = tickerAgree_late+1
  }
  if(OnlineData$AgreeLate[i]==0){
    Y_disagree_late[tickerDisagree_late] <- OnlineData$LateResponse[i]
    Prior_late[tickerDisagree_late] <- OnlineData$LatePrior[i]
    
    X_disagree_late_underfull[tickerDisagree_late,2] <- 1
    Participant_list_disagree_late[tickerDisagree_late] <- OnlineData$Participant[i]
    
    tickerDisagree_late = tickerDisagree_late+1
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

## X and Y for novices. X is always just a 1 in the first column because these
## are all novice participants

X_nov_late_underfull <- matrix(0, nrow = Novice_late, ncol = nBeta)
X_nov_late_underfull[,1] <- rep(1,Novice_late)
ticker = 1
Y_nov_late <- rep(NaN, Novice_late)

for(i in 1:nrow(SONAData)){
  if(SONAData$SpeakerLate[i]==1){
    Y_nov_late[ticker] <- SONAData$LateResponse[i]
    ticker = ticker + 1
  }
}

# Run late condition in JAGS (underfull) -------------------------------------------
nSamples = 4000

X_agree <- X_agree_late_underfull
Y_agree <- Y_agree_late
## Which participant produced each data point
Participant_list_agree <- Participant_list_agree_late
## How many generalizations did each participant agree with
nAgree <- nAgree_late

X_disagree <- X_disagree_late_underfull
Y_disagree <- Y_disagree_late
## Which participant produced each data point
Participant_list_disagree <- Participant_list_disagree_late
## Prior estimates for each disagreement
Prior <- Prior_late
## How many generalizations did each participant disagree with
nDisagree <- nDisagree_late

X_nov <- X_nov_late_underfull
Y_nov <- Y_nov_late
nNovices <- nNovices
nTrials <- length(which(SONAData$SpeakerLate[1:12]==1))

bPrec <- bPrec
nBeta <- nBeta
## Total experts
nExperts <- nExperts
## Number of experts who agreed and disagreed with any generalizations
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
                "phi_disagree",
                "psi_disagree","psi_disagree_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(2,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_contaminant.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

SigmaExp_Late_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_Late_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_Late_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$sigma_nov

Tau_Late_add_phi_samples_underfull <- samples$BUGSoutput$sims.list$tau

summary_add_phi_Late_underfull <- samples$BUGSoutput$summary

DIC[6] <- samples$BUGSoutput$DIC

## Test the fancy indexing matrix to make sure it's built properly
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

# Compare B0 and Bj -------------------------------------------------------

prior_density = dnorm(0,0,.3537)
posterior_density_yesEarly <- rep(NaN,6)
posterior_density_yesLate <- rep(NaN,6)

for(i in 1:6){
  tempdif <- Beta_Early_add_phi_samples[,1]-Beta_Early_add_phi_samples[,i+1]
  posterior_density_yesEarly[i] <- dlogspline(0,logspline(tempdif))

  tempdif <- Beta_Late_add_phi_samples[,1]-Beta_Late_add_phi_samples[,i+1]
  posterior_density_yesLate[i] <- dlogspline(0,logspline(tempdif))
}

SavDick <- data.frame(matrix(NaN,nrow=6,ncol=2))
names(SavDick) <- c("Early","Late")

SavDick[,1] <- posterior_density_yesEarly/prior_density
SavDick[,2] <- posterior_density_yesLate/prior_density

# Test for expected order in "Yes" conditions -----------------------------

OrderTestEarly <- matrix(NaN,nrow=nrow(Beta_Early_add_phi_samples),ncol = 7)

OrderTestEarly[,1:6] <- Beta_Early_add_phi_samples[,2:7]

##Early hypothesis: 1 > 2 > 3, 6 > 5 > 4
for(i in 1:nrow(Beta_Early_add_phi_samples)){
  firstcheck <- OrderTestEarly[i,1] > OrderTestEarly[i,2] && 
    OrderTestEarly[i,2] > OrderTestEarly[i,3]
  secondcheck <- OrderTestEarly[i,6] > OrderTestEarly[i,5] && 
    OrderTestEarly[i,5] > OrderTestEarly[i,4]
  OrderTestEarly[i,7] <- firstcheck && secondcheck
}

OrderTestLate <- matrix(NaN,nrow=nrow(Beta_Late_add_phi_samples),ncol = 7)

OrderTestLate[,1:6] <- Beta_Late_add_phi_samples[,2:7]

##Late hypothesis: 3 > 2 > 1, 4 > 5 > 6
for(i in 1:nrow(Beta_Late_add_phi_samples)){
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

# # By participant (watch out) -----------------------------------------------
# par(mfrow=c(6,4))
# 
# for(i in 1:nExperts){
#   if(length(which(Participant_list_disagree_early==i))>1){
#     temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
#     plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
#          xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
#     
#     text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
#     abline(a=0,b=1,lty = 3)
#   }
#   if(length(which(Participant_list_disagree_early==i))==1){
#     temp <- disagree_plot[which(Participant_list_disagree_early==i),1:5]
#     plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
#          xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
#     
#     text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
#     abline(a=0,b=1,lty = 3)
#   }
# }
# 
# for(i in 1:nExperts){
#   if(length(which(Participant_list_disagree_late==i))>1){
#     temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
#     plot(temp[,1],temp[,2],col=rgb(temp[,3],0,0), 
#          xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
#     
#     text(temp[,1],temp[,2],label = temp[,4],col=rgb(temp[,3],0,0))
#     abline(a=0,b=1,lty = 3)
#   }
#   if(length(which(Participant_list_disagree_late==i))==1){
#     temp <- disagree_plot_late[which(Participant_list_disagree_late==i),1:5]
#     plot(temp[1],temp[2],col=rgb(temp[3],0,0), 
#          xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ",type = 'n',main = i)
#     
#     text(temp[1],temp[2],label = temp[4],col=rgb(temp[3],0,0))
#     abline(a=0,b=1,lty = 3)
#   }
# }
# 
# # Run early condition in JAGS (basic) -------------------------------------------
# nSamples = 4000
# 
# X_agree <- X_agree_early
# Y_agree <- Y_agree_early
# ## Which participant produced each data point
# Participant_list_agree <- Participant_list_agree_early
# ## How many generalizations did each participant agree with
# nAgree <- nAgree_early
# 
# X_disagree <- X_disagree_early
# Y_disagree <- Y_disagree_early
# ## Which participant produced each data point
# Participant_list_disagree <- Participant_list_disagree_early
# ## Prior estimates for each disagreement
# Prior <- Prior_early
# ## How many generalizations did each participant disagree with
# nDisagree <- nDisagree_early
# 
# X_nov <- X_nov_early
# Y_nov <- Y_nov_early
# nNovices <- nNovices
# nTrials <- length(which(SONAData$SpeakerEarly[1:12]==1))
# 
# bPrec <- bPrec
# nBeta <- 7
# ## Total experts
# nExperts <- nExperts
# ## Number of experts who agreed and disagreed with any generalizations
# nExperts_agree <- nExperts_agree_early
# nExperts_disagree <- nExperts_disagree_early
# index_agree <- index_agree_early
# index_disagree <- index_disagree_early
# 
# data <- list("X_agree","Y_agree","Participant_list_agree","nAgree", 
#              "X_disagree","Y_disagree","Participant_list_disagree","Prior","nDisagree",
#              "bPrec","nBeta","index_agree","index_disagree",
#              "nExperts","nExperts_agree","nExperts_disagree",
#              "nNovices","X_nov","Y_nov","nTrials") # to be passed on to JAGS
# # parameters to be monitored:	
# parameters <- c("BETA","tau","sigma_exp","sigma_con","sigma_nov","alpha",
#                 "Phi_disagree","psi_disagree_bern","sigma_psi")
# 
# #Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# myinits <- list(list("BETA" = runif(7,-.5,.5)))
# 
# # The following command calls JAGS with specific options.
# #This is for running 1 chain (use code below for faster multiple chains)
# samples <- jags(data, #inits=myinits, 
#                 parameters.to.save=parameters,
#                 model.file="LoLComp_contaminant_basic.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)
# 
# # Now the values for the monitored parameters are in the "samples" object, 
# # ready for inspection.
# 
# Beta_Early_basic_samples <- samples$BUGSoutput$sims.list$BETA
# 
# alpha_Early_basic_samples <- samples$BUGSoutput$sims.list$alpha
# 
# psi_Early_basic_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern
# 
# Phi_Early_basic_samples <- samples$BUGSoutput$sims.list$Phi_disagree
# 
# SigmaExp_Early_basic_samples <- samples$BUGSoutput$sims.list$sigma_exp
# 
# SigmaCon_Early_basic_samples <- samples$BUGSoutput$sims.list$sigma_con
# 
# SigmaNov_Early_basic_samples <- samples$BUGSoutput$sims.list$sigma_nov
# 
# Tau_Early_basic_samples <- samples$BUGSoutput$sims.list$tau
# 
# summary_basic_Early <- samples$BUGSoutput$summary
# 
# DIC[1] <- samples$BUGSoutput$DIC
# 
# 
# ## Test the fancy indexing matrix to make sure it's built properly
# test_vector_agree <- rep(NaN, sum(nAgree_early))
# 
# ticker = 1
# for(i in 1:nExperts_agree){
#   for(j in 1:nAgree[i]){
#     test_vector_agree[ticker] <- index_agree[i,j]
#     ticker = ticker+1
#   }
# }
# 
# test_vector_disagree <- rep(NaN, sum(nDisagree_early))
# 
# ticker = 1
# for(i in 1:nExperts_disagree){
#   for(j in 1:nDisagree[i]){
#     test_vector_disagree[ticker] <- index_disagree[i,j]
#     ticker = ticker+1
#   }
# }
# 
# 
# # Run late condition in JAGS (basic) ----------------------------------------------
# 
# nSamples = 4000
# 
# X_agree <- X_agree_late
# Y_agree <- Y_agree_late
# Participant_list_agree <- Participant_list_agree_late
# nAgree <- nAgree_late
# 
# X_disagree <- X_disagree_late
# Y_disagree <- Y_disagree_late
# Participant_list_disagree <- Participant_list_disagree_late
# Prior <- Prior_late
# nDisagree <- nDisagree_late
# 
# X_nov <- X_nov_late
# Y_nov <- Y_nov_late
# nNovices <- nNovices
# nTrials <- length(which(SONAData$SpeakerLate[1:12]==1))
# 
# bPrec <- bPrec
# nBeta <- nBeta
# nExperts <- nExperts
# nExperts_agree <- nExperts_agree_late
# nExperts_disagree <- nExperts_disagree_late
# index_agree <- index_agree_late
# index_disagree <- index_disagree_late
# 
# data <- list("X_agree","Y_agree","Participant_list_agree","nAgree", 
#              "X_disagree","Y_disagree","Participant_list_disagree","Prior","nDisagree",
#              "bPrec","nBeta","index_agree","index_disagree",
#              "nExperts","nExperts_agree","nExperts_disagree",
#              "nNovices","X_nov","Y_nov","nTrials") # to be passed on to JAGS
# # parameters to be monitored:	
# parameters <- c("BETA","tau","sigma_exp","sigma_con","sigma_nov","alpha",
#                 "Phi_disagree","psi_disagree_bern","sigma_psi")
# 
# #Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# myinits <- list(list("BETA" = runif(7,-.5,.5)))
# 
# # The following command calls JAGS with specific options.
# #This is for running 1 chain (use code below for faster multiple chains)
# samples <- jags(data, #inits=myinits, 
#                 parameters.to.save=parameters,
#                 model.file="LoLComp_contaminant_basic.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)
# 
# # Now the values for the monitored parameters are in the "samples" object, 
# # ready for inspection.
# 
# Beta_Late_basic_samples <- samples$BUGSoutput$sims.list$BETA
# 
# alpha_Late_basic_samples <- samples$BUGSoutput$sims.list$alpha
# 
# Phi_Late_basic_samples <- samples$BUGSoutput$sims.list$Phi_disagree
# 
# psi_Late_basic_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern
# 
# SigmaExp_Late_basic_samples <- samples$BUGSoutput$sims.list$sigma_exp
# 
# SigmaCon_Late_basic_samples <- samples$BUGSoutput$sims.list$sigma_con
# 
# SigmaNov_Late_basic_samples <- samples$BUGSoutput$sims.list$sigma_nov
# 
# Tau_Late_basic_samples <- samples$BUGSoutput$sims.list$tau
# 
# summary_basic_Late <- samples$BUGSoutput$summary
# 
# DIC[3] <- samples$BUGSoutput$DIC
# 
# test_vector_agree <- rep(NaN, sum(nAgree_late))
# 
# ticker = 1
# for(i in 1:nExperts_agree){
#   for(j in 1:nAgree[i]){
#     test_vector_agree[ticker] <- index_agree[i,j]
#     ticker = ticker+1
#   }
# }
# 
# test_vector_disagree <- rep(NaN, sum(nDisagree_late))
# 
# ticker = 1
# for(i in 1:nExperts_disagree){
#   for(j in 1:nDisagree[i]){
#     test_vector_disagree[ticker] <- index_disagree[i,j]
#     ticker = ticker+1
#   }
# }
# 
# # Prep data for unified phi (everything else separated) -------------------
# nSamples = 4000
# 
# nBeta <- 14
# 
# X_agree_early14 <- matrix(0,nrow = sum(nAgree_early), ncol = nBeta)
# X_agree_late14 <- matrix(0,nrow = sum(nAgree_late), ncol = nBeta)
# 
# X_agree_early14[1:nrow(X_agree_early14),1:7] <- X_agree_early
# X_agree_late14[1:nrow(X_agree_late14),8:14] <- X_agree_late
# Y_agree_early <- Y_agree_early
# Y_agree_late <- Y_agree_late
# Participant_list_agree_early <- Participant_list_agree_early
# Participant_list_agree_late <- Participant_list_agree_late
# nAgree_early <- nAgree_early
# nAgree_late <- nAgree_late
# 
# X_disagree_early14 <- matrix(0,nrow = nrow(X_disagree_early), ncol = nBeta)
# X_disagree_late14 <- matrix(0,nrow = nrow(X_disagree_late), ncol = nBeta)
# 
# X_disagree_early14[1:nrow(X_disagree_early),1:7] <- X_disagree_early
# X_disagree_late14[1:nrow(X_disagree_late),8:14] <- X_disagree_late
# Y_disagree_early <- Y_disagree_early
# Y_disagree_late <- Y_disagree_late
# Participant_list_disagree_early <- Participant_list_disagree_early
# Participant_list_disagree_late <- Participant_list_disagree_late
# nDisagree_early <- nDisagree_early
# nDisagree_late <- nDisagree_late
# Prior_early <- Prior_early
# Prior_late <- Prior_late
# 
# X_nov_early14 <- matrix(0,nrow = nrow(X_nov_early), ncol = nBeta)
# X_nov_early14[1:nrow(X_nov_early),1:7] <- X_nov_early
# 
# X_nov_late14 <- matrix(0,nrow = nrow(X_nov_late), ncol = nBeta)
# X_nov_late14[1:nrow(X_nov_late),8:14] <- X_nov_late
# nNovices <- nNovices
# nTrials <- length(which(SONAData$SpeakerLate[1:12]==1))
# 
# bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)
# 
# for(i in 1:nBeta){
#   bPrec[i,i] <- 16
# }
# 
# 
# # Run all in JAGS (unified phi) -------------------------------------------
# 
# bPrec <- bPrec
# nBeta <- nBeta
# nExperts <- nExperts
# nExperts_agree_early <- nExperts_agree_early
# nExperts_disagree_early <- nExperts_disagree_early
# nExperts_agree_late <- nExperts_agree_late
# nExperts_disagree_late <- nExperts_disagree_late
# index_agree_early <- index_agree_early
# index_disagree_early <- index_disagree_early
# index_agree_late <- index_agree_late
# index_disagree_late <- index_disagree_late
# 
# data <- list("X_agree_early14","X_agree_late14","Y_agree_early","Y_agree_late",
#              "Participant_list_agree_early", "Participant_list_agree_late",
#              "X_disagree_early14","X_disagree_late14","Y_disagree_early","Y_disagree_late",
#              "Participant_list_disagree_early","Participant_list_disagree_late",
#              "Prior_early","Prior_late","bPrec","nBeta",
#              "nAgree_early","nAgree_late","nDisagree_early","nDisagree_late",
#              "index_agree_early","index_disagree_early","index_agree_late","index_disagree_late",
#              "nExperts","nExperts_agree_early","nExperts_disagree_early","nExperts_agree_late","nExperts_disagree_late",
#              "nNovices","X_nov_early14","X_nov_late14","Y_nov_early","Y_nov_late","nTrials") # to be passed on to JAGS
# # parameters to be monitored:	
# parameters <- c("BETA","tau","sigma_exp_early","sigma_con_early","sigma_nov_early",
#                 "sigma_exp_late","sigma_con_late","sigma_nov_late",
#                 "Phi_disagree","phi_disagree","psi_disagree_bern_early","psi_disagree_bern_late",
#                 "alpha_early","alpha_late","sigma_phi")
# 
# #Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# myinits <- list(list("BETA" = runif(14,-.5,.5)))
# 
# # The following command calls JAGS with specific options.
# #This is for running 1 chain (use code below for faster multiple chains)
# samples <- jags(data, #inits=myinits, 
#                 parameters.to.save=parameters,
#                 model.file="LoLComp_contaminant_14_splitAlpha.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)
# 
# # Now the values for the monitored parameters are in the "samples" object, 
# # ready for inspection.
# 
# Beta_uniphi_samples <- samples$BUGSoutput$sims.list$BETA
# 
# alpha_early_uniphi_samples <- samples$BUGSoutput$sims.list$alpha_early
# alpha_late_uniphi_samples <- samples$BUGSoutput$sims.list$alpha_late
# 
# Phi_uniphi_samples <- samples$BUGSoutput$sims.list$Phi_disagree
# 
# phi_uniphi_samples <- samples$BUGSoutput$sims.list$phi_disagree
# 
# psi_uniphi_samples_early <- samples$BUGSoutput$sims.list$psi_disagree_bern_early
# psi_uniphi_samples_late <- samples$BUGSoutput$sims.list$psi_disagree_bern_late
# 
# SigmaExp_uniphi_early_samples <- samples$BUGSoutput$sims.list$sigma_exp_early
# SigmaCon_uniphi_early_samples <- samples$BUGSoutput$sims.list$sigma_con_early
# SigmaNov_uniphi_early_samples <- samples$BUGSoutput$sims.list$sigma_nov_early
# 
# SigmaExp_uniphi_late_samples <- samples$BUGSoutput$sims.list$sigma_exp_late
# SigmaCon_uniphi_late_samples <- samples$BUGSoutput$sims.list$sigma_con_late
# SigmaNov_uniphi_late_samples <- samples$BUGSoutput$sims.list$sigma_nov_late
# 
# Tau_uniphi_samples <- samples$BUGSoutput$sims.list$tau
# 
# summary_uniphi <- samples$BUGSoutput$summary
# 
# # Make a fancy plot where color indicates psi_bern (basic) ---------------------------
# par(mar = c(4.5, 4, 1, 1))
# 
# X_prior_graph_latentMix_psi <- cbind(Prior_early,
#                                      colMeans(psi_Early_basic_samples))
# 
# disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi),ncol = 5)
# 
# disagree_plot[,1] = X_prior_graph_latentMix_psi[,1]
# disagree_plot[,2] = Y_disagree_early
# disagree_plot[,3] = X_prior_graph_latentMix_psi[,2]
# 
# for(i in 1:nrow(X_prior_graph_latentMix_psi)){
#   disagree_plot[i,4] = which(X_disagree_early[i,]==1)-1
#   disagree_plot[i,5] = which(X_disagree_early[i,]==1)-1
#   if(disagree_plot[i,4]>3){
#     disagree_plot[i,5] <- 7-disagree_plot[i,4]
#   }
# }
# 
# 
# plot(disagree_plot[,1],disagree_plot[,2],col=rgb(disagree_plot[,3],0,0), 
#      xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation", main= "early basic")
# 
# abline(a=0,b=1,lty = 3)
# 
# text(disagree_plot[,1],disagree_plot[,2],label = disagree_plot[,4],col=rgb(disagree_plot[,3],0,0))
# 
# 
# X_prior_graph_latentMix_psi_late <- cbind(Prior_late,
#                                           colMeans(psi_Late_basic_samples))
# 
# disagree_plot_late = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi_late),ncol = 5)
# 
# disagree_plot_late[,1] = X_prior_graph_latentMix_psi_late[,1]
# disagree_plot_late[,2] = Y_disagree_late
# disagree_plot_late[,3] = X_prior_graph_latentMix_psi_late[,2]
# 
# for(i in 1:nrow(X_prior_graph_latentMix_psi_late)){
#   disagree_plot_late[i,4] = which(X_disagree_late[i,]==1)-1
#   disagree_plot_late[i,5] = which(X_disagree_late[i,]==1)-1
#   if(disagree_plot_late[i,4]>3){
#     disagree_plot_late[i,5] <- 7-disagree_plot_late[i,4]
#   }
# }
# 
# 
# plot(disagree_plot_late[,1],disagree_plot_late[,2],col=rgb(disagree_plot_late[,3],0,0), 
#      xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation",type = 'n', main = "late basic")
# 
# text(disagree_plot_late[,1],disagree_plot_late[,2],label = disagree_plot_late[,4],col=rgb(disagree_plot_late[,3],0,0))
# 
# abline(a=0,b=1,lty = 3)
# 
# # Make a fancy plot where color indicates psi_bern (uniphi) ---------------------------
# par(mar = c(4.5, 4, 1, 1))
# 
# X_prior_graph_latentMix_psi <- cbind(Prior_early,
#                                      colMeans(psi_uniphi_samples_early))
# 
# disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi),ncol = 5)
# 
# disagree_plot[,1] = X_prior_graph_latentMix_psi[,1]
# disagree_plot[,2] = Y_disagree_early
# disagree_plot[,3] = X_prior_graph_latentMix_psi[,2]
# 
# for(i in 1:nrow(X_prior_graph_latentMix_psi)){
#   disagree_plot[i,4] = which(X_disagree_early[i,]==1)-1
#   disagree_plot[i,5] = which(X_disagree_early[i,]==1)-1
#   if(disagree_plot[i,4]>3){
#     disagree_plot[i,5] <- 7-disagree_plot[i,4]
#   }
# }
# 
# 
# plot(disagree_plot[,1],disagree_plot[,2],col=rgb(disagree_plot[,3],0,0), 
#      xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation", main= "early uniphi")
# 
# abline(a=0,b=1,lty = 3)
# 
# text(disagree_plot[,1],disagree_plot[,2],label = disagree_plot[,4],col=rgb(disagree_plot[,3],0,0))
# 
# 
# X_prior_graph_latentMix_psi_late <- cbind(Prior_late,
#                                           colMeans(psi_uniphi_samples_late))
# 
# disagree_plot_late = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi_late),ncol = 5)
# 
# disagree_plot_late[,1] = X_prior_graph_latentMix_psi_late[,1]
# disagree_plot_late[,2] = Y_disagree_late
# disagree_plot_late[,3] = X_prior_graph_latentMix_psi_late[,2]
# 
# for(i in 1:nrow(X_prior_graph_latentMix_psi_late)){
#   disagree_plot_late[i,4] = which(X_disagree_late[i,]==1)-1
#   disagree_plot_late[i,5] = which(X_disagree_late[i,]==1)-1
#   if(disagree_plot_late[i,4]>3){
#     disagree_plot_late[i,5] <- 7-disagree_plot_late[i,4]
#   }
# }
# 
# 
# plot(disagree_plot_late[,1],disagree_plot_late[,2],col=rgb(disagree_plot_late[,3],0,0), 
#      xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation",type = 'n', main = "late uniphi")
# 
# text(disagree_plot_late[,1],disagree_plot_late[,2],label = disagree_plot_late[,4],col=rgb(disagree_plot_late[,3],0,0))
# 
# abline(a=0,b=1,lty = 3)
# 
# # Investigating Results ---------------------------------------------------
# 
# latent_results <- data.frame(matrix(NaN,nrow = (length(Y_disagree_early)+length(Y_disagree_late)), ncol = 19))
# 
# col_names <- c("Prior","Interpretation","Participant","Composition","Early",
#                "Phi_basic","psi_basic","alpha_basic","beta_basic","Phi","phi",
#                "psi","alpha","beta","Phi_uniphi","phi_uniphi","psi_uniphi",
#                "alpha_uniphi","beta_uniphi")
# 
# colnames(latent_results) <- col_names
# 
# latent_results$Prior[1:sum(nDisagree_early)] <- Prior_early
# latent_results$Prior[(sum(nDisagree_early)+1):nrow(latent_results)] <- Prior_late
# 
# latent_results$Interpretation[1:sum(nDisagree_early)] <- Y_disagree_early
# latent_results$Interpretation[(sum(nDisagree_early)+1):nrow(latent_results)] <- Y_disagree_late
# 
# latent_results$Participant[1:sum(nDisagree_early)] <- Participant_list_disagree_early
# latent_results$Participant[(sum(nDisagree_early)+1):nrow(latent_results)] <- Participant_list_disagree_late
# 
# latent_results$Composition[1:sum(nDisagree_early)] <- disagree_plot[,4]
# latent_results$Composition[(sum(nDisagree_early)+1):nrow(latent_results)] <- disagree_plot_late[,4]
# 
# latent_results$psi_basic[1:sum(nDisagree_early)] <- colMeans(psi_Early_basic_samples)
# latent_results$psi_basic[(sum(nDisagree_early)+1):nrow(latent_results)] <- colMeans(psi_Late_basic_samples)
# 
# latent_results$psi[1:sum(nDisagree_early)] <- colMeans(psi_Early_add_phi_samples)
# latent_results$psi[(sum(nDisagree_early)+1):nrow(latent_results)] <- colMeans(psi_Late_add_phi_samples)
# 
# 
# latent_results$psi_uniphi[1:sum(nDisagree_early)] <- colMeans(psi_uniphi_samples_early)
# latent_results$psi_uniphi[(sum(nDisagree_early)+1):nrow(latent_results)] <- colMeans(psi_uniphi_samples_late)
# 
# Betas_early_basic <- colMeans(Beta_Early_basic_samples[,2:7])
# Betas_late_basic <- colMeans(Beta_Late_basic_samples[,2:7])
# alphas_early_basic <- colMeans(alpha_Early_basic_samples[,1:nExperts])
# alphas_late_basic <-colMeans(alpha_Late_basic_samples[,1:nExperts])
# 
# Betas_early_add_phi <- colMeans(Beta_Early_add_phi_samples[,2:7])
# Betas_late_add_phi <- colMeans(Beta_Late_add_phi_samples[,2:7])
# alphas_early_add_phi <- colMeans(alpha_Early_add_phi_samples[,1:nExperts])
# alphas_late_add_phi <-colMeans(alpha_Late_add_phi_samples[,1:nExperts])
# phis_early_add_phi <- colMeans(phi_Early_add_phi_samples)
# phis_late_add_phi <- colMeans(phi_Late_add_phi_samples)
# 
# Betas_early_uniphi <- colMeans(Beta_uniphi_samples[,2:7])
# Betas_late_uniphi <- colMeans(Beta_uniphi_samples[,9:14])
# alphas_early_uniphi <- colMeans(alpha_early_uniphi_samples[,1:nExperts])
# alphas_late_uniphi <- colMeans(alpha_late_uniphi_samples[,1:nExperts])
# phis_uniphi <- colMeans(phi_uniphi_samples)
# 
# 
# for(i in 1:sum(nDisagree_early)){
#   latent_results$beta_basic[i] <- Betas_early_basic[latent_results$Composition[i]]
#   latent_results$alpha_basic[i] <- alphas_early_basic[latent_results$Participant[i]]
#   
#   latent_results$beta[i] <- Betas_early_add_phi[latent_results$Composition[i]]
#   latent_results$alpha[i] <- alphas_early_add_phi[latent_results$Participant[i]]
#   latent_results$phi[i] <- phis_early_add_phi[latent_results$Participant[i]]
#   
#   latent_results$beta_uniphi[i] <- Betas_early_uniphi[latent_results$Composition[i]]
#   latent_results$alpha_uniphi[i] <- alphas_early_uniphi[latent_results$Participant[i]]
#   latent_results$phi_uniphi[i] <- phis_uniphi[latent_results$Participant[i]]
# }
# 
# for(i in 1:sum(nDisagree_late)){
#   latent_results$beta_basic[i] <- Betas_late_basic[latent_results$Composition[i]]
#   latent_results$alpha_basic[i] <- alphas_late_basic[latent_results$Participant[i]]
#   
#   latent_results$beta[i] <- Betas_late_add_phi[latent_results$Composition[i]]
#   latent_results$alpha[i] <- alphas_late_add_phi[latent_results$Participant[i]]
#   latent_results$phi[i] <- phis_late_add_phi[latent_results$Participant[i]]
#   
#   latent_results$beta_uniphi[i] <- Betas_early_uniphi[latent_results$Composition[i]]
#   latent_results$alpha_uniphi[i] <- alphas_late_uniphi[latent_results$Participant[i]]
#   latent_results$phi_uniphi[i] <- phis_uniphi[latent_results$Participant[i]]
# }
# 
# latent_results$Early[1:sum(nDisagree_early)] <- 1
# 
# latent_results$Early[(sum(nDisagree_early)+1):nrow(latent_results)] <- 0
# 
# latent_results$Phi_basic[1:sum(nDisagree_early)] <- mean(Phi_Early_basic_samples)
# latent_results$Phi[1:sum(nDisagree_early)] <- mean(Phi_Early_add_phi_samples)
# 
# latent_results$Phi_basic[(sum(nDisagree_early)+1):nrow(latent_results)] <- mean(Phi_Late_basic_samples)
# latent_results$Phi[(sum(nDisagree_early)+1):nrow(latent_results)] <- mean(Phi_Late_add_phi_samples)
# 
# latent_results$Phi_uniphi <- mean(Phi_uniphi_samples)
# 
# write_csv(latent_results,"latent_result.csv")
# 
# # Plot beta inferences (basic) -------------------------------------------
# plotValues_early <- matrix(NaN, nrow = 7,ncol = 3)
# plotValues_late <- matrix(NaN, nrow = 7,ncol = 3)
# 
# loconf = nrow(Beta_Early_basic_samples) * .025
# hiconf = nrow(Beta_Early_basic_samples) * .975
# 
# tempsort <- sort(Beta_Early_basic_samples[,1])
# plotValues_early[1,1] <- mean(tempsort[1500:1001])
# plotValues_early[1,2] <- tempsort[loconf]
# plotValues_early[1,3] <- tempsort[hiconf]
# 
# tempsort <- sort(Beta_Late_basic_samples[,1])
# plotValues_late[1,1] <- mean(tempsort[1500:1001])
# plotValues_late[1,2] <- tempsort[loconf]
# plotValues_late[1,3] <- tempsort[hiconf]
# 
# for(i in 2:7){
#   
#   tempsort <- sort(Beta_Early_basic_samples[,i])
#   plotValues_early[i,1] <- mean(tempsort[1500:1001])
#   plotValues_early[i,2] <- tempsort[loconf]
#   plotValues_early[i,3] <- tempsort[hiconf]
#   
#   
#   tempsort <- sort(Beta_Late_basic_samples[,i])
#   plotValues_late[i,1] <- mean(tempsort[1500:1001])
#   plotValues_late[i,2] <- tempsort[loconf]
#   plotValues_late[i,3] <- tempsort[hiconf]
# }
# 
# 
# nNovResponses = nNovices*nTrials
# SONA_empirical_Early <- matrix(NaN,nrow = nNovResponses, ncol = 2)
# SONA_empirical_Late <- matrix(NaN,nrow = nNovResponses, ncol = 2)
# 
# SONA_empirical_Early[,1] <- SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)]
# SONA_empirical_Early[,2] <- SONAData$CompType[which(SONAData$SpeakerEarly==1)]
# 
# SONA_empirical_Late[,1] <- SONAData$LateResponse[which(SONAData$SpeakerLate==1)]
# SONA_empirical_Late[,2] <- SONAData$CompType[which(SONAData$SpeakerLate==1)]
# 
# SONA_empirical_Early_mean <- mean(SONA_empirical_Early[,1])
# SONA_empirical_Late_mean <- mean(SONA_empirical_Late[,1])
# 
# SONA_empirical_bycomp_Early <- rep(NaN,6)
# SONA_empirical_bycomp_Late <- rep(NaN,6)
# 
# for(i in 1:max(SONAData$CompType)){
#   SONA_empirical_bycomp_Early[i] <- mean(SONA_empirical_Early[which(SONA_empirical_Early[,2]==i),1])
#   SONA_empirical_bycomp_Late[i] <- mean(SONA_empirical_Late[which(SONA_empirical_Late[,2]==i),1])
# }
# 
# nExpResponses = nExperts*nTrials
# Online_empirical_Early <- matrix(NaN,nrow = nExpResponses, ncol = 2)
# Online_empirical_Late <- matrix(NaN,nrow = nExpResponses, ncol = 2)
# 
# Online_empirical_Early[,1] <- OnlineData$EarlyResponse[which(OnlineData$SpeakerEarly==1)]
# Online_empirical_Early[,2] <- OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]
# 
# Online_empirical_Late[,1] <- OnlineData$LateResponse[which(OnlineData$SpeakerLate==1)]
# Online_empirical_Late[,2] <- OnlineData$CompType[which(OnlineData$SpeakerLate==1)]
# 
# Online_empirical_Early_mean <- mean(Online_empirical_Early[,1])
# Online_empirical_Late_mean <- mean(Online_empirical_Late[,1])
# 
# Online_empirical_bycomp_Early <- rep(NaN,6)
# Online_empirical_bycomp_Late <- rep(NaN,6)
# 
# for(i in 1:max(OnlineData$CompType)){
#   Online_empirical_bycomp_Early[i] <- mean(Online_empirical_Early[which(Online_empirical_Early[,2]==i),1])
#   Online_empirical_bycomp_Late[i] <- mean(Online_empirical_Late[which(Online_empirical_Late[,2]==i),1])
# }
# 
# 
# par(mfrow=c(1,1))
# par(mar = c(3.5, 3, 1, 2))
# 
# plot1 = MainPlot(plotValues_early[1,1:3],SONA_empirical_Early_mean,SONA_empirical_bycomp_Early,
#                  plotValues_early[2:7,1:3],Online_empirical_bycomp_Early,'',
#                  if_legend = 1,leftmost = 1)
# par(mar = c(3.5, 3, 1, 2))
# plot2 = MainPlot(plotValues_late[1,1:3],SONA_empirical_Late_mean,SONA_empirical_bycomp_Late,
#                  plotValues_late[2:7,1:3],Online_empirical_bycomp_Late,'',
#                  if_legend = 0, leftmost = 0)
# 
# # Plot beta inferences (uniphi) -------------------------------------------
# plotValues_early <- matrix(NaN, nrow = 7,ncol = 3)
# plotValues_late <- matrix(NaN, nrow = 7,ncol = 3)
# 
# loconf = nrow(Beta_uniphi_samples) * .025
# hiconf = nrow(Beta_uniphi_samples) * .975
# 
# tempsort <- sort(Beta_uniphi_samples[,1])
# plotValues_early[1,1] <- mean(tempsort[1500:1001])
# plotValues_early[1,2] <- tempsort[loconf]
# plotValues_early[1,3] <- tempsort[hiconf]
# 
# tempsort <- sort(Beta_uniphi_samples[,8])
# plotValues_late[1,1] <- mean(tempsort[1500:1001])
# plotValues_late[1,2] <- tempsort[loconf]
# plotValues_late[1,3] <- tempsort[hiconf]
# 
# for(i in 2:7){
#   
#   tempsort <- sort(Beta_uniphi_samples[,i])
#   plotValues_early[i,1] <- mean(tempsort[1500:1001])
#   plotValues_early[i,2] <- tempsort[loconf]
#   plotValues_early[i,3] <- tempsort[hiconf]
#   
#   j = i + 7
#   
#   tempsort <- sort(Beta_uniphi_samples[,j])
#   plotValues_late[i,1] <- mean(tempsort[1500:1001])
#   plotValues_late[i,2] <- tempsort[loconf]
#   plotValues_late[i,3] <- tempsort[hiconf]
# }
# 
# 
# nNovResponses = nNovices*nTrials
# SONA_empirical_Early <- matrix(NaN,nrow = nNovResponses, ncol = 2)
# SONA_empirical_Late <- matrix(NaN,nrow = nNovResponses, ncol = 2)
# 
# SONA_empirical_Early[,1] <- SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)]
# SONA_empirical_Early[,2] <- SONAData$CompType[which(SONAData$SpeakerEarly==1)]
# 
# SONA_empirical_Late[,1] <- SONAData$LateResponse[which(SONAData$SpeakerLate==1)]
# SONA_empirical_Late[,2] <- SONAData$CompType[which(SONAData$SpeakerLate==1)]
# 
# SONA_empirical_Early_mean <- mean(SONA_empirical_Early[,1])
# SONA_empirical_Late_mean <- mean(SONA_empirical_Late[,1])
# 
# SONA_empirical_bycomp_Early <- rep(NaN,6)
# SONA_empirical_bycomp_Late <- rep(NaN,6)
# 
# for(i in 1:max(SONAData$CompType)){
#   SONA_empirical_bycomp_Early[i] <- mean(SONA_empirical_Early[which(SONA_empirical_Early[,2]==i),1])
#   SONA_empirical_bycomp_Late[i] <- mean(SONA_empirical_Late[which(SONA_empirical_Late[,2]==i),1])
# }
# 
# nExpResponses = nExperts*nTrials
# Online_empirical_Early <- matrix(NaN,nrow = nExpResponses, ncol = 2)
# Online_empirical_Late <- matrix(NaN,nrow = nExpResponses, ncol = 2)
# 
# Online_empirical_Early[,1] <- OnlineData$EarlyResponse[which(OnlineData$SpeakerEarly==1)]
# Online_empirical_Early[,2] <- OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]
# 
# Online_empirical_Late[,1] <- OnlineData$LateResponse[which(OnlineData$SpeakerLate==1)]
# Online_empirical_Late[,2] <- OnlineData$CompType[which(OnlineData$SpeakerLate==1)]
# 
# Online_empirical_Early_mean <- mean(Online_empirical_Early[,1])
# Online_empirical_Late_mean <- mean(Online_empirical_Late[,1])
# 
# Online_empirical_bycomp_Early <- rep(NaN,6)
# Online_empirical_bycomp_Late <- rep(NaN,6)
# 
# for(i in 1:max(OnlineData$CompType)){
#   Online_empirical_bycomp_Early[i] <- mean(Online_empirical_Early[which(Online_empirical_Early[,2]==i),1])
#   Online_empirical_bycomp_Late[i] <- mean(Online_empirical_Late[which(Online_empirical_Late[,2]==i),1])
# }
# 
# 
# par(mfrow=c(1,1))
# par(mar = c(3.5, 3, 1, 2))
# 
# plot1 = MainPlot(plotValues_early[1,1:3],SONA_empirical_Early_mean,SONA_empirical_bycomp_Early,
#                  plotValues_early[2:7,1:3],Online_empirical_bycomp_Early,'',
#                  if_legend = 1,leftmost = 1)
# par(mar = c(3.5, 3, 1, 2))
# plot2 = MainPlot(plotValues_late[1,1:3],SONA_empirical_Late_mean,SONA_empirical_bycomp_Late,
#                  plotValues_late[2:7,1:3],Online_empirical_bycomp_Late,'',
#                  if_legend = 0, leftmost = 0)
