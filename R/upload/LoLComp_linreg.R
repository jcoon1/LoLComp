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

SONAData <- read_csv("SONAData.csv",col_names=TRUE)

nNovices = max(SONAData$Participant)

totalNovice = nItems*nNovices

# Prep data for regression analysis of prior to interpretation, separate intercepts ----------------------------
total_interest <- length(which(OnlineData$SpeakerEarly==1)) + length(which(OnlineData$SpeakerLate==1))
X_prior <- matrix(0,nrow = total_interest, ncol = 4)
Y_prior <- rep(NaN,total_interest)
Participant_prior <- rep(NaN,total_interest)

ticker = 1
for(i in 1:nrow(OnlineData)){
  if(OnlineData$SpeakerEarly[i]==1){
    if(OnlineData$EarlyTVJ[i]==1){
      X_prior[ticker,1] <- 1
      X_prior[ticker,3] <- OnlineData$EarlyPrior[i]
    }
    else{
      X_prior[ticker,2] <- 1
      X_prior[ticker,4] <- OnlineData$EarlyPrior[i]
    }
    Y_prior[ticker] <- OnlineData$EarlyResponse[i]
    Participant_prior[ticker] <- OnlineData$Participant[i]
    ticker = ticker+1
  }
  
  if(OnlineData$SpeakerLate[i]==1){
    if(OnlineData$LateTVJ[i]==1){
      X_prior[ticker,1] <- 1
      X_prior[ticker,3] <- OnlineData$LatePrior[i]
    }
    else{
      X_prior[ticker,2] <- 1
      X_prior[ticker,4] <- OnlineData$LatePrior[i]
    }
    Y_prior[ticker] <- OnlineData$LateResponse[i]
    Participant_prior[ticker] <- OnlineData$Participant[i]
    ticker = ticker+1
  }
}

# Run JAGS model for regression analysis of prior to interpretation (basic version) -------

nSamples = 4000

nBeta = 4

X <- X_prior
Y<- Y_prior
Participant_list <- Participant_prior
nParticipants <- max(Participant_list)
total_priors <- length(Y_prior)

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

data <- list("X", "Y","bPrec","nBeta","total_priors","Participant_list","nParticipants") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(nBeta,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_prior_to_int.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

priorBeta_samples <- samples$BUGSoutput$sims.list$BETA

priorsigma_samples <- samples$BUGSoutput$sims.list$sigma

prioralpha_samples <- samples$BUGSoutput$sims.list$alpha

priortau_samples <- samples$BUGSoutput$sims.list$tau

prior_summary <- samples$BUGSoutput$summary

prior_DIC <- samples$BUGSoutput$DIC

# # Graph prior to interpretation results -----------------------------------
# 
# agree_plot = matrix(NaN,nrow = length(which(X_prior[,1]==1)),ncol = 2)
# 
# agree_plot[,1] = X_prior[which(X_prior[,1]==1),3]
# agree_plot[,2] = Y_prior[which(X_prior[,1]==1)]
# 
# disagree_plot = matrix(NaN,nrow = length(which(X_prior[,2]==1)),ncol = 2)
# 
# disagree_plot[,1] = X_prior[which(X_prior[,2]==1),4]
# disagree_plot[,2] = Y_prior[which(X_prior[,2]==1)]
# 
# par(mar = c(5, 4, 1, 2))
# plot(agree_plot[,1],agree_plot[,2],xlim = c(0,1), ylim = c(0,1))
# 
# points(disagree_plot[,1],disagree_plot[,2],col="red")
# 
# abline(a=0,b=1,lty = 3)
# 
# plot_index <- floor(runif(4000, min=1, max=12001))
# 
# q = seq(0,1,by = .01)
# 
# priorBeta1_reg <- rep(NaN,length(q))
# priorBeta2_reg <- rep(NaN,length(q))
# 
# for(i in 1:nrow(priorBeta_samples)){
# 
#   for(j in 1:length(q)){
#     priorBeta1_reg[j] <- priorBeta_samples[i,1] + priorBeta_samples[i,3]*q[j]
#   }
# 
#   for(j in 1:length(q)){
#     priorBeta2_reg[j] <- priorBeta_samples[i,2] + priorBeta_samples[i,4]*q[j]
#   }
# 
#   lines(q,priorBeta1_reg,col=c(rgb(0,0,1,0.005)))
#   lines(q,priorBeta2_reg,col=c(rgb(1,0,0,0.005)))
# }

# Latent mixture version 1 (2 disagreers) ---------------------------------

nSamples = 4000

nBeta = 2

X_agree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)
X_disagree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)

X_agree[,1] <- X_prior[,1]
X_agree[,2] <- X_prior[,3]
X_disagree[,1] <- X_prior[,2]
X_disagree[,2] <- X_prior[,4]

Y<- Y_prior
nParticipants <- max(Participant_list)
nTrials <- length(which(Participant_prior==1))

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)
IDmat <- matrix(0,nrow = nBeta, ncol = nBeta)
IDmat0 <- matrix(0,ncol = nBeta, nrow = 1)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
  IDmat[i,i] <- 1
}

data <- list("X_agree","X_disagree","Y","bPrec","nBeta",
             "nParticipants","nTrials") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA_agree","BETA_disagree",
                "tau","sigma","alpha","Ypred",
                "phi_disagree","psi_disagree","psi_disagree_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA_agree" = runif(nBeta,-.5,.5),
                     "BETA_disagree" = runif(nBeta,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_latentMix.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

latentMix1_Beta_agree_samples <- samples$BUGSoutput$sims.list$BETA_agree

latentMix1_Beta_disagree_samples <- samples$BUGSoutput$sims.list$BETA_disagree

latentMix1_phi_disagree_samples <- samples$BUGSoutput$sims.list$phi_disagree

latentMix1_psi_disagree_samples <- samples$BUGSoutput$sims.list$psi_disagree

latentMix1_psi_disagree_bern_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

latentMix1_sigma_samples <- samples$BUGSoutput$sims.list$sigma

latentMix1_alpha_samples <- samples$BUGSoutput$sims.list$alpha

latentMix1_tau_samples <- samples$BUGSoutput$sims.list$tau

latentMix1_summary <- samples$BUGSoutput$summary

latentMix1_DIC <- samples$BUGSoutput$DIC

# Latent mixture version 2 (2 agreers) ---------------------------------

nSamples = 4000

nBeta = 2

X_agree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)
X_disagree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)

X_agree[,1] <- X_prior[,1]
X_agree[,2] <- X_prior[,3]
X_disagree[,1] <- X_prior[,2]
X_disagree[,2] <- X_prior[,4]

Y<- Y_prior
nParticipants <- max(Participant_list)
nTrials <- length(which(Participant_prior==1))

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)
IDmat <- matrix(0,nrow = nBeta, ncol = nBeta)
IDmat0 <- matrix(0,ncol = nBeta, nrow = 1)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
  IDmat[i,i] <- 1
}

data <- list("X_agree","X_disagree","Y","bPrec","nBeta",
             "nParticipants","nTrials") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA_agree1","BETA_agree2","BETA_disagree",
                "tau","sigma","alpha","Ypred",
                "phi_agree","psi_agree","psi_agree_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA_agree1_temp" = runif(nBeta,-.5,.5),
                     "BETA_agree2_temp" = runif(nBeta,-.5,.5),
                     "BETA_disagree_temp" = runif(nBeta,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_prior_to_int_latentMix2.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

latentMix2_Beta_agree1_samples <- samples$BUGSoutput$sims.list$BETA_agree1
latentMix2_Beta_agree2_samples <- samples$BUGSoutput$sims.list$BETA_agree2

latentMix2_Beta_disagree_samples <- samples$BUGSoutput$sims.list$BETA_disagree

latentMix2_phi_agree_samples <- samples$BUGSoutput$sims.list$phi_agree

latentMix2_psi_agree_samples <- samples$BUGSoutput$sims.list$psi_agree

latentMix2_sigma_samples <- samples$BUGSoutput$sims.list$sigma

latentMix2_alpha_samples <- samples$BUGSoutput$sims.list$alpha

latentMix2_tau_samples <- samples$BUGSoutput$sims.list$tau

latentMix2_summary <- samples$BUGSoutput$summary

latentMix2_DIC <- samples$BUGSoutput$DIC

# Latent mixture version 3 (2 agreers, 2 disagreers; 4 combos) ------------

nSamples = 4000

nBeta = 2

X_agree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)
X_disagree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)

X_agree[,1] <- X_prior[,1]
X_agree[,2] <- X_prior[,3]
X_disagree[,1] <- X_prior[,2]
X_disagree[,2] <- X_prior[,4]

Y<- Y_prior
nParticipants <- max(Participant_list)
nTrials <- length(which(Participant_prior==1))

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

data <- list("X_agree","X_disagree","Y","bPrec","nBeta","nParticipants","nTrials") # to be passed on to JAGS
# parameters to be monitored:
parameters <- c("BETA_agree1","BETA_agree2","BETA_disagree1","BETA_disagree2",
                "tau","sigma","alpha","Ypred",
                "phi_agree","psi_agree","psi_agree_bern","sigma_psi_agree",
                "phi_disagree","psi_disagree","psi_disagree_bern","sigma_psi_disagree")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA_agree1_temp" = runif(nBeta,-.5,.5),
                     "BETA_agree2_temp" = runif(nBeta,-.5,.5),
                     "BETA_disagree1_temp" = runif(nBeta,-.5,.5),
                     "BETA_disagree2_temp" = runif(nBeta,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits,
                parameters.to.save=parameters,
                model.file="LoLComp_prior_to_int_latentMix3.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object,
# ready for inspection.

latentMix3_Beta_agree1_samples <- samples$BUGSoutput$sims.list$BETA_agree1
latentMix3_Beta_agree2_samples <- samples$BUGSoutput$sims.list$BETA_agree2

latentMix3_Beta_disagree1_samples <- samples$BUGSoutput$sims.list$BETA_disagree1
latentMix3_Beta_disagree2_samples <- samples$BUGSoutput$sims.list$BETA_disagree2

latentMix3_phi_agree_samples <- samples$BUGSoutput$sims.list$phi_agree
latentMix3_phi_disagree_samples <- samples$BUGSoutput$sims.list$phi_disagree

latentMix3_psi_agree_samples <- samples$BUGSoutput$sims.list$psi_agree
latentMix3_psi_disagree_samples <- samples$BUGSoutput$sims.list$psi_disagree

latentMix3_psi_agree_bern_samples <- samples$BUGSoutput$sims.list$psi_agree_bern
latentMix3_psi_disagree_bern_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

latentMix3_sigma_samples <- samples$BUGSoutput$sims.list$sigma

latentMix3_alpha_samples <- samples$BUGSoutput$sims.list$alpha

latentMix3_tau_samples <- samples$BUGSoutput$sims.list$tau

latentMix3_summary <- samples$BUGSoutput$summary

latentMix3_DIC <- samples$BUGSoutput$DIC

# Latent mixture version 4 (2 respoders) -------------------------

nSamples = 4000

nBeta = 2

X_agree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)
X_disagree <- matrix(NaN, nrow = length(Y_prior),ncol = 2)

X_agree[,1] <- X_prior[,1]
X_agree[,2] <- X_prior[,3]
X_disagree[,1] <- X_prior[,2]
X_disagree[,2] <- X_prior[,4]

Y<- Y_prior
nParticipants <- max(Participant_list)
nTrials <- length(which(Participant_prior==1))

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

data <- list("X_agree","X_disagree","Y","bPrec","nBeta","nParticipants","nTrials") # to be passed on to JAGS
# parameters to be monitored:
parameters <- c("BETA_agree1","BETA_agree2","BETA_disagree1","BETA_disagree2",
                "tau","sigma","alpha","Ypred","phi","psi","psi_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA_agree1" = runif(nBeta,-.5,.5),
                     "BETA_agree2" = runif(nBeta,-.5,.5),
                     "BETA_disagree1" = runif(nBeta,-.5,.5),
                     "BETA_disagree2" = runif(nBeta,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits,
                parameters.to.save=parameters,
                model.file="LoLComp_prior_to_int_latentMix4.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object,
# ready for inspection.

latentMix4_Beta_agree1_samples <- samples$BUGSoutput$sims.list$BETA_agree1
latentMix4_Beta_agree2_samples <- samples$BUGSoutput$sims.list$BETA_agree2

latentMix4_Beta_disagree1_samples <- samples$BUGSoutput$sims.list$BETA_disagree1
latentMix4_Beta_disagree2_samples <- samples$BUGSoutput$sims.list$BETA_disagree2

latentMix4_phi_samples <- samples$BUGSoutput$sims.list$phi

latentMix4_psi_samples <- samples$BUGSoutput$sims.list$psi

latentMix4_sigma_samples <- samples$BUGSoutput$sims.list$sigma

latentMix4_alpha_samples <- samples$BUGSoutput$sims.list$alpha

latentMix4_tau_samples <- samples$BUGSoutput$sims.list$tau

latentMix4_summary <- samples$BUGSoutput$summary

latentMix4_DIC <- samples$BUGSoutput$DIC

# Graph latent mixture model 1 (2 disagreers) -----------------------------------

par(mfrow=c(1,1))

X_prior_graph_latentMix <- matrix(NaN,nrow = nrow(X_prior),ncol = 5)
X_prior_graph_latentMix[,1:4] <- X_prior

for(i in 1:nParticipants){
  hi = i*nTrials
  lo = hi-(nTrials-1)
  X_prior_graph_latentMix[lo:hi,5] <- rep(i,nTrials)
}

Participant_psi <- matrix(nrow = max(Participant_list),ncol = 2)
Participant_psi[,1] <- 1:max(Participant_list)
for(i in 1:max(Participant_list)){
  Participant_psi[i,2] <- mean(latentMix1_psi_disagree_samples[,i])
}

agree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,1]==1)),ncol = 2)

agree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,1]==1),3]
agree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,1]==1)]

disagree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,2]==1)),ncol = 2)

disagree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,2]==1),4]
disagree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,2]==1)]

par(mar = c(5, 4, 1, 2))
plot(agree_plot[,1],agree_plot[,2],xlim = c(0,1), ylim = c(0,1))

points(disagree_plot[,1],disagree_plot[,2],col="red")

abline(a=0,b=1,lty = 3)

plot_index <- floor(runif(4000, min=1, max=12001))

q = seq(0,1,by = .01)

latentMix1_Beta_agree_reg <- rep(NaN,length(q))
latentMix1_Beta_disagree1_reg <- rep(NaN,length(q))
latentMix1_Beta_disagree2_reg <- rep(NaN,length(q))

for(i in 1:nrow(latentMix1_Beta_agree_samples)){

  for(j in 1:length(q)){
    latentMix1_Beta_agree_reg[j] <- latentMix1_Beta_agree_samples[i,1] + latentMix1_Beta_agree_samples[i,2]*q[j]
  }

  for(j in 1:length(q)){
    latentMix1_Beta_disagree1_reg[j] <- latentMix1_Beta_disagree1_samples[i,1] + latentMix1_Beta_disagree1_samples[i,2]*q[j]
  }

  for(j in 1:length(q)){
    latentMix1_Beta_disagree2_reg[j] <- latentMix1_Beta_disagree2_samples[i,1] + latentMix1_Beta_disagree2_samples[i,2]*q[j]
  }

  lines(q,latentMix1_Beta_agree_reg,col=c(rgb(0,0,1,0.005)))
  lines(q,latentMix1_Beta_disagree1_reg,col=c(rgb(.2,0,0,0.005)))
  lines(q,latentMix1_Beta_disagree2_reg,col=c(rgb(.8,.2,0,0.005)))

}

# # Graph latent mixture model 2 (2 agreers) -----------------------------------
# 
# X_prior_graph_latentMix <- matrix(NaN,nrow = nrow(X_prior),ncol = 5)
# X_prior_graph_latentMix[,1:4] <- X_prior
# 
# for(i in 1:nParticipants){
#   hi = i*nTrials
#   lo = hi-(nTrials-1)
#   X_prior_graph_latentMix[lo:hi,5] <- rep(i,nTrials)
# }
# 
# Participant_psi <- matrix(nrow = max(Participant_list),ncol = 2)
# Participant_psi[,1] <- 1:max(Participant_list)
# for(i in 1:max(Participant_list)){
#   Participant_psi[i,2] <- mean(latentMix2_psi_agree_samples[,i])
# }
# 
# agree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,1]==1)),ncol = 2)
# 
# agree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,1]==1),3]
# agree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,1]==1)]
# 
# disagree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,2]==1)),ncol = 2)
# 
# disagree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,2]==1),4]
# disagree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,2]==1)]
# 
# par(mar = c(5, 4, 1, 2))
# plot(agree_plot[,1],agree_plot[,2],xlim = c(0,1), ylim = c(0,1))
# 
# points(disagree_plot[,1],disagree_plot[,2],col="red")
# 
# abline(a=0,b=1,lty = 3)
# 
# plot_index <- floor(runif(4000, min=1, max=12001))
# 
# q = seq(0,1,by = .01)
# 
# latentMix2_Beta_agree1_reg <- rep(NaN,length(q))
# latentMix2_Beta_agree2_reg <- rep(NaN,length(q))
# latentMix2_Beta_disagree_reg <- rep(NaN,length(q))
# 
# for(i in 1:nrow(latentMix2_Beta_agree1_samples)){
#   
#   for(j in 1:length(q)){
#     latentMix2_Beta_agree1_reg[j] <- latentMix2_Beta_agree1_samples[i,1] + latentMix2_Beta_agree1_samples[i,2]*q[j]
#   }
#   
#   for(j in 1:length(q)){
#     latentMix2_Beta_agree2_reg[j] <- latentMix2_Beta_agree2_samples[i,1] + latentMix2_Beta_agree2_samples[i,2]*q[j]
#   }
#   
#   for(j in 1:length(q)){
#     latentMix2_Beta_disagree_reg[j] <- latentMix2_Beta_disagree_samples[i,1] + latentMix2_Beta_disagree_samples[i,2]*q[j]
#   }
#   
#   lines(q,latentMix2_Beta_agree1_reg,col=c(rgb(0,0,.2,0.005)))
#   lines(q,latentMix2_Beta_agree2_reg,col=c(rgb(0,.2,.8,0.005)))
#   lines(q,latentMix2_Beta_disagree_reg,col=c(rgb(1,0,0,0.005)))
#   
# }

# Graph latent mixture model 4 (2 agreers, 2 disagreers; 2 types of participant) -----------------------------------
X_prior_graph_latentMix <- matrix(NaN,nrow = nrow(X_prior),ncol = 5)
X_prior_graph_latentMix[,1:4] <- X_prior

agree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,1]==1)),ncol = 2)

agree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,1]==1),3]
agree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,1]==1)]

disagree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,2]==1)),ncol = 2)

disagree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,2]==1),4]
disagree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,2]==1)]

par(mar = c(5, 4, 1, 2))
plot(agree_plot[,1],agree_plot[,2],xlim = c(0,1), ylim = c(0,1))

points(disagree_plot[,1],disagree_plot[,2],col="red")

abline(a=0,b=1,lty = 3)

plot_index <- floor(runif(4000, min=1, max=12001))

q = seq(0,1,by = .01)

latentMix4_Beta_agree1_reg <- rep(NaN,length(q))
latentMix4_Beta_agree2_reg <- rep(NaN,length(q))
latentMix4_Beta_disagree1_reg <- rep(NaN,length(q))
latentMix4_Beta_disagree2_reg <- rep(NaN,length(q))

for(i in 1:nrow(latentMix4_Beta_agree1_samples)){
  
  for(j in 1:length(q)){
    latentMix4_Beta_agree1_reg[j] <- latentMix4_Beta_agree1_samples[i,1] + latentMix4_Beta_agree1_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix4_Beta_agree2_reg[j] <- latentMix4_Beta_agree2_samples[i,1] + latentMix4_Beta_agree2_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix4_Beta_disagree1_reg[j] <- latentMix4_Beta_disagree1_samples[i,1] + latentMix4_Beta_disagree1_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix4_Beta_disagree2_reg[j] <- latentMix4_Beta_disagree2_samples[i,1] + latentMix4_Beta_disagree2_samples[i,2]*q[j]
  }
  
  lines(q,latentMix4_Beta_agree1_reg,col=c(rgb(0,0,.2,0.005)))
  lines(q,latentMix4_Beta_agree2_reg,col=c(rgb(0,0,1,0.005)))
  lines(q,latentMix4_Beta_disagree1_reg,col=c(rgb(.2,0,0,0.005)))
  lines(q,latentMix4_Beta_disagree2_reg,col=c(rgb(1,0,0,0.005)))
  
}

# Graph latent mixture model 3 (2 agreers, 2 disagreers; 4 combos) -----------------------------------

X_prior_graph_latentMix <- matrix(NaN,nrow = nrow(X_prior),ncol = 5)
X_prior_graph_latentMix[,1:4] <- X_prior

for(i in 1:nParticipants){
  hi = i*nTrials
  lo = hi-(nTrials-1)
  X_prior_graph_latentMix[lo:hi,5] <- rep(i,nTrials)
}

Participant_psi <- matrix(nrow = max(Participant_list),ncol = 2)
Participant_psi[,1] <- 1:max(Participant_list)
for(i in 1:max(Participant_list)){
  Participant_psi[i,2] <- mean(latentMix3_psi_disagree_samples[,i])
}

agree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,1]==1)),ncol = 2)

agree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,1]==1),3]
agree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,1]==1)]

disagree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,2]==1)),ncol = 2)

disagree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,2]==1),4]
disagree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,2]==1)]

par(mar = c(5, 4, 1, 2))
plot(agree_plot[,1],agree_plot[,2],xlim = c(0,1), ylim = c(0,1))

points(disagree_plot[,1],disagree_plot[,2],col="red")

abline(a=0,b=1,lty = 3)

plot_index <- floor(runif(4000, min=1, max=12001))

q = seq(0,1,by = .01)

latentMix3_Beta_agree1_reg <- rep(NaN,length(q))
latentMix3_Beta_agree2_reg <- rep(NaN,length(q))
latentMix3_Beta_disagree1_reg <- rep(NaN,length(q))
latentMix3_Beta_disagree2_reg <- rep(NaN,length(q))

for(i in 1:nrow(latentMix3_Beta_agree1_samples)){
  
  for(j in 1:length(q)){
    latentMix3_Beta_agree1_reg[j] <- latentMix3_Beta_agree1_samples[i,1] + latentMix3_Beta_agree1_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix3_Beta_agree2_reg[j] <- latentMix3_Beta_agree2_samples[i,1] + latentMix3_Beta_agree2_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix3_Beta_disagree1_reg[j] <- latentMix3_Beta_disagree1_samples[i,1] + latentMix3_Beta_disagree1_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix3_Beta_disagree2_reg[j] <- latentMix3_Beta_disagree2_samples[i,1] + latentMix3_Beta_disagree2_samples[i,2]*q[j]
  }
  
  lines(q,latentMix3_Beta_agree1_reg,col=c(rgb(0,0,.2,0.005)))
  lines(q,latentMix3_Beta_agree2_reg,col=c(rgb(0,0,1,0.005)))
  lines(q,latentMix3_Beta_disagree1_reg,col=c(rgb(.2,0,0,0.005)))
  lines(q,latentMix3_Beta_disagree2_reg,col=c(rgb(1,0,0,0.005)))
  
}

# Graph latent mixture model 1, split agree vs disagree -----------------------------------
par(mfrow=c(2,1))

X_prior_graph_latentMix <- matrix(NaN,nrow = nrow(X_prior),ncol = 5)
X_prior_graph_latentMix[,1:4] <- X_prior

for(i in 1:nParticipants){
  hi = i*nTrials
  lo = hi-(nTrials-1)
  X_prior_graph_latentMix[lo:hi,5] <- rep(i,nTrials)
}

Participant_psi <- matrix(nrow = max(Participant_list),ncol = 2)
Participant_psi[,1] <- 1:max(Participant_list)
for(i in 1:max(Participant_list)){
  Participant_psi[i,2] <- mean(latentMix3_psi_disagree_samples[,i])
}

agree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,1]==1)),ncol = 2)

agree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,1]==1),3]
agree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,1]==1)]

disagree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,2]==1)),ncol = 2)

disagree_plot[,1] = X_prior_graph_latentMix[which(X_prior_graph_latentMix[,2]==1),4]
disagree_plot[,2] = Y_prior[which(X_prior_graph_latentMix[,2]==1)]

par(mar = c(5, 4, 1, 2))
plot(agree_plot[,1],agree_plot[,2],xlim = c(0,1), ylim = c(0,1),
     xlab = "prior",ylab="interpretation")


abline(a=0,b=1,lty = 3)

plot_index <- floor(runif(4000, min=1, max=12001))

q = seq(0,1,by = .01)

latentMix3_Beta_agree1_reg <- rep(NaN,length(q))
latentMix3_Beta_agree2_reg <- rep(NaN,length(q))
latentMix3_Beta_disagree1_reg <- rep(NaN,length(q))
latentMix3_Beta_disagree2_reg <- rep(NaN,length(q))

for(i in 1:nrow(latentMix3_Beta_agree1_samples)){
  
  for(j in 1:length(q)){
    latentMix3_Beta_agree1_reg[j] <- latentMix3_Beta_agree1_samples[i,1] + latentMix3_Beta_agree1_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix3_Beta_agree2_reg[j] <- latentMix3_Beta_agree2_samples[i,1] + latentMix3_Beta_agree2_samples[i,2]*q[j]
  }
  
  lines(q,latentMix3_Beta_agree1_reg,col=c(rgb(0,0,.2,0.005)))
  lines(q,latentMix3_Beta_agree2_reg,col=c(rgb(0,0,1,0.005)))
}


plot(disagree_plot[,1],disagree_plot[,2],col="red",xlim = c(0,1), ylim = c(0,1),
     xlab = "prior",ylab="interpretation")

for(i in 1:nrow(latentMix3_Beta_disagree1_samples)){
  
  for(j in 1:length(q)){
    latentMix3_Beta_disagree1_reg[j] <- latentMix3_Beta_disagree1_samples[i,1] + latentMix3_Beta_disagree1_samples[i,2]*q[j]
  }
  
  for(j in 1:length(q)){
    latentMix3_Beta_disagree2_reg[j] <- latentMix3_Beta_disagree2_samples[i,1] + latentMix3_Beta_disagree2_samples[i,2]*q[j]
  }
  
  lines(q,latentMix3_Beta_disagree1_reg,col=c(rgb(.2,0,0,0.005)))
  lines(q,latentMix3_Beta_disagree2_reg,col=c(rgb(1,0,0,0.005)))
  
}

# Agree plots by individual -----------------------------------------------
par(mfrow=c(5,4))

ind_latentMix <- cbind(X_prior_graph_latentMix,Y_prior)

agree_ind <- ind_latentMix[which(ind_latentMix[,1]==1),]

for(i in 1:nParticipants){
  
  plot(agree_ind[which(agree_ind[,5]==i),3],agree_ind[which(agree_ind[,5]==i),6],
       main = i, xlab = 'Prior', ylab = 'Interpretation',xlim = c(0,1),ylim = c(0,1))
}

# Disagree plots by individual -----------------------------------------------
par(mfrow=c(5,4))

disagree_ind <- ind_latentMix[which(ind_latentMix[,2]==1),]

for(i in 1:nParticipants){
  
  plot(disagree_ind[which(disagree_ind[,5]==i),4],disagree_ind[which(disagree_ind[,5]==i),6],
       main = i, xlab = 'Prior', ylab = 'Interpretation',xlim = c(0,1),ylim = c(0,1),col = 'red')
  abline(a=0,b=1)
}

# Participant-level tendencies --------------------------------------------
par(mfrow=c(1,1))
Participant_psi <- matrix(nrow = max(Participant_list),ncol = 2)
Participant_psi[,1] <- 1:max(Participant_list)
for(i in 1:max(Participant_list)){
  Participant_psi[i,2] <- mean(latentMix1_psi_disagree_samples[,i])
}

hist(Participant_psi[,2])

plot(Participant_psi[,1],Participant_psi[,2],type = 'n')

for(i in 1:nrow(Participant_psi)){
  text(Participant_psi[i,1],Participant_psi[i,2],label=Participant_psi[i,1])
}

# Which data with which beta ----------------------------------------------

psi_disagree_bern_plot <- rep(NaN,length(Y_prior))
for(i in 1:length(Y_prior)){
  psi_disagree_bern_plot[i] <-  mean(latentMix1_psi_disagree_bern_samples[,i])
}

par(mfrow=c(1,1))

hist(psi_disagree_bern_plot)


X_prior_good_latentMix <- cbind(X_prior_graph_latentMix,Y_prior,
                                psi_disagree_bern_plot)

#X_prior_good_latentMix[which(X_prior_good_latentMix[,1]==1),7] <- 0

hist(X_prior_good_latentMix[which(X_prior_good_latentMix[,1]==0),7])


# Make csv to check impact on main result ---------------------------------

OnlineData$Early_psi_bern <- rep(0, nrow(OnlineData))

OnlineData$Late_psi_bern <- rep(0, nrow(OnlineData))
ticker=1

for(i in 1:nrow(OnlineData)){
  if(OnlineData$SpeakerEarly[i]==1){
    OnlineData$Early_psi_bern[i] <- psi_disagree_bern_plot[ticker]
    ticker = ticker+1
  }
  if(OnlineData$SpeakerLate[i]==1){
    OnlineData$Late_psi_bern[i] <- psi_disagree_bern_plot[ticker]
    ticker = ticker+1
  }
  
  if(OnlineData$EarlyTVJ[i]==1){
    OnlineData$Early_psi_bern[i] <- 0
  }
  
  if(OnlineData$LateTVJ[i]==1){
    OnlineData$Late_psi_bern[i] <- 0
  }
}

write.csv(OnlineData, "OnlineData_psi.csv")

