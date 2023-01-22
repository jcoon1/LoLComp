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

# Latent mixture version 1 (2 disagreers) ---------------------------------

nSamples = 4000

nBeta = 2
Participant_list <- Participant_prior
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

# Latent mixture version 2 (split for loops) ---------------------------------

nSamples = 4000

X_prior_par <- cbind(X_prior,Participant_list,Y_prior)

X_agree_temp <- X_prior_par[which(X_prior_par[,1]==1),]
X_disagree_temp <- X_prior_par[which(X_prior_par[,2]==1),]

nBeta = 2
X_agree <- matrix(NaN, nrow = nrow(X_agree_temp),ncol = 2)
X_disagree <- matrix(NaN, nrow = nrow(X_disagree_temp),ncol = 2)

X_agree[,1] <- X_agree_temp[,1]
X_agree[,2] <- X_agree_temp[,3]
Participant_list_agree <- X_agree_temp[,5]
Y_agree <- X_agree_temp[,6]

X_disagree[,1] <- X_disagree_temp[,2]
X_disagree[,2] <- X_disagree_temp[,4]
Participant_list_disagree <- X_disagree_temp[,5]
Y_disagree <- X_disagree_temp[,6]

nParticipants <- max(Participant_list_agree,Participant_list_disagree)
nAgree <- length(Y_agree)
nDisagree <- length(Y_disagree)
nTrials <- length(which(Participant_prior==1))

bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)
IDmat <- matrix(0,nrow = nBeta, ncol = nBeta)
IDmat0 <- matrix(0,ncol = nBeta, nrow = 1)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
  IDmat[i,i] <- 1
}

data <- list("X_agree","X_disagree","Y_agree","Y_disagree",
             "Participant_list_agree","Participant_list_disagree",
             "bPrec","nBeta","nParticipants","nAgree","nDisagree") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA_agree","BETA_disagree",
                "tau","sigma","alpha",
                "phi_disagree","psi_disagree","psi_disagree_bern","sigma_psi")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA_agree" = runif(nBeta,-.5,.5),
                     "BETA_disagree" = runif(nBeta,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_latentMix_splitfor.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

latentMix2_Beta_agree_samples <- samples$BUGSoutput$sims.list$BETA_agree

latentMix2_Beta_disagree_samples <- samples$BUGSoutput$sims.list$BETA_disagree

latentMix2_phi_disagree_samples <- samples$BUGSoutput$sims.list$phi_disagree

latentMix2_psi_disagree_samples <- samples$BUGSoutput$sims.list$psi_disagree

latentMix2_psi_disagree_bern_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

latentMix2_sigma_samples <- samples$BUGSoutput$sims.list$sigma

latentMix2_alpha_samples <- samples$BUGSoutput$sims.list$alpha

latentMix2_tau_samples <- samples$BUGSoutput$sims.list$tau

latentMix2_summary <- samples$BUGSoutput$summary

latentMix2_DIC <- samples$BUGSoutput$DIC

# Graph latent mixture model (2 disagreers) -----------------------------------

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
latentMix1_Beta_disagree_reg <- rep(NaN,length(q))

for(i in 1:nrow(latentMix1_Beta_agree_samples)){

  for(j in 1:length(q)){
    latentMix1_Beta_agree_reg[j] <- latentMix1_Beta_agree_samples[i,1] + latentMix1_Beta_agree_samples[i,2]*q[j]
  }

  for(j in 1:length(q)){
    latentMix1_Beta_disagree_reg[j] <- latentMix1_Beta_disagree_samples[i,1] + latentMix1_Beta_disagree_samples[i,2]*q[j]
  }

  lines(q,latentMix1_Beta_agree_reg,col=c(rgb(0,0,1,0.005)))
  lines(q,latentMix1_Beta_disagree_reg,col=c(rgb(.2,0,0,0.005)))

}

# Graph latent mixture model, split agree vs disagree -----------------------------------
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
plot(agree_plot[,1],agree_plot[,2],xlim = c(0,1), ylim = c(0,1),
     xlab = "prior",ylab="interpretation")


plot_index <- floor(runif(4000, min=1, max=12001))

q = seq(0,1,by = .01)

latentMix1_Beta_agree_reg <- rep(NaN,length(q))
latentMix1_Beta_disagree_reg <- rep(NaN,length(q))

for(i in 1:nrow(latentMix1_Beta_agree_samples)){
  
  for(j in 1:length(q)){
    latentMix1_Beta_agree_reg[j] <- latentMix1_Beta_agree_samples[i,1] + latentMix1_Beta_agree_samples[i,2]*q[j]
  }
  
  lines(q,latentMix1_Beta_agree_reg,col=c(rgb(0,0,1,0.005)))
}


plot(disagree_plot[,1],disagree_plot[,2],col="red",xlim = c(0,1), ylim = c(0,1),
     xlab = "prior",ylab="interpretation")

abline(a=0,b=1,lty = 3)

for(i in 1:nrow(latentMix1_Beta_disagree_samples)){
  
  for(j in 1:length(q)){
    latentMix1_Beta_disagree_reg[j] <- latentMix1_Beta_disagree_samples[i,1] + latentMix1_Beta_disagree_samples[i,2]*q[j]
  }
  
  lines(q,latentMix1_Beta_disagree_reg,col=c(rgb(.2,0,0,0.005)))
  
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


# Make a fancy plot where color indicates n.eff ---------------------------

X_prior_graph_latentMix_neff <- cbind(X_prior_graph_latentMix,
                                      latentMix1_summary[693:1280,9])

disagree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,2]==1)),ncol = 3)

disagree_plot[,1] = X_prior_graph_latentMix_neff[which(X_prior_graph_latentMix_neff[,2]==1),4]
disagree_plot[,2] = Y_prior[which(X_prior_graph_latentMix_neff[,2]==1)]
disagree_plot[,3] = X_prior_graph_latentMix_neff[which(X_prior_graph_latentMix_neff[,2]==1),6]


plot(disagree_plot[,1],disagree_plot[,2],col="red",xlim = c(0,1), ylim = c(0,1),
     xlab = "prior",ylab="interpretation")

# Make a fancy plot where color indicates psi_bern ---------------------------
X_prior_graph_latentMix_psi <- cbind(X_prior_graph_latentMix,
                                      latentMix1_summary[693:1280,1])

disagree_plot = matrix(NaN,nrow = length(which(X_prior_graph_latentMix[,2]==1)),ncol = 3)

disagree_plot[,1] = X_prior_graph_latentMix_psi[which(X_prior_graph_latentMix_psi[,2]==1),4]
disagree_plot[,2] = Y_prior[which(X_prior_graph_latentMix_psi[,2]==1)]
disagree_plot[,3] = X_prior_graph_latentMix_psi[which(X_prior_graph_latentMix_psi[,2]==1),6]


plot(disagree_plot[,1],disagree_plot[,2],col=rgb(disagree_plot[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation")

abline(a=0,b=1,lty = 3)

for(i in 1:nrow(latentMix1_Beta_disagree_samples)){
  
  for(j in 1:length(q)){
    latentMix1_Beta_disagree_reg[j] <- latentMix1_Beta_disagree_samples[i,1] + latentMix1_Beta_disagree_samples[i,2]*q[j]
  }
  
  lines(q,latentMix1_Beta_disagree_reg,col=c(rgb(.2,0,0,0.005)))
  
}

# Compare psi_bern for 1 and 2 --------------------------------------------
X_prior_graph_latentMix_psi2 <- 
  cbind(X_prior_graph_latentMix_psi,rep(NaN,nrow(X_prior_graph_latentMix_psi)),rep(NaN,nrow(X_prior_graph_latentMix_psi)))

psi_bern_disagree2 <- rep(NaN,ncol(latentMix2_psi_disagree_bern_samples))
for(i in 1:ncol(latentMix2_psi_disagree_bern_samples)){
  psi_bern_disagree2[i] <- mean(latentMix2_psi_disagree_bern_samples[,i])
}
ticker=1
for (i in 1:nrow(X_prior_graph_latentMix_psi)){
  if(X_prior_graph_latentMix_psi2[i,2]==1){
    X_prior_graph_latentMix_psi2[i,7] <- psi_bern_disagree2[ticker]
    X_prior_graph_latentMix_psi2[i,8] <- X_prior_graph_latentMix_psi2[i,6]-X_prior_graph_latentMix_psi2[i,7]
    ticker = ticker+1
  }
}

disagree_plot2 = matrix(NaN,nrow = length(which(X_prior_graph_latentMix_psi2[,2]==1)),ncol = 3)

disagree_plot2[,1] = X_prior_graph_latentMix_psi[which(X_prior_graph_latentMix_psi2[,2]==1),4]
disagree_plot2[,2] = Y_prior[which(X_prior_graph_latentMix_psi2[,2]==1)]
disagree_plot2[,3] = X_prior_graph_latentMix_psi2[which(X_prior_graph_latentMix_psi2[,2]==1),7]


plot(disagree_plot2[,1],disagree_plot2[,2],col=rgb(disagree_plot2[,3],0,0), 
     xlim = c(0,1), ylim = c(0,1), xlab = "prior",ylab="interpretation")

abline(a=0,b=1,lty = 3)

latentMix2_Beta_disagree_reg <- rep(NaN,length(q))

for(i in 1:nrow(latentMix2_Beta_disagree_samples)){
  
  for(j in 1:length(q)){
    latentMix2_Beta_disagree_reg[j] <- latentMix2_Beta_disagree_samples[i,1] + latentMix2_Beta_disagree_samples[i,2]*q[j]
  }
  
  lines(q,latentMix2_Beta_disagree_reg,col=c(rgb(.2,0,0,0.005)))
  
}