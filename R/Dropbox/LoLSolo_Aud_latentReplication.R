# clears workspace:  
rm(list=ls()) 

setwd('C:/Users/Jeff/Documents/R/upload/LoLSolo_Aud')

load("latentReplication.RData")

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)
library(fMultivar)

source("BootstrapMean95.R")


# Get experienced participants' data --------------------------------------

Experienced_replication <- read_csv("Experienced_replication_csv.csv",col_names=TRUE)

nReps=12

Participant <- seq(1:nrow(Experienced_replication))

Experienced_replication <- cbind(Experienced_replication,Participant)

# Get inexperienced participants' data ------------------------------------

Inexperienced_replication <- read_csv("Inexperienced_replication_csv.csv",col_names=TRUE)

nParticipants <- nrow(Inexperienced_replication)

Participant <- 1:nParticipants

Inexperienced_replication <- cbind(Inexperienced_replication,Participant)

nResponses_SONA <- nParticipants * nReps

nParticipants <- nrow(Inexperienced_replication)

nResponses_SONA <- nParticipants * nReps

SONA_vector <- rep(NaN,nResponses_SONA)

for(i in 1:nParticipants){
  hi = i*nReps
  lo = hi-(nReps-1)
  hi_grab = 1 + 49
  lo_grab = nReps + 49
  SONA_vector[lo:hi] <- Inexperienced_replication[i,lo_grab:hi_grab]
}

SONA_vector <- unlist(SONA_vector)

SONA_descriptives <- matrix(NaN,nrow = 3,ncol = 1)

SONA_descriptives[1,1] <- mean(SONA_vector)

nIterations = 1000
nSamples = 100

SONA_descriptives[2:3,1] <- BootstrapMean95(SONA_vector,nSamples,nIterations)

### Inexperienced_replication is main inexperienced participants' data source moving forward


# Expected model (only experienced varies by condition): prep data ----------------

iterations = 50

DIC_standard <- rep(NaN,50)

for(k in 1:iterations){

Agree_index <- matrix(NaN,nrow = nrow(Experienced_replication),ncol = 12)

for(i in 1:nrow(Agree_index)){
  for(j in 1:ncol(Agree_index)){
    TVJ_ticker = j + 49
    Agree <- Experienced_replication[i,TVJ_ticker] == 1
    Agree_index[i,j] <- Agree
  }
}

Agree_data <- matrix(NaN,nrow = length(which(Agree_index==1)),ncol = 4)
Disagree_data <- matrix(NaN,nrow = length(which(Agree_index==0)),ncol = 4)

Agree_ticker = 1
Disagree_ticker = 1

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:ncol(Agree_index)){
    typeTicker = j + 110
    prior_ticker = j+65
    interpret_ticker = j+77
    if(Agree_index[i,j]==1){
      Agree_data[Agree_ticker,1] <- unlist(Experienced_replication[i,prior_ticker])
      Agree_data[Agree_ticker,2] <- unlist(Experienced_replication[i,interpret_ticker])
      Agree_data[Agree_ticker,3] <- i
      Agree_data[Agree_ticker,4] <- unlist(Experienced_replication[i,typeTicker])+3
      Agree_ticker = Agree_ticker+1
    }
    if(Agree_index[i,j]==0){
      Disagree_data[Disagree_ticker,1] <- unlist(Experienced_replication[i,prior_ticker])
      Disagree_data[Disagree_ticker,2] <- unlist(Experienced_replication[i,interpret_ticker])
      Disagree_data[Disagree_ticker,3] <- i
      Disagree_data[Disagree_ticker,4] <- unlist(Experienced_replication[i,typeTicker])+3
      Disagree_ticker = Disagree_ticker+1
    }
  }
}

nExperts = nrow(Experienced_replication)
nNovices = nrow(Inexperienced_replication)

nBeta = 4
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}


## Set up a count of how many times each expert participant agreed and disagreed
## with the speaker
nAgree <- rep(NaN,nExperts)
nDisagree <- rep(NaN,nExperts)

for(i in 1:nExperts){
  nAgree[i] <- sum(Agree_index[i,1:12])
  nDisagree[i] <- 12 - nAgree[i]
}

nExperts_agree <- length(which(nAgree>0))
nExperts_disagree <- length(which(nDisagree>0))

X_Main_agree <- matrix(0, nrow = sum(nAgree),ncol = nBeta)
for(i in 1:nrow(X_Main_agree)){
  X_Main_agree[i,Agree_data[i,4]] <- 1
}
Y_agree <- Agree_data[,2]
Participant_list_agree <- Agree_data[,3]

X_Main_disagree <- matrix(0, nrow = sum(nDisagree),ncol = nBeta)
for(i in 1:nrow(X_Main_disagree)){
  X_Main_disagree[i,Disagree_data[i,4]] <- 1
}

Y_disagree <- Disagree_data[,2]
Participant_list_disagree <- Disagree_data[,3]
Prior <- Disagree_data[,1]

tickerAgree = 1
tickerDisagree = 1

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

index_agree <- matrix(NaN, nrow = nExperts_agree, ncol = (nReps))
index_disagree <- matrix(NaN, nrow = nExperts_disagree, ncol = (nReps))

for(i in 1:nExperts_agree){
  for(j in 1:nAgree[i]){
    index_agree[i,j] <- j + sum(nAgree[1:i]) - nAgree[i]
  }
}

for(i in 1:nExperts_disagree){
  for(j in 1:nDisagree[i]){
    index_disagree[i,j] <- j + sum(nDisagree[1:i]) - nDisagree[i]
  }
}

Novice_length <- nNovices * nReps

## X and Y for novices. X is always just a 1 in the first column because these
## are all novice participants

X_nov <- matrix(0, nrow = Novice_length, ncol = nBeta)
X_nov[,1] <- rep(1,Novice_length)

Y_nov <- rep(NaN, Novice_length)
ticker = 1
for(i in 1:nrow(Inexperienced_replication)){
  for(j in 1:nReps){
    trialTicker = j + 49
    Y_nov[ticker] <- Inexperienced_replication[i,trialTicker]
    ticker = ticker+1
  }
}

# Expected model (only experienced varies by condition): run JAGS ----------------

nSamples = 4000

DIC <- rep(NaN,3)

X_agree <- X_Main_agree
Y_agree <- Y_agree
## Which participant produced each data point
Participant_list_agree <- Participant_list_agree
## How many generalizations did each participant agree with
nAgree <- nAgree

X_disagree <- X_Main_disagree
Y_disagree <- Y_disagree
## Which participant produced each data point
Participant_list_disagree <- Participant_list_disagree
## Prior estimates for each disagreement
Prior <- Prior
## How many generalizations did each participant disagree with
nDisagree <- nDisagree

X_nov <- X_nov
Y_nov <- Y_nov
nNovices <- nNovices
nTrials <- nReps

bPrec <- bPrec
nBeta <- nBeta
## Total experts
nExperts <- nExperts
## Number of experts who agreed and disagreed with any generalizations
nExperts_agree <- nExperts_agree
nExperts_disagree <- nExperts_disagree
index_agree <- index_agree
index_disagree <- index_disagree

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

Beta_add_phi_samples <- samples$BUGSoutput$sims.list$BETA

alpha_add_phi_samples <- samples$BUGSoutput$sims.list$alpha

phi_add_phi_samples <- samples$BUGSoutput$sims.list$psi_disagree

psi_add_phi_samples <- samples$BUGSoutput$sims.list$psi_disagree_bern

Phi_add_phi_samples <- samples$BUGSoutput$sims.list$phi_disagree

SigmaExp_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_add_phi_samples <- samples$BUGSoutput$sims.list$sigma_nov

Tau_add_phi_samples <- samples$BUGSoutput$sims.list$tau

summary_add_phi <- samples$BUGSoutput$summary

DIC[1] <- samples$BUGSoutput$DIC

DIC_standard[k] <- DIC[1]
}

## Test the fancy indexing matrix to make sure it's built properly
test_vector_agree <- rep(NaN, sum(nAgree))

ticker = 1
for(i in 1:nExperts_agree){
  for(j in 1:nAgree[i]){
    test_vector_agree[ticker] <- index_agree[i,j]
    ticker = ticker+1
  }
}

test_vector_disagree <- rep(NaN, sum(nDisagree))

ticker = 1
for(i in 1:nExperts_disagree){
  for(j in 1:nDisagree[i]){
    test_vector_disagree[ticker] <- index_disagree[i,j]
    ticker = ticker+1
  }
}

# Compare expert and novice sigmas ----------------------------------------

SigmaTest <- matrix(NaN,nrow = length(SigmaExp_add_phi_samples),ncol = 2)

SigmaDiff <- rep(NaN,length(SigmaExp_add_phi_samples))

for(i in 1:length(SigmaExp_add_phi_samples)){
  SigmaTest[i,1] <- Tau_add_phi_samples[i]^2/
    (Tau_add_phi_samples[i]^2 + SigmaExp_add_phi_samples[i]^2)
  SigmaTest[i,2] <- Tau_add_phi_samples[i]^2/
    (Tau_add_phi_samples[i]^2 + SigmaNov_add_phi_samples[i]^2)
  
  SigmaDiff[i] <- SigmaExp_add_phi_samples[i]-SigmaNov_add_phi_samples[i]
}

sortSigmaDiff <- sort(SigmaDiff)

ind975 = length(sortSigmaDiff)*.975
ind075 = length(sortSigmaDiff)*.075

mean(SigmaDiff)
sortSigmaDiff[ind975]
sortSigmaDiff[ind075]

sortSigmaTest <- matrix(NaN,nrow = length(SigmaExp_add_phi_samples),ncol = 2)

sortSigmaTest[,1] <- sort(SigmaTest[,1])
sortSigmaTest[,2] <- sort(SigmaTest[,2])

mean(SigmaTest[,1])
sortSigmaTest[ind975,1]
sortSigmaTest[ind075,1]

# Savage Dickey for sigma difference --------------------------------------

SigmaSimPulls1 <- rgamma(1000000, shape = 1.5, rate = 2)
SigmaSimPulls2 <- rgamma(1000000, shape = 1.5, rate = 2)

SigmaSimDiff <- SigmaSimPulls1-SigmaSimPulls2

hist(SigmaSimDiff)

SigmaDiffPrior_mean = mean(SigmaSimDiff)
SigmaDiffPrior_sd = sd(SigmaSimDiff)

prior_density_SigDiff = dnorm(0,0,.867)

posterior_density_sigDiff <- dlogspline(0,logspline(SigmaDiff))

SavDick_sig <- posterior_density_sigDiff/prior_density_SigDiff

# Fancy latent mixture plot for publication -------------------------------

X_prior_graph_latentMix_psi <- cbind(Prior, colMeans(psi_add_phi_samples))

disagree_plot = matrix(NaN,nrow = nrow(X_prior_graph_latentMix_psi),ncol = 5)

disagree_plot[,1] = X_prior_graph_latentMix_psi[,1]
disagree_plot[,2] = Y_disagree
disagree_plot[,3] = X_prior_graph_latentMix_psi[,2]

for(i in 1:nrow(X_prior_graph_latentMix_psi)){
  disagree_plot[i,4] = which(X_disagree[i,]==1)-1
  disagree_plot[i,5] = which(X_disagree[i,]==1)-1
  if(disagree_plot[i,4]>3){
    disagree_plot[i,5] <- 7-disagree_plot[i,4]
  }
}

par(mfrow=c(1,3))
par(mar = c(3.3, 3.3, 1, 1))

titles <- c("1A","1B","1C")

reverse_order <- c(3,2,1)

for(i in 1:(max(Agree_data[,4])-1)){
  id_comp <- i+1
  plot_index <- which(disagree_plot[,4]==reverse_order[i])
  
  plot(disagree_plot[plot_index,1],disagree_plot[plot_index,2],col = rgb(0,0,0),pch = 21,
       bg=rgb(disagree_plot[plot_index,3],disagree_plot[plot_index,3],disagree_plot[plot_index,3]),
       xlim = c(0,1), ylim = c(0,1), xlab = " ",ylab=" ", main= " ",
       axes = FALSE,cex = 2,type="n")
  
  lines(x=c(0,1),y=c(0,1),lty = 3,lwd =2)
  
  temp_lines <- Beta_add_phi_samples[,id_comp]
  
  text(x = .85, y = 1,labels = titles[i], cex = 2)
  
  for(j in 1:length(temp_lines)){
    lines(x = c(0,1), y = c(temp_lines[j],temp_lines[j]),col = rgb(0, 0, 0, max = 255, alpha = 1.2))
  }
  
  points(disagree_plot[plot_index,1],disagree_plot[plot_index,2],col = rgb(0,0,0),pch = 21,
         bg=rgb(disagree_plot[plot_index,3],disagree_plot[plot_index,3],disagree_plot[plot_index,3]),cex = 2)
  
  axis(side=1,pos=0,cex.axis = 1.6,lwd = 2)
  axis(side = 2, pos=0,cex.axis = 1.6,lwd = 2)
  
  title(ylab="Interpretation ( /100 games)", line=1.7, cex.lab=1.8)
  
  
  if(i > 2){
    title(xlab="Prior ( /100 games)", line=1.7, cex.lab=1.8)
  }
  
  if(i == 3){
    legend(x = .6, y = .3,col = c("black","black","darkgray"), cex = 1.7,lty = c(NaN,NaN,1),
           lwd = c(NaN,NaN,3), pch = c(21,21,NaN), pt.bg = c("black","white","white"), bty = 'n',
           legend = c(expression(paste(psi, " = ", 0)),
                      expression(paste(psi, " = ", 1)),
                      expression(beta[j])))
  }
  
}

# Bad model 1 (neither vary by condition): prep data -----------------------
nParticipants = nExperts+nNovices
nTrials = 12
nBeta = 2

X_One_agree <- matrix(data=rep(0,sum(nAgree)*nBeta),nrow = sum(nAgree)*nTrials,ncol = nBeta)

X_One_disagree <- matrix(data=rep(0,sum(nDisagree)*nBeta),nrow = sum(nDisagree)*nTrials,ncol = nBeta)

X_One_nov <- matrix(data=rep(0,nNovices*nTrials*nBeta),nrow = nNovices*nTrials,ncol = nBeta)


X_One_nov[,1] <- 1

X_One_agree[,2] <- 1
X_One_disagree[,2] <- 1

### MVN Precision Matrix (for raw JAGS betas)
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Bad model 1 (neither vary by condition): run JAGS -------------------------

X_agree <- X_One_agree
X_disagree <- X_One_disagree
X_nov <- X_One_nov

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

DIC[2] <- samples$BUGSoutput$DIC

# Bad model 2 (both vary by condition): prep data ------------------------

DIC_weird <- rep(NaN,iterations)

betas_weird <- matrix(NaN,nrow=iterations,ncol = 6)

for(k in 1:iterations){

nBeta = 6

X_nov_overfull <- matrix(data=rep(0,nNovices*nTrials*nBeta),nrow = nNovices*nTrials,ncol = nBeta)

tickerBothSplit = 1

for(i in 1:nrow(Inexperienced_replication)){
  for(j in 1:nTrials){
    typeTicker = j + 82
    temp = unlist(Inexperienced_replication[i,typeTicker])+2
    X_nov_overfull[tickerBothSplit,temp] <- 1
    
    tickerBothSplit = tickerBothSplit+1
  }
}

X_agree_overfull <- matrix(data=rep(0,sum(nAgree)*nBeta),nrow = sum(nAgree),ncol = nBeta)
X_agree_overfull[,4:6] <- X_Main_agree[,2:4]

X_disagree_overfull <- matrix(data=rep(0,sum(nDisagree)*nBeta),nrow = sum(nDisagree),ncol = nBeta)
X_disagree_overfull[,4:6] <- X_Main_disagree[,2:4]

### MVN Precision Matrix (for raw JAGS betas)
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Bad model 2 (both vary by condition): run JAGS ------------------------

nSamples = 4000

X_agree <- X_agree_overfull
X_disagree <- X_disagree_overfull
X_nov <- X_nov_overfull

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

Beta_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$BETA

SigmaExp_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_exp

SigmaCon_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_con

SigmaNov_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$sigma_nov

Tau_add_phi_samples_overfull <- samples$BUGSoutput$sims.list$tau

summary_add_phi_overfull <- samples$BUGSoutput$summary

DIC[3] <- samples$BUGSoutput$DIC

DIC_weird[k] <- DIC[3]

betas_weird[k,1:6] <- colMeans(Beta_add_phi_samples_overfull)

}

### Model comparison of (1) experienced varies, inexperienced doesn't 
### (2) neither varies (3) both do

# TVJ graph (collapse lane, combine with interpretation, split by audience expertise) ------------------------------------------------------------

par(mfrow=c(1,2))

final_TVJ_graph = matrix(NaN, 6, 3)

nGood = length(which(Experienced_replication[,111:122]==1))
nOk = length(which(Experienced_replication[,111:122]==0))
nBad = length(which(Experienced_replication[,111:122]==-1))

good_TVJ_split = matrix(NaN, nrow = nGood, ncol = 2)
ok_TVJ_split = matrix(NaN, nrow = nOk, ncol = 2)
bad_TVJ_split = matrix(NaN, nrow = nBad, ncol = 2)

good_ticker = 1
ok_ticker = 1
bad_ticker = 1

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:nReps){
    check_ticker = 110+j
    grab_ticker = 49+j
    if(Experienced_replication[i,check_ticker]==1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        good_TVJ_split[good_ticker,1] <- 0
      } else{good_TVJ_split[good_ticker,1] <- unlist(Experienced_replication[i,grab_ticker])}
      good_TVJ_split[good_ticker,2] <- unlist(Experienced_replication$Aud_expertise_index[i])
      good_ticker = good_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==0){
      if(Experienced_replication[i,grab_ticker]>1.1){
        ok_TVJ_split[ok_ticker,1] <- 0
      } else{ok_TVJ_split[ok_ticker,1] <- unlist(Experienced_replication[i,grab_ticker])}
      ok_TVJ_split[ok_ticker,2] <-  unlist(Experienced_replication$Aud_expertise_index[i])
      ok_ticker = ok_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==-1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        bad_TVJ_split[bad_ticker,1] <- 0
      } else{bad_TVJ_split[bad_ticker,1] <- unlist(Experienced_replication[i,grab_ticker])}
      bad_TVJ_split[bad_ticker,2] <-  unlist(Experienced_replication$Aud_expertise_index[i])
      bad_ticker = bad_ticker+1
    }
  }
}

good_TVJ_exp <- good_TVJ_split[which(good_TVJ_split[,2]==1),1]
good_TVJ_nov <- good_TVJ_split[which(good_TVJ_split[,2]==0),1]

ok_TVJ_exp <- ok_TVJ_split[which(ok_TVJ_split[,2]==1),1]
ok_TVJ_nov <- ok_TVJ_split[which(ok_TVJ_split[,2]==0),1]

bad_TVJ_exp <- bad_TVJ_split[which(bad_TVJ_split[,2]==1),1]
bad_TVJ_nov <- bad_TVJ_split[which(bad_TVJ_split[,2]==0),1]

nIterations = 1000
nSamples = 100

bootRange = BootstrapMean95(good_TVJ_exp,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[1,1] <- mean(good_TVJ_exp)
final_TVJ_graph[1,2:3] <- bootRange

bootRange = BootstrapMean95(ok_TVJ_exp,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[2,1] <- mean(ok_TVJ_exp)
final_TVJ_graph[2,2:3] <- bootRange

bootRange = BootstrapMean95(bad_TVJ_exp,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[3,1] <- mean(bad_TVJ_exp)
final_TVJ_graph[3,2:3] <- bootRange

bootRange = BootstrapMean95(good_TVJ_nov,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[4,1] <- mean(good_TVJ_nov)
final_TVJ_graph[4,2:3] <- bootRange

bootRange = BootstrapMean95(ok_TVJ_nov,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[5,1] <- mean(ok_TVJ_nov)
final_TVJ_graph[5,2:3] <- bootRange

bootRange = BootstrapMean95(bad_TVJ_nov,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[6,1] <- mean(bad_TVJ_nov)
final_TVJ_graph[6,2:3] <- bootRange

par(mar = c(3.5, 3, 1, 2))

exp_X <- c(.9,1.9,2.9)

plot(exp_X,final_TVJ_graph[1:3,1],ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'darkblue',bty = 'n',xlim = c(.75,3.25),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 16)

for(i in 1:3){
  arrows(x0 = exp_X[i], y0 = final_TVJ_graph[i,2],
         x1 = exp_X[i], y1 = final_TVJ_graph[i,3],code = 3,
         col = 'darkblue',angle = 90, length = .1,lty = 1,lwd = 2)
}

nov_X <- c(1.1,2.1,3.1)

points(nov_X,final_TVJ_graph[4:6,1], col = 'firebrick',cex = 2.2,lwd =2,pch = 15)

for(j in 1:3){
  i = j+3
  arrows(x0 = nov_X[j], y0 = final_TVJ_graph[i,2],
         x1 = nov_X[j], y1 = final_TVJ_graph[i,3],code = 3,
         col = 'firebrick',angle = 90, length = .1,lty = 1,lwd = 2)
}

xlabels = c('+','0','-')
xat = c(1:3)

title(main=" ", line=.5, cex.main=1.8)
title(ylab="Endorsement Rate", line=1.6, cex.lab=1.8)
title(xlab="Matchup Type", line=1.8, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,3.25),y=c(0,0),lwd =2)
lines(x = c(.75,3.25),y=c(.5,.5),lwd =2,lty = 2)
abline(v=0)

# Interpret graph (collapse lane, combine with TVJ) ------------------------------------------------------------

final_Interpret_graph = matrix(NaN, 3, 3)

good_Interpret = rep(NaN,nGood)
ok_Interpret = rep(NaN,nOk)
bad_Interpret = rep(NaN,nBad)

good_ticker = 1
ok_ticker = 1
bad_ticker = 1

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:nReps){
    check_ticker = 110+j
    grab_ticker = 77+j
    if(Experienced_replication[i,check_ticker]==1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        good_Interpret[good_ticker] <- 0
      } else{good_Interpret[good_ticker] <- Experienced_replication[i,grab_ticker]}
      good_ticker = good_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==0){
      if(Experienced_replication[i,grab_ticker]>1.1){
        ok_Interpret[ok_ticker] <- 0
      } else{ok_Interpret[ok_ticker] <- Experienced_replication[i,grab_ticker]}
      ok_ticker = ok_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==-1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        bad_Interpret[bad_ticker] <- 0
      } else{bad_Interpret[bad_ticker] <- Experienced_replication[i,grab_ticker]}
      bad_ticker = bad_ticker+1
    }
  }
}

good_Interpret <- unlist(good_Interpret)
ok_Interpret <- unlist(ok_Interpret)
bad_Interpret <- unlist(bad_Interpret)

nIterations = 1000
nSamples = 100

bootRange = BootstrapMean95(good_Interpret,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph[1,1] <- mean(good_Interpret)
final_Interpret_graph[1,2:3] <- bootRange

bootRange = BootstrapMean95(ok_Interpret,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph[2,1] <- mean(ok_Interpret)
final_Interpret_graph[2,2:3] <- bootRange

bootRange = BootstrapMean95(bad_Interpret,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph[3,1] <- mean(bad_Interpret)
final_Interpret_graph[3,2:3] <- bootRange

par(mar = c(3.5, 3, 1, 2))

exp_X <- c(.9,1.9,2.9)

plot(exp_X,final_Interpret_graph[,1],ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'darkblue',bty = 'n',xlim = c(.75,3.25),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 16)

for(i in 1:3){
  arrows(x0 = exp_X[i], y0 = final_Interpret_graph[i,2],
         x1 = exp_X[i], y1 = final_Interpret_graph[i,3],code = 3,
         col = 'darkblue',angle = 90, length = .1,lty = 1,lwd = 2)
}

xlabels = c('+','0','-')
xat = c(1:3)

title(main=" ", line=.5, cex.main=1.8)
title(ylab="Interpretation (/100 games)", line=1.6, cex.lab=1.8)
title(xlab="Matchup Type", line=1.8, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,3.25),y=c(0,0),lwd =2)
abline(v=0)

lines(x = c(.75,3.25),y = c(SONA_descriptives[1,1],SONA_descriptives[1,1]),lty = 2,lwd = 2, col = 'firebrick')
lines(x = c(.75,3.25),y = c(SONA_descriptives[2,1],SONA_descriptives[2,1]),lty = 3,lwd = 2, col = 'firebrick')
lines(x = c(.75,3.25),y = c(SONA_descriptives[3,1],SONA_descriptives[3,1]),lty = 3,lwd = 2, col = 'firebrick')

legend(1,.3,legend = c('Experienced','Inexperienced'),
       pch = c(16,15),lty = c(NA,2),lwd = c(NA,2),pt.cex = c(2,2),
       col = c('darkblue','firebrick'),
       cex = 1.5,bty = 'n')

##Add SONA responses by condition

final_Interpret_graph_SONA = matrix(NaN, 3, 3)

nGood_SONA = length(which(Inexperienced_replication[,83:94]==1))
nOk_SONA = length(which(Inexperienced_replication[,83:94]==0))
nBad_SONA = length(which(Inexperienced_replication[,83:94]==-1))

good_Interpret_SONA = rep(NaN,nGood_SONA)
ok_Interpret_SONA = rep(NaN,nOk_SONA)
bad_Interpret_SONA = rep(NaN,nBad_SONA)

good_ticker = 1
ok_ticker = 1
bad_ticker = 1

for(i in 1:nrow(Inexperienced_replication)){
  for(j in 1:nReps){
    check_ticker = 82+j
    grab_ticker = 49+j
    if(Inexperienced_replication[i,check_ticker]==1){
      if(Inexperienced_replication[i,grab_ticker]>1.1){
        good_Interpret_SONA[good_ticker] <- 0
      } else{good_Interpret_SONA[good_ticker] <- Inexperienced_replication[i,grab_ticker]}
      good_ticker = good_ticker+1
    }
    if(Inexperienced_replication[i,check_ticker]==0){
      if(Inexperienced_replication[i,grab_ticker]>1.1){
        ok_Interpret_SONA[ok_ticker] <- 0
      } else{ok_Interpret_SONA[ok_ticker] <- Inexperienced_replication[i,grab_ticker]}
      ok_ticker = ok_ticker+1
    }
    if(Inexperienced_replication[i,check_ticker]==-1){
      if(Inexperienced_replication[i,grab_ticker]>1.1){
        bad_Interpret_SONA[bad_ticker] <- 0
      } else{bad_Interpret_SONA[bad_ticker] <- Inexperienced_replication[i,grab_ticker]}
      bad_ticker = bad_ticker+1
    }
  }
}

good_Interpret_SONA <- unlist(good_Interpret_SONA)
ok_Interpret_SONA <- unlist(ok_Interpret_SONA)
bad_Interpret_SONA <- unlist(bad_Interpret_SONA)

nIterations = 1000
nSamples = 100

bootRange = BootstrapMean95(good_Interpret_SONA,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph_SONA[1,1] <- mean(good_Interpret_SONA)
final_Interpret_graph_SONA[1,2:3] <- bootRange

bootRange = BootstrapMean95(ok_Interpret_SONA,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph_SONA[2,1] <- mean(ok_Interpret_SONA)
final_Interpret_graph_SONA[2,2:3] <- bootRange

bootRange = BootstrapMean95(bad_Interpret_SONA,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph_SONA[3,1] <- mean(bad_Interpret_SONA)
final_Interpret_graph_SONA[3,2:3] <- bootRange

nov_X <- c(1.1,2.1,3.1)

points(nov_X,final_Interpret_graph_SONA[,1], col = 'firebrick',cex = 2.2,lwd =2,pch = 15)

for(i in 1:3){
  arrows(x0 = nov_X[i], y0 = final_Interpret_graph_SONA[i,2],
         x1 = nov_X[i], y1 = final_Interpret_graph_SONA[i,3],code = 3,
         col = 'firebrick',angle = 90, length = .1,lty = 1,lwd = 2)
}

# TVJ graph (collapse lane, combine with interpretation) ------------------------------------------------------------

par(mfrow=c(1,2))

final_TVJ_graph = matrix(NaN, 3, 3)

nGood = length(which(Experienced_replication[,111:122]==1))
nOk = length(which(Experienced_replication[,111:122]==0))
nBad = length(which(Experienced_replication[,111:122]==-1))

good_TVJ = rep(NaN,nGood)
ok_TVJ = rep(NaN,nOk)
bad_TVJ = rep(NaN,nBad)

good_ticker = 1
ok_ticker = 1
bad_ticker = 1

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:nReps){
    check_ticker = 110+j
    grab_ticker = 49+j
    if(Experienced_replication[i,check_ticker]==1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        good_TVJ[good_ticker] <- 0
      } else{good_TVJ[good_ticker] <- Experienced_replication[i,grab_ticker]}
      good_ticker = good_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==0){
      if(Experienced_replication[i,grab_ticker]>1.1){
        ok_TVJ[ok_ticker] <- 0
      } else{ok_TVJ[ok_ticker] <- Experienced_replication[i,grab_ticker]}
      ok_ticker = ok_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==-1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        bad_TVJ[bad_ticker] <- 0
      } else{bad_TVJ[bad_ticker] <- Experienced_replication[i,grab_ticker]}
      bad_ticker = bad_ticker+1
    }
  }
}

good_TVJ <- unlist(good_TVJ)
ok_TVJ <- unlist(ok_TVJ)
bad_TVJ <- unlist(bad_TVJ)

nIterations = 1000
nSamples = 100

bootRange = BootstrapMean95(good_TVJ,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[1,1] <- mean(good_TVJ)
final_TVJ_graph[1,2:3] <- bootRange

bootRange = BootstrapMean95(ok_TVJ,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[2,1] <- mean(ok_TVJ)
final_TVJ_graph[2,2:3] <- bootRange

bootRange = BootstrapMean95(bad_TVJ,nIterations = nIterations,nSamples = nSamples)
final_TVJ_graph[3,1] <- mean(bad_TVJ)
final_TVJ_graph[3,2:3] <- bootRange

par(mar = c(3.5, 3, 1, 2))

plot(1:3,final_TVJ_graph[,1],ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'black',bty = 'n',xlim = c(.75,3.25),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 16)

for(i in 1:3){
  arrows(x0 = i, y0 = final_TVJ_graph[i,2],
         x1 = i, y1 = final_TVJ_graph[i,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 1,lwd = 2)
}

xlabels = c('+','0','-')
xat = c(1:3)

title(main=" ", line=.5, cex.main=1.8)
title(ylab="Endorsement Rate", line=1.6, cex.lab=1.8)
title(xlab="Matchup Type", line=1.8, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,3.25),y=c(0,0),lwd =2)
lines(x = c(.75,3.25),y=c(.5,.5),lwd =2,lty = 2)
abline(v=0)

# Interpret graph (collapse lane, combine with TVJ) ------------------------------------------------------------

final_Interpret_graph = matrix(NaN, 3, 3)

good_Interpret = rep(NaN,nGood)
ok_Interpret = rep(NaN,nOk)
bad_Interpret = rep(NaN,nBad)

good_ticker = 1
ok_ticker = 1
bad_ticker = 1

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:nReps){
    check_ticker = 110+j
    grab_ticker = 77+j
    if(Experienced_replication[i,check_ticker]==1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        good_Interpret[good_ticker] <- 0
      } else{good_Interpret[good_ticker] <- Experienced_replication[i,grab_ticker]}
      good_ticker = good_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==0){
      if(Experienced_replication[i,grab_ticker]>1.1){
        ok_Interpret[ok_ticker] <- 0
      } else{ok_Interpret[ok_ticker] <- Experienced_replication[i,grab_ticker]}
      ok_ticker = ok_ticker+1
    }
    if(Experienced_replication[i,check_ticker]==-1){
      if(Experienced_replication[i,grab_ticker]>1.1){
        bad_Interpret[bad_ticker] <- 0
      } else{bad_Interpret[bad_ticker] <- Experienced_replication[i,grab_ticker]}
      bad_ticker = bad_ticker+1
    }
  }
}

good_Interpret <- unlist(good_Interpret)
ok_Interpret <- unlist(ok_Interpret)
bad_Interpret <- unlist(bad_Interpret)

nIterations = 1000
nSamples = 100

bootRange = BootstrapMean95(good_Interpret,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph[1,1] <- mean(good_Interpret)
final_Interpret_graph[1,2:3] <- bootRange

bootRange = BootstrapMean95(ok_Interpret,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph[2,1] <- mean(ok_Interpret)
final_Interpret_graph[2,2:3] <- bootRange

bootRange = BootstrapMean95(bad_Interpret,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph[3,1] <- mean(bad_Interpret)
final_Interpret_graph[3,2:3] <- bootRange

par(mar = c(3.5, 3, 1, 2))

exp_X <- c(.9,1.9,2.9)

plot(exp_X,final_Interpret_graph[,1],ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'darkblue',bty = 'n',xlim = c(.75,3.25),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 16)

for(i in 1:3){
  arrows(x0 = exp_X[i], y0 = final_Interpret_graph[i,2],
         x1 = exp_X[i], y1 = final_Interpret_graph[i,3],code = 3,
         col = 'darkblue',angle = 90, length = .1,lty = 1,lwd = 2)
}

xlabels = c('+','0','-')
xat = c(1:3)

title(main=" ", line=.5, cex.main=1.8)
title(ylab="Interpretation (/100 games)", line=1.6, cex.lab=1.8)
title(xlab="Matchup Type", line=1.8, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,3.25),y=c(0,0),lwd =2)
abline(v=0)

lines(x = c(.75,3.25),y = c(SONA_descriptives[1,1],SONA_descriptives[1,1]),lty = 2,lwd = 2, col = 'firebrick')
lines(x = c(.75,3.25),y = c(SONA_descriptives[2,1],SONA_descriptives[2,1]),lty = 3,lwd = 2, col = 'firebrick')
lines(x = c(.75,3.25),y = c(SONA_descriptives[3,1],SONA_descriptives[3,1]),lty = 3,lwd = 2, col = 'firebrick')

legend(1,.3,legend = c('Experienced','Inexperienced'),
       pch = c(16,15),lty = c(NA,2),lwd = c(NA,2),pt.cex = c(2,2),
       col = c('darkblue','firebrick'),
       cex = 1.5,bty = 'n')

##Add SONA responses by condition

final_Interpret_graph_SONA = matrix(NaN, 3, 3)

nGood_SONA = length(which(Inexperienced_replication[,83:94]==1))
nOk_SONA = length(which(Inexperienced_replication[,83:94]==0))
nBad_SONA = length(which(Inexperienced_replication[,83:94]==-1))

good_Interpret_SONA = rep(NaN,nGood_SONA)
ok_Interpret_SONA = rep(NaN,nOk_SONA)
bad_Interpret_SONA = rep(NaN,nBad_SONA)

good_ticker = 1
ok_ticker = 1
bad_ticker = 1

for(i in 1:nrow(Inexperienced_replication)){
  for(j in 1:nReps){
    check_ticker = 82+j
    grab_ticker = 49+j
    if(Inexperienced_replication[i,check_ticker]==1){
      if(Inexperienced_replication[i,grab_ticker]>1.1){
        good_Interpret_SONA[good_ticker] <- 0
      } else{good_Interpret_SONA[good_ticker] <- Inexperienced_replication[i,grab_ticker]}
      good_ticker = good_ticker+1
    }
    if(Inexperienced_replication[i,check_ticker]==0){
      if(Inexperienced_replication[i,grab_ticker]>1.1){
        ok_Interpret_SONA[ok_ticker] <- 0
      } else{ok_Interpret_SONA[ok_ticker] <- Inexperienced_replication[i,grab_ticker]}
      ok_ticker = ok_ticker+1
    }
    if(Inexperienced_replication[i,check_ticker]==-1){
      if(Inexperienced_replication[i,grab_ticker]>1.1){
        bad_Interpret_SONA[bad_ticker] <- 0
      } else{bad_Interpret_SONA[bad_ticker] <- Inexperienced_replication[i,grab_ticker]}
      bad_ticker = bad_ticker+1
    }
  }
}

good_Interpret_SONA <- unlist(good_Interpret_SONA)
ok_Interpret_SONA <- unlist(ok_Interpret_SONA)
bad_Interpret_SONA <- unlist(bad_Interpret_SONA)

nIterations = 1000
nSamples = 100

bootRange = BootstrapMean95(good_Interpret_SONA,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph_SONA[1,1] <- mean(good_Interpret_SONA)
final_Interpret_graph_SONA[1,2:3] <- bootRange

bootRange = BootstrapMean95(ok_Interpret_SONA,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph_SONA[2,1] <- mean(ok_Interpret_SONA)
final_Interpret_graph_SONA[2,2:3] <- bootRange

bootRange = BootstrapMean95(bad_Interpret_SONA,nIterations = nIterations,nSamples = nSamples)
final_Interpret_graph_SONA[3,1] <- mean(bad_Interpret_SONA)
final_Interpret_graph_SONA[3,2:3] <- bootRange

nov_X <- c(1.1,2.1,3.1)

points(nov_X,final_Interpret_graph_SONA[,1], col = 'firebrick',cex = 2.2,lwd =2,pch = 15)

for(i in 1:3){
  arrows(x0 = nov_X[i], y0 = final_Interpret_graph_SONA[i,2],
         x1 = nov_X[i], y1 = final_Interpret_graph_SONA[i,3],code = 3,
         col = 'firebrick',angle = 90, length = .1,lty = 1,lwd = 2)
}


# Compare B0 and Bj -------------------------------------------------------

prior_density = dnorm(0,0,.3537)
posterior_density <- rep(NaN,3)

for(i in 1:3){
  tempdif <- Beta_add_phi_samples[,1]-Beta_add_phi_samples[,i+1]
  posterior_density[i] <- dlogspline(0,logspline(tempdif))
}

SavDick <- data.frame(matrix(NaN,nrow=3,ncol=1))

SavDick <- posterior_density/prior_density
