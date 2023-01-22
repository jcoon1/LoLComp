# clears workspace:  
rm(list=ls()) 

setwd('C:/Users/Jeff/Documents/R/upload/LoLSolo_Aud')

#load(".RData")

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

# Get inexperienced participants' data ------------------------------------

Inexperienced_replication <- read_csv("Inexperienced_replication_csv.csv",col_names=TRUE)

nParticipants <- nrow(Inexperienced_replication)

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

nExperts = nrow(Experienced_replication)
nNovices = nrow(Inexperienced_replication)

### Prep data for main JAGS inputs
nParticipants = nExperts+nNovices
nTrials = 12
nExpResponses = nExperts*nTrials
nNovResponses = nNovices*nTrials
nBeta = 4

XSplit <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YSplit <- rep(NaN,nParticipants*nTrials)

XSplit[(nExpResponses+1):(nParticipants*nTrials),1] <- 1

tickerSplit = 1

temp <- rep(NaN,nrow(Experienced_replication))

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:nTrials){
    trialTicker = j + 77
    typeTicker = j + 110
    YSplit[tickerSplit] <- unlist(Experienced_replication[i,trialTicker])
    temp[i] = unlist(Experienced_replication[i,typeTicker])+3
    XSplit[tickerSplit,temp[i]] <- 1
    
    tickerSplit = tickerSplit+1
  }
}

for(i in 1:nrow(Inexperienced_replication)){
  for(j in 1:nTrials){
    trialTicker = j + 49
    
    YSplit[tickerSplit] <- unlist(Inexperienced_replication[i,trialTicker])
    
    tickerSplit = tickerSplit+1
  }
}

### MVN Precision Matrix (for raw JAGS betas)
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Expected model (only experienced varies by condition): run JAGS ----------------

nSamples = 4000

betaSplitsamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
compareDIC <- rep(NaN,3)

X <- XSplit
Y<- YSplit

data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

betaSplitsamples <- samples$BUGSoutput$sims.list$BETA

sigma_Splitsamples <- samples$BUGSoutput$sims.list$sigma

tau_Splitsamples <- samples$BUGSoutput$sims.list$tau

summarySplit <- samples$BUGSoutput$summary

compareDIC[1] <- samples$BUGSoutput$DIC


# Bad model 1 (neither vary by condition): prep data -----------------------
nParticipants = nExperts+nNovices
nTrials = 12
nExpResponses = nExperts*nTrials
nNovResponses = nNovices*nTrials
nBeta = 2

XOne <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)

XOne[1:nExpResponses,1] <- 1

XOne[(nExpResponses+1):(nParticipants*nTrials),1] <- 1

### MVN Precision Matrix (for raw JAGS betas)
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Bad model 1 (neither vary by condition): run JAGS -------------------------
betaOnesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XOne
Y<- YSplit

data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

betaOnesamples <- samples$BUGSoutput$sims.list$BETA

sigma_Onesamples <- samples$BUGSoutput$sims.list$sigma

tau_Onesamples <- samples$BUGSoutput$sims.list$tau

summaryOne <- samples$BUGSoutput$summary

compareDIC[2] <- samples$BUGSoutput$DIC



# Bad model 2 (both vary by condition): prep data ------------------------

nParticipants = nExperts+nNovices
nTrials = 12
nExpResponses = nExperts*nTrials
nNovResponses = nNovices*nTrials
nBeta = 6

XBothSplit <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)

tickerBothSplit = 1

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:nTrials){
    trialTicker = j + 77
    typeTicker = j + 110
    temp[i] = unlist(Experienced_replication[i,typeTicker])+3
    XBothSplit[tickerBothSplit,temp[i]] <- 1
    
    tickerBothSplit = tickerBothSplit+1
  }
}

temp <- rep(NaN, nrow(Inexperienced_replication))

for(i in 1:nrow(Inexperienced_replication)){
  for(j in 1:nTrials){
    trialTicker = j + 49
    typeTicker = j + 82
    temp[i] = unlist(Inexperienced_replication[i,typeTicker])+3
    XBothSplit[tickerBothSplit,temp[i]] <- 1
    
    tickerBothSplit = tickerBothSplit+1
  }
}

### MVN Precision Matrix (for raw JAGS betas)
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Bad model 2 (both vary by condition): run JAGS ------------------------

nSamples = 4000

betaBothSplitsamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XBothSplit
Y<- YSplit

data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_jags.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

betaBothSplitsamples <- samples$BUGSoutput$sims.list$BETA

sigma_BothSplitsamples <- samples$BUGSoutput$sims.list$sigma

tau_BothSplitsamples <- samples$BUGSoutput$sims.list$tau

summaryBothSplit <- samples$BUGSoutput$summary

compareDIC[3] <- samples$BUGSoutput$DIC

### Model comparison of (1) experienced varies, inexperienced doesn't 
### (2) neither varies (3) both do
compareDIC

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

# Agreement vs. Disagreement Graph ----------------------------------------
Agree_index <- matrix(NaN,nrow = nrow(Experienced_replication),ncol = 12)
Disagree_index <- matrix(NaN,nrow = nrow(Experienced_replication),ncol = 12)

for(i in 1:nrow(Agree_index)){
  for(j in 1:ncol(Agree_index)){
    TVJ_ticker = j + 49
    ref_ticker = j*4
    Agree <- Experienced_replication[i,TVJ_ticker] == 1
    Disagree <- Experienced_replication[i,TVJ_ticker] != 1
    Agree_index[i,j] <- Agree
    Disagree_index[i,j] <- Disagree
  }
}

Agree_data <- matrix(NaN,nrow = length(which(Agree_index==1)),ncol = 3)
Disagree_data <- matrix(NaN,nrow = length(which(Disagree_index==1)),ncol = 3)

Agree_ticker = 1
Disagree_ticker = 1

for(i in 1:nrow(Experienced_replication)){
  for(j in 1:ncol(Agree_index)){
    prior_ticker = j+65
    interpret_ticker = j+77
    if(Agree_index[i,j]==1){
      Agree_data[Agree_ticker,1] <- unlist(Experienced_replication[i,prior_ticker])
      Agree_data[Agree_ticker,2] <- unlist(Experienced_replication[i,interpret_ticker])
      Agree_data[Agree_ticker,3] <- i
      Agree_ticker = Agree_ticker+1
    }
    if(Disagree_index[i,j]==1){
      Disagree_data[Disagree_ticker,1] <- unlist(Experienced_replication[i,prior_ticker])
      Disagree_data[Disagree_ticker,2] <- unlist(Experienced_replication[i,interpret_ticker])
      Disagree_data[Disagree_ticker,3] <- i
      Disagree_ticker = Disagree_ticker+1
    }
  }
}

Diff_Agree <- Agree_data[,2] - Agree_data[,1]
Diff_Disagree <- Disagree_data[,2] - Disagree_data[,1]

c1 <- rgb(173,216,230,max = 255, alpha = 130, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 100, names = "lt.pink")

b <- min(c(Diff_Agree,Diff_Disagree)) - 0.001 # Set the minimum for the breakpoints
e <- max(c(Diff_Agree,Diff_Disagree)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 15) # Make a neat vector for the breakpoints
hist1 <- hist(Diff_Agree)
hist2 <- hist(Diff_Disagree)
par(mfrow=c(1,1))
#par(mar = c(5, 4, 1, 2))
plot(hist1,col = c1,ylim = c(0,200),xlim = c(-1,1),xlab = 'Interpretation - Prior Estimate',main = "Interpretation Shift",
     axes = FALSE)
plot(hist2,col = c2,add=TRUE)
abline(v=0,lty=3)
axis(side = 2, pos= -1)
axis(side = 1, pos=0)
legend(-.9,100,legend = c('Agreement','Disagreement'),pch = c(15,15),
       col = c(c1,c2),cex = 1.3)


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
  tempdif <- betaSplitsamples[,1]-betaSplitsamples[,i+1]
  posterior_density[i] <- dlogspline(0,logspline(tempdif))
}

SavDick <- data.frame(matrix(NaN,nrow=3,ncol=1))

SavDick <- posterior_density/prior_density

