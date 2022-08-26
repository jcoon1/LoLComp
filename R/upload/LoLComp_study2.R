# clears workspace:  
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load("LoLComp_study2.RData")

load("splitsig.RData")
library(purrr)

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
nItems = 24

SONAData <- read_csv("SONAData.csv",col_names=TRUE)

nNovices = max(SONAData$Participant)

totalNovice = nItems*nNovices


# Prep data for JAGS ------------------------------------------------------

nExp <- nExperts
nReps <- nItems
p <- c(OnlineData$EarlyPrior,OnlineData$LatePrior)
y <- c(OnlineData$EarlyTVJ,OnlineData$LateTVJ)

# Run JAGS model ----------------------------------------------------------

data <- list("nReps","y","p","nExp") # to be passed on to JAGS
# myinits <- list(
#   list("muGrand" = runif(1,0,1), "v" = 3,"prior_v" = 3))

# parameters to be monitored:	
parameters <- c("alphaGrand","betaGrand","alpha_sigma","beta_sigma","alphaDelta","betaDelta","alpha","beta")


#Now use jags.parallel to run multiple chains much quicker, adjust chains in n.chains

samples =jags.parallel(data,
                       #inits=inits.jagsParallel(), 
                       parameters.to.save = parameters,
                       model.file="hierarchical_aud_LoLComp.txt",n.chains=4,n.thin=1,n.burnin=1000,n.iter=4000)

# End parallel running section

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

samples_linking <- samples

alphaGrand_linking    <- samples_linking$BUGSoutput$sims.list$alphaGrand
betaGrand_linking <- samples_linking$BUGSoutput$sims.list$betaGrand

alphaDelta_linking    <- samples_linking$BUGSoutput$sims.list$alphaDelta
betaDelta_linking <- samples_linking$BUGSoutput$sims.list$betaDelta

alpha_linking    <- samples_linking$BUGSoutput$sims.list$alpha
beta_linking <- samples_linking$BUGSoutput$sims.list$beta

sigmaAlpha_linking <- samples_linking$BUGSoutput$sims.list$alpha_sigma

sigmaBeta_linking <- samples_linking$BUGSoutput$sims.list$beta_sigma

chainlength = nrow(alphaGrand_linking)/4

DIC_linking = samples_linking$BUGSoutput$DIC

# Calculations for inference table (95% CI's) -------------------------------------

lo = 12000*.05
hi = 12000*.95
inference_table <- matrix(NaN, nrow = 4, ncol = 3)

inference_table[1,1] <- mean(alphaGrand_linking)
tempsort <- sort(alphaGrand_linking)
inference_table[1,2] <- tempsort[lo]
inference_table[1,3] <- tempsort[hi]

inference_table[2,1] <- mean(betaGrand_linking)
tempsort <- sort(betaGrand_linking)
inference_table[2,2] <- tempsort[lo]
inference_table[2,3] <- tempsort[hi]

inference_table[3,1] <- mean(sigmaAlpha_linking)
tempsort <- sort(sigmaAlpha_linking)
inference_table[3,2] <- tempsort[lo]
inference_table[3,3] <- tempsort[hi]

inference_table[4,1] <- mean(sigmaBeta_linking)
tempsort <- sort(sigmaBeta_linking)
inference_table[4,2] <- tempsort[lo]
inference_table[4,3] <- tempsort[hi]

mean(sigmaAlpha_linking)

# Set up binned overlay for main plot -------------------------------------

nReps_nov = 12

exp_data <- matrix(NaN,nrow = (nReps*length(p)),ncol = 4)
for(i in 1:length(p)){
  for(j in 1:nReps){
    long_ticker = (i-1)*nReps_nov+j
    exp_data[long_ticker,1] <- p[long_ticker]
    exp_data[long_ticker,2] <- y[long_ticker]
    if(exp_data[long_ticker,2]==2){
      exp_data[long_ticker,2]<- 0
    }
  }
}

exp_overlay <- matrix(NaN,nrow = 10,ncol = 4)
exp_overlay[,1] <- seq(.05,.95,by = .1)

bin_size_exp = rep(NaN,10)

for(i in 1:nrow(exp_overlay)){
  tol = .0001
  lo = (i*.1)-.1
  hi = i * .1
  
  templo <- exp_data
  templo[,1] <- templo[,1]+tol
  temp_bin1 <- exp_data[which(templo[,1]>lo),1:2]
  
  temphi <- temp_bin1
  temphi[,1] <- temphi[,1]+tol
  if(i<10){
    temp_bin2 <- temp_bin1[which(temphi[,1]<hi),2]
  }
  else{temp_bin2 <- temp_bin1[,2]}
  bin_size_exp[i] = length(temp_bin2)
  
  exp_overlay[i,2] <- mean(temp_bin2)
  
  nSamples = 10
  nIterations = 1000
  
  if(i==1){nSamples = 2}
  
  bootRange <- BootstrapMean95(temp_bin2,nSamples,nIterations)
  
  exp_overlay[i,3:4] <- bootRange
}


# Calculate inexperienced interpretation for comparison -------------------

nParticipants <- max(SONAData$Participant)

nResponses_SONA <- nParticipants * nReps

SONA_vector <- c(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)],
                 SONAData$LateResponse[which(SONAData$SpeakerLate==1)])

SONA_descriptives <- matrix(NaN,nrow = 3,ncol = 1)

SONA_descriptives[1,1] <- mean(SONA_vector)

nIterations = 1000
nSamples = 100

SONA_descriptives[2:3,1] <- BootstrapMean95(SONA_vector,nSamples,nIterations)

# Alex plot for hierarchical by participant --------------------------------------------------

par(mfrow=c(1,1))

plot_index <- floor(runif(4000, min=1, max=12001))

q = seq(0,1,by = .01)
q_odds <- q/(1-q)
q_logodds <- log(q_odds)

temp_hierarchical <- q


plot(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,0,0.002)),frame.plot = FALSE,
     axes = FALSE,xlab = '',ylab = '')

xlabels = c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
xat = c(.1,.2,.3,.4,.5,.6,.7,.8,.9)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=0,cex.axis = 1.5,lwd = 2)

lines(x = c(0,1),y=c(0,0),lwd =2)

lines(x = c(0,1),y = c(.5,.5),lty = 4)
lines(x = c(.5,.5),y = c(0,1),lty = 4)


for(i in 2:length(alphaGrand_linking)){
  temp_hierarchical <- 1/(1 + exp(-betaGrand_linking[plot_index[i]]*(q_logodds-alphaGrand_linking[plot_index[i]])))
  lines(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,0,0.002)))
}

# plot_index <- floor(runif(4000, min=1, max=12001))

points(x=exp_overlay[,1],y=exp_overlay[,2],col = "black",pch = 16, cex = 2.2)

for(i in 1:nrow(exp_overlay)){
  if(exp_overlay[i,3]-exp_overlay[i,4]!=0){
    if(exp_overlay[i,3]==0){
      exp_overlay[i,3] <- exp_overlay[i,3]+.005
    }
    arrows(x0 = exp_overlay[i,1], y0 = exp_overlay[i,3],
           x1 = exp_overlay[i,1], y1 = exp_overlay[i,4],code = 3,
           col = "black",angle = 90, length = .1,lty = 2, lwd = 2)
  }
}

points(x = SONA_descriptives[1,1],y = c(.05),pch = 18,cex = 2.2,lwd =2)
arrows(x0 = SONA_descriptives[2,1], y0 = c(.05),
       x1 = SONA_descriptives[3,1], y1 = c(.05),code = 3,
       lwd = 2, angle = 90, length = .1,lty=1)

legend(.005,.98,legend = c('Empirical','Inferred',"Naive Interpretation"),
       pch = c(16,NA,18),lty = c(NA,1,NA),lwd = c(NA,2,2),pt.cex = c(2,NA,2),
       col = c('black',"black","black"),
       cex = 1.2,bty = 'n',
       y.intersp = .7)



title(main=" ", line=.5, cex.main=1.8)
title(ylab="Endorsement Rate", line=1.6, cex.lab=1.8)
title(xlab="Estimated Applicability", line=1.8, cex.lab=1.8)


# Marginal Simulation -----------------------------------------------------

nSim <- 10000

beta_x <- seq(.01,4,.01)

y <- (1.083*beta_x)/((4.339*beta_x^2)*(2.083*beta_x+1))

y_diff <- abs(y - .044)

beta_answer <- beta_x[which(y_diff==min(y_diff))]

alpha_answer <- beta_answer * 1.083

hist(c(p), freq = F,breaks = 20)
curve(dbeta(x, alpha_answer, beta_answer),from = 0, to =1,add=T)

curve(dbeta(x, 3.2, 2.8),from = 0, to =1,add=T)

input_prob <- rbeta(nSim,alpha_answer,beta_answer)

input_logodds <- log(input_prob/(1-input_prob))

output_prob <- rep(NaN,nSim*nSim)

ticker = 1
for(i in 1:nSim){
  for(j in 1:nSim){
    parameter_index <- round(runif(1,min = 1, max = nrow(alphaGrand_linking)))
    
    alpha_temp <- alphaGrand_linking[parameter_index,1]
    beta_temp <- betaGrand_linking[parameter_index,1]
    
    prob_gen <- 1/(1 + exp(-beta_temp*(input_logodds[i]-alpha_temp)))
    outcome <- rbernoulli(1,p = prob_gen)
    
    if(outcome==1){
      output_prob[ticker] <- input_prob[i]
      ticker = ticker+1
    }
  }
}

output_prob_clean_index <- min(which(is.na(output_prob)))-1

output_prob_clean <- output_prob[1:output_prob_clean_index]


# Marginal Sim part 2 -----------------------------------------------------

betaSplitsamples <- c(rawBetayesEarlysamples[,1],rawBetayesLatesamples[,1])

c1 <- rgb(173,216,230,max = 255, alpha = 130, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 100, names = "lt.pink")

SONA_comp <- c(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)],SONAData$LateResponse[which(SONAData$SpeakerLate==1)])

pick_output <- round(runif(length(SONA_comp),min=1,max=length(output_prob_clean)))

output_prob_plot <- output_prob_clean[pick_output]

b <- min(c(output_prob_plot,betaSplitsamples)) - 0.001 # Set the minimum for the breakpoints
e <- max(c(output_prob_plot,betaSplitsamples)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 100) # Make a neat vector for the breakpoints
hist1 <- hist(output_prob_plot)
# hist2 <- hist(betaSplitsamples[,1])
hist2 <- hist(c(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)],SONAData$LateResponse[which(SONAData$SpeakerLate==1)]))
par(mfrow=c(1,1))
#par(mar = c(5, 4, 1, 2))
plot(hist1,col = c1,ylim = c(0,200),xlim = c(0,1),xlab = 'Interpretation',main = "Speaker vs. Listener",
     axes = FALSE)
plot(hist2,col = c2,add=TRUE)
abline(v=0,lty=3)
axis(side = 2, pos= 0)
axis(side = 1, pos=0)

nSim = 10000

counter = 1

for(i in 1:nSim){
  speaker_sample <- round(runif(1,min=1,max=length(output_prob_plot)))
  listener_sample <- round(runif(1,min=1,max=length(betaSplitsamples)))
  
  if(output_prob_plot[speaker_sample]<betaSplitsamples[listener_sample]){
    counter = counter+1
  }
}

prob_marginal = counter/nSim

odds_marginal = prob_marginal/(1-prob_marginal)

odds_prior = .5/(1-.5)

BF_marginal = odds_marginal/odds_prior


# Inexperienced distance simulation (all items) --------------------------------------------

temp_yes_early <- OnlineData[OnlineData$SpeakerEarly==1,]

temp_yes_late <- OnlineData[OnlineData$SpeakerLate==1,]

exp_statement_pool <- c(temp_yes_early$EarlyPrior[which(temp_yes_early$EarlyTVJ==1)],
                        temp_yes_late$LatePrior[which(temp_yes_late$LateTVJ==1)])

exp_interpret_pool <- c(temp_yes_early$EarlyResponse[which(temp_yes_early$EarlyTVJ==1)],
                        temp_yes_late$LateResponse[which(temp_yes_late$LateTVJ==1)])

plot(x=exp_statement_pool,y=exp_interpret_pool,xlim= c(0,1),ylim=c(0,1),xlab='Prior',ylab='Interpretation',main='Experienced Belief to Interpretation',col = 'darkgreen')

abline(a=0,b=1)

abline(h=.5,lty=3)
abline(v=.5,lty=3)

exp_statement_pool_no <- c(temp_yes_early$EarlyPrior[which(temp_yes_early$EarlyTVJ==0)],
                           temp_yes_late$LatePrior[which(temp_yes_late$LateTVJ==0)])

exp_interpret_pool_no <- c(temp_yes_early$EarlyResponse[which(temp_yes_early$EarlyTVJ==0)],
                           temp_yes_late$LateResponse[which(temp_yes_late$LateTVJ==0)])

points(x=exp_statement_pool_no,y=exp_interpret_pool_no,col = 'red')

nov_interpret_pool <- c(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)],
                        SONAData$LateResponse[which(SONAData$SpeakerLate==1)])

nPulls <- 300

diff_raw <- matrix(nrow = length(exp_statement_pool),ncol = nPulls)

for(i in 1:length(exp_statement_pool)){
  for(j in 1:nPulls){
    nov_index <- round(runif(1,min = 1, max = length(nov_interpret_pool)))
    diff_raw[i,j] <- nov_interpret_pool[nov_index] - exp_statement_pool[i]
  }
}

diff_plot <- matrix(nrow = length(exp_statement_pool),ncol = 4)

diff_plot[,1] <- exp_statement_pool

nSamples = 50

for(i in 1:length(exp_statement_pool)){
  diff_plot[i,2] <- mean(diff_raw[i,])
  diff_plot[i,3:4] <- BootstrapMean95(diff_raw[i,],nSamples,nIterations)
}
plot(diff_plot[,1],diff_plot[,2],xlim = c(0,1),ylim = c(-1,1),xlab = "Speaker Belief",ylab="Listener - Speaker", col = rgb(0,0,0,.08))

for(i in 1:nrow(diff_plot)){
  arrows(x0 = diff_plot[i,1], y0 = diff_plot[i,3],
        x1 = diff_plot[i,1], y1 = diff_plot[i,4],code = 3,
        lwd = 1, angle = 90, length = .1,lty=1, col = rgb(0,0,0,.08))
}

abline(v = .5,lty = 2)
abline(h = 0, lty = 2)

abline(a = .7, b = -1)

abs_diff_raw <- matrix(nrow = length(exp_statement_pool),ncol = nPulls)

for(i in 1:length(exp_statement_pool)){
  for(j in 1:nPulls){
    nov_index <- round(runif(1,min = 1, max = length(nov_interpret_pool)))
    abs_diff_raw[i,j] <- abs(nov_interpret_pool[nov_index] - exp_statement_pool[i])
  }
}

abs_diff_plot <- matrix(nrow = length(exp_statement_pool),ncol = 4)

abs_diff_plot[,1] <- exp_statement_pool

nSamples = 50

for(i in 1:length(exp_statement_pool)){
  abs_diff_plot[i,2] <- mean(abs_diff_raw[i,])
  abs_diff_plot[i,3:4] <- BootstrapMean95(abs_diff_raw[i,],nSamples,nIterations)
}
plot(abs_diff_plot[,1],abs_diff_plot[,2],xlim = c(0,1),ylim = c(0,1),xlab = "Speaker Belief",ylab="abs(Listener - Speaker)", col = rgb(0,0,0,.08))

for(i in 1:nrow(abs_diff_plot)){
  arrows(x0 = abs_diff_plot[i,1], y0 = abs_diff_plot[i,3],
         x1 = abs_diff_plot[i,1], y1 = abs_diff_plot[i,4],code = 3,
         lwd = 1, angle = 90, length = .1,lty=1, col = rgb(0,0,0,.08))
}

abline(v = .5,lty = 2)
                           
hist(diff_raw)
abline(v = mean(diff_raw),lty = 2)


# Inexperienced Distance Simulation (by item) -----------------------------

sanity_check_early <- matrix(NaN, ncol = 6,nrow = 8)
sanity_check_late <- matrix(NaN, ncol = 6,nrow = 8)

comp1_speaker_pool_temp <- OnlineData[which(OnlineData$CompType==1),]
comp1_speaker_pool_early <- comp1_speaker_pool_temp$EarlyPrior[which(comp1_speaker_pool_temp$EarlyTVJ==1)]
comp1_speaker_pool_late <- comp1_speaker_pool_temp$LatePrior[which(comp1_speaker_pool_temp$LateTVJ==1)]
comp1_speaker_par_early <- comp1_speaker_pool_temp$Participant[which(comp1_speaker_pool_temp$EarlyTVJ==1)]
comp1_speaker_par_late <- comp1_speaker_pool_temp$Participant[which(comp1_speaker_pool_temp$LateTVJ==1)]
comp1_full_pool_early <- cbind(comp1_speaker_pool_early,comp1_speaker_par_early)
comp1_full_pool_late <- cbind(comp1_speaker_pool_late,comp1_speaker_par_late)
sanity_check_early[1,1]<- mean(comp1_full_pool_early[,1])
sanity_check_late[1,1]<- mean(comp1_full_pool_late[,1])
sanity_check_early[4,1] <- length(comp1_speaker_pool_early)
sanity_check_late[4,1] <- length(comp1_speaker_pool_late)

comp2_speaker_pool_temp <- OnlineData[which(OnlineData$CompType==2),]
comp2_speaker_pool_early <- comp2_speaker_pool_temp$EarlyPrior[which(comp2_speaker_pool_temp$EarlyTVJ==1)]
comp2_speaker_pool_late <- comp2_speaker_pool_temp$LatePrior[which(comp2_speaker_pool_temp$LateTVJ==1)]
comp2_speaker_par_early <- comp2_speaker_pool_temp$Participant[which(comp2_speaker_pool_temp$EarlyTVJ==1)]
comp2_speaker_par_late <- comp2_speaker_pool_temp$Participant[which(comp2_speaker_pool_temp$LateTVJ==1)]
comp2_full_pool_early <- cbind(comp2_speaker_pool_early,comp2_speaker_par_early)
comp2_full_pool_late <- cbind(comp2_speaker_pool_late,comp2_speaker_par_late)
sanity_check_early[1,2]<- mean(comp2_full_pool_early[,1])
sanity_check_late[1,2]<- mean(comp2_full_pool_late[,1])
sanity_check_early[4,2] <- length(comp2_speaker_pool_early)
sanity_check_late[4,2] <- length(comp2_speaker_pool_late)

comp3_speaker_pool_temp <- OnlineData[which(OnlineData$CompType==3),]
comp3_speaker_pool_early <- comp3_speaker_pool_temp$EarlyPrior[which(comp3_speaker_pool_temp$EarlyTVJ==1)]
comp3_speaker_pool_late <- comp3_speaker_pool_temp$LatePrior[which(comp3_speaker_pool_temp$LateTVJ==1)]
comp3_speaker_par_early <- comp3_speaker_pool_temp$Participant[which(comp3_speaker_pool_temp$EarlyTVJ==1)]
comp3_speaker_par_late <- comp3_speaker_pool_temp$Participant[which(comp3_speaker_pool_temp$LateTVJ==1)]
comp3_full_pool_early <- cbind(comp3_speaker_pool_early,comp3_speaker_par_early)
comp3_full_pool_late <- cbind(comp3_speaker_pool_late,comp3_speaker_par_late)
sanity_check_early[1,3]<- mean(comp3_full_pool_early[,1])
sanity_check_late[1,3]<- mean(comp3_full_pool_late[,1])
sanity_check_early[4,3] <- length(comp3_speaker_pool_early)
sanity_check_late[4,3] <- length(comp3_speaker_pool_late)

comp4_speaker_pool_temp <- OnlineData[which(OnlineData$CompType==4),]
comp4_speaker_pool_early <- comp4_speaker_pool_temp$EarlyPrior[which(comp4_speaker_pool_temp$EarlyTVJ==1)]
comp4_speaker_pool_late <- comp4_speaker_pool_temp$LatePrior[which(comp4_speaker_pool_temp$LateTVJ==1)]
comp4_speaker_par_early <- comp4_speaker_pool_temp$Participant[which(comp4_speaker_pool_temp$EarlyTVJ==1)]
comp4_speaker_par_late <- comp4_speaker_pool_temp$Participant[which(comp4_speaker_pool_temp$LateTVJ==1)]
comp4_full_pool_early <- cbind(comp4_speaker_pool_early,comp4_speaker_par_early)
comp4_full_pool_late <- cbind(comp4_speaker_pool_late,comp4_speaker_par_late)
sanity_check_early[1,4]<- mean(comp4_full_pool_early[,1])
sanity_check_late[1,4]<- mean(comp4_full_pool_late[,1])
sanity_check_early[4,4] <- length(comp4_speaker_pool_early)
sanity_check_late[4,4] <- length(comp4_speaker_pool_late)

comp5_speaker_pool_temp <- OnlineData[which(OnlineData$CompType==5),]
comp5_speaker_pool_early <- comp5_speaker_pool_temp$EarlyPrior[which(comp5_speaker_pool_temp$EarlyTVJ==1)]
comp5_speaker_pool_late <- comp5_speaker_pool_temp$LatePrior[which(comp5_speaker_pool_temp$LateTVJ==1)]
comp5_speaker_par_early <- comp5_speaker_pool_temp$Participant[which(comp5_speaker_pool_temp$EarlyTVJ==1)]
comp5_speaker_par_late <- comp5_speaker_pool_temp$Participant[which(comp5_speaker_pool_temp$LateTVJ==1)]
comp5_full_pool_early <- cbind(comp5_speaker_pool_early,comp5_speaker_par_early)
comp5_full_pool_late <- cbind(comp5_speaker_pool_late,comp5_speaker_par_late)
sanity_check_early[1,5]<- mean(comp5_full_pool_early[,1])
sanity_check_late[1,5]<- mean(comp5_full_pool_late[,1])
sanity_check_early[4,5] <- length(comp5_speaker_pool_early)
sanity_check_late[4,5] <- length(comp5_speaker_pool_late)

comp6_speaker_pool_temp <- OnlineData[which(OnlineData$CompType==6),]
comp6_speaker_pool_early <- comp6_speaker_pool_temp$EarlyPrior[which(comp6_speaker_pool_temp$EarlyTVJ==1)]
comp6_speaker_pool_late <- comp6_speaker_pool_temp$LatePrior[which(comp6_speaker_pool_temp$LateTVJ==1)]
comp6_speaker_par_early <- comp6_speaker_pool_temp$Participant[which(comp6_speaker_pool_temp$EarlyTVJ==1)]
comp6_speaker_par_late <- comp6_speaker_pool_temp$Participant[which(comp6_speaker_pool_temp$LateTVJ==1)]
comp6_full_pool_early <- cbind(comp6_speaker_pool_early,comp6_speaker_par_early)
comp6_full_pool_late <- cbind(comp6_speaker_pool_late,comp6_speaker_par_late)
sanity_check_early[1,6]<- mean(comp6_full_pool_early[,1])
sanity_check_late[1,6]<- mean(comp6_full_pool_late[,1])
sanity_check_early[4,6] <- length(comp6_speaker_pool_early)
sanity_check_late[4,6] <- length(comp6_speaker_pool_late)

nPulls = 1000
comp1_early_diffs_nov <- matrix(NaN, nrow = nrow(comp1_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp1_listener_pool_temp<- SONAData[which(SONAData$CompType==1),]
listener_pool <- comp1_listener_pool_temp$EarlyResponse[which(comp1_listener_pool_temp$SpeakerEarly==1)]
sanity_check_early[2,1]<- mean(listener_pool)
for(i in 1:nrow(comp1_full_pool_early)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    comp1_early_diffs_nov[i,j] <- listener_pool[pull_index] - comp1_full_pool_early[i]
  }
}

comp1_late_diffs_nov <- matrix(NaN, nrow = nrow(comp1_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp1_listener_pool_temp<- SONAData[which(SONAData$CompType==1),]
listener_pool <- comp1_listener_pool_temp$LateResponse[which(comp1_listener_pool_temp$SpeakerLate==1)]
sanity_check_late[2,1]<- mean(listener_pool)
for(i in 1:nrow(comp1_full_pool_late)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp1_late_diffs_nov[i,] <- temp_pulls-comp1_full_pool_late[i]
}

comp2_early_diffs_nov <- matrix(NaN, nrow = nrow(comp2_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp2_listener_pool_temp<- SONAData[which(SONAData$CompType==2),]
listener_pool <- comp2_listener_pool_temp$EarlyResponse[which(comp2_listener_pool_temp$SpeakerEarly==1)]
sanity_check_early[2,2]<- mean(listener_pool)
for(i in 1:nrow(comp2_full_pool_early)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp2_early_diffs_nov[i,] <- temp_pulls-comp2_full_pool_early[i]
}

comp2_late_diffs_nov <- matrix(NaN, nrow = nrow(comp2_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp2_listener_pool_temp<- SONAData[which(SONAData$CompType==2),]
listener_pool <- comp2_listener_pool_temp$LateResponse[which(comp2_listener_pool_temp$SpeakerLate==1)]
sanity_check_late[2,2]<- mean(listener_pool)
for(i in 1:nrow(comp2_full_pool_late)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp2_late_diffs_nov[i,] <- temp_pulls-comp2_full_pool_late[i]
}

comp3_early_diffs_nov <- matrix(NaN, nrow = nrow(comp3_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp3_listener_pool_temp<- SONAData[which(SONAData$CompType==3),]
listener_pool <- comp3_listener_pool_temp$EarlyResponse[which(comp3_listener_pool_temp$SpeakerEarly==1)]
sanity_check_early[2,3]<- mean(listener_pool)
for(i in 1:nrow(comp3_full_pool_early)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp3_early_diffs_nov[i,] <- temp_pulls-comp3_full_pool_early[i]
}

comp3_late_diffs_nov <- matrix(NaN, nrow = nrow(comp3_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp3_listener_pool_temp<- SONAData[which(SONAData$CompType==3),]
listener_pool <- comp3_listener_pool_temp$LateResponse[which(comp3_listener_pool_temp$SpeakerLate==1)]
sanity_check_late[2,3]<- mean(listener_pool)
for(i in 1:nrow(comp3_full_pool_late)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp3_late_diffs_nov[i,] <- temp_pulls-comp3_full_pool_late[i]
}

comp4_early_diffs_nov <- matrix(NaN, nrow = nrow(comp4_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp4_listener_pool_temp<- SONAData[which(SONAData$CompType==4),]
listener_pool <- comp4_listener_pool_temp$EarlyResponse[which(comp4_listener_pool_temp$SpeakerEarly==1)]
sanity_check_early[2,4]<- mean(listener_pool)
for(i in 1:nrow(comp4_full_pool_early)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp4_early_diffs_nov[i,] <- temp_pulls-comp4_full_pool_early[i]
}

comp4_late_diffs_nov <- matrix(NaN, nrow = nrow(comp4_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp4_listener_pool_temp<- SONAData[which(SONAData$CompType==4),]
listener_pool <- comp4_listener_pool_temp$LateResponse[which(comp4_listener_pool_temp$SpeakerLate==1)]
sanity_check_late[2,4]<- mean(listener_pool)
for(i in 1:nrow(comp4_full_pool_late)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp4_late_diffs_nov[i,] <- temp_pulls-comp4_full_pool_late[i]
}

comp5_early_diffs_nov <- matrix(NaN, nrow = nrow(comp5_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp5_listener_pool_temp<- SONAData[which(SONAData$CompType==5),]
listener_pool <- comp5_listener_pool_temp$EarlyResponse[which(comp5_listener_pool_temp$SpeakerEarly==1)]
sanity_check_early[2,5]<- mean(listener_pool)
for(i in 1:nrow(comp5_full_pool_early)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp5_early_diffs_nov[i,] <- temp_pulls-comp5_full_pool_early[i]
}

comp5_late_diffs_nov <- matrix(NaN, nrow = nrow(comp5_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp5_listener_pool_temp<- SONAData[which(SONAData$CompType==5),]
listener_pool <- comp5_listener_pool_temp$LateResponse[which(comp5_listener_pool_temp$SpeakerLate==1)]
sanity_check_late[2,5]<- mean(listener_pool)
for(i in 1:nrow(comp5_full_pool_late)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp5_late_diffs_nov[i,] <- temp_pulls-comp5_full_pool_late[i]
}

comp6_early_diffs_nov <- matrix(NaN, nrow = nrow(comp6_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp6_listener_pool_temp<- SONAData[which(SONAData$CompType==6),]
listener_pool <- comp6_listener_pool_temp$EarlyResponse[which(comp6_listener_pool_temp$SpeakerEarly==1)]
sanity_check_early[2,6]<- mean(listener_pool)
for(i in 1:nrow(comp6_full_pool_early)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp6_early_diffs_nov[i,] <- temp_pulls-comp6_full_pool_early[i]
}

comp6_late_diffs_nov <- matrix(NaN, nrow = nrow(comp6_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
comp6_listener_pool_temp<- SONAData[which(SONAData$CompType==6),]
listener_pool <- comp6_listener_pool_temp$LateResponse[which(comp6_listener_pool_temp$SpeakerLate==1)]
sanity_check_late[2,6]<- mean(listener_pool)
for(i in 1:nrow(comp6_full_pool_late)){
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp6_late_diffs_nov[i,] <- temp_pulls-comp6_full_pool_late[i]
}

sanity_check_early[3,] <- sanity_check_early[2,]-sanity_check_early[1,]
sanity_check_late[3,] <- sanity_check_late[2,]-sanity_check_late[1,]

sanity_check_early[5,] <- sanity_check_early[3,] * sanity_check_early[4,]
sum(sanity_check_early[5,])/sum(sanity_check_early[4,])

sanity_check_late[5,] <- sanity_check_late[3,] * sanity_check_late[4,]
sum(sanity_check_late[5,])/sum(sanity_check_late[4,])

sanity_check_early[6,] <- mean(c(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)],
                                 SONAData$LateResponse[which(SONAData$SpeakerLate==1)]))

sanity_check_late[6,] <- mean(c(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)],
                                 SONAData$LateResponse[which(SONAData$SpeakerLate==1)]))

sanity_check_early[7,] <- sanity_check_early[6,]-sanity_check_early[1,]
sanity_check_late[7,] <- sanity_check_late[6,]-sanity_check_late[1,]

sanity_check_early[8,] <- sanity_check_early[7,] * sanity_check_early[4,]
sum(sanity_check_early[8,])/sum(sanity_check_early[4,])

sanity_check_late[8,] <- sanity_check_late[7,] * sanity_check_late[4,]
sum(sanity_check_late[8,])/sum(sanity_check_late[4,])


par(mfrow=c(3,2))

hist(comp1_early_diffs_nov)
hist(comp2_early_diffs_nov)
hist(comp3_early_diffs_nov)
hist(comp4_early_diffs_nov)
hist(comp5_early_diffs_nov)
hist(comp6_early_diffs_nov)

hist(comp1_late_diffs_nov)
hist(comp2_late_diffs_nov)
hist(comp3_late_diffs_nov)
hist(comp4_late_diffs_nov)
hist(comp5_late_diffs_nov)
hist(comp6_late_diffs_nov)

all_nov_diff <- c(comp1_late_diffs_nov,comp2_late_diffs_nov,comp3_late_diffs_nov,
                 comp4_late_diffs_nov,comp5_late_diffs_nov,comp6_late_diffs_nov,
                 comp1_early_diffs_nov,comp2_early_diffs_nov,comp3_early_diffs_nov,
                 comp4_early_diffs_nov,comp5_early_diffs_nov,comp6_early_diffs_nov)

all_nov_diff_desc <- c(mean(all_nov_diff),sd(all_nov_diff))

# Experienced Distance Simulation -----------------------------------------

nPulls = 1000
comp1_early_diffs_exp <- matrix(NaN, nrow = nrow(comp1_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp1_full_pool_early)){
  temp_listener_pool <- comp1_speaker_pool_temp[which(comp1_speaker_pool_temp$Participant!=comp1_full_pool_early[i,2]),]
  listener_pool <- temp_listener_pool$EarlyResponse[which(temp_listener_pool$SpeakerEarly==1)]
    for(j in 1:nPulls){
      pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
      temp_pulls[j] <- listener_pool[pull_index]
    }
  comp1_early_diffs_exp[i,] <- temp_pulls-comp1_full_pool_early[i]
}

comp1_late_diffs_exp <- matrix(NaN, nrow = nrow(comp1_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp1_full_pool_late)){
  temp_listener_pool <- comp1_speaker_pool_temp[which(comp1_speaker_pool_temp$Participant!=comp1_full_pool_late[i,2]),]
  listener_pool <- temp_listener_pool$LateResponse[which(temp_listener_pool$SpeakerLate==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp1_late_diffs_exp[i,] <- temp_pulls-comp1_full_pool_late[i]
}

comp2_early_diffs_exp <- matrix(NaN, nrow = nrow(comp2_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp2_full_pool_early)){
  temp_listener_pool <- comp2_speaker_pool_temp[which(comp2_speaker_pool_temp$Participant!=comp2_full_pool_early[i,2]),]
  listener_pool <- temp_listener_pool$EarlyResponse[which(temp_listener_pool$SpeakerEarly==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp2_early_diffs_exp[i,] <- temp_pulls-comp2_full_pool_early[i]
}

comp2_late_diffs_exp <- matrix(NaN, nrow = nrow(comp2_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp2_full_pool_late)){
  temp_listener_pool <- comp2_speaker_pool_temp[which(comp2_speaker_pool_temp$Participant!=comp2_full_pool_late[i,2]),]
  listener_pool <- temp_listener_pool$LateResponse[which(temp_listener_pool$SpeakerLate==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp2_late_diffs_exp[i,] <- temp_pulls-comp2_full_pool_late[i]
}

comp3_early_diffs_exp <- matrix(NaN, nrow = nrow(comp3_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp3_full_pool_early)){
  temp_listener_pool <- comp3_speaker_pool_temp[which(comp3_speaker_pool_temp$Participant!=comp3_full_pool_early[i,2]),]
  listener_pool <- temp_listener_pool$EarlyResponse[which(temp_listener_pool$SpeakerEarly==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp3_early_diffs_exp[i,] <- temp_pulls-comp3_full_pool_early[i]
}

comp3_late_diffs_exp <- matrix(NaN, nrow = nrow(comp3_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp3_full_pool_late)){
  temp_listener_pool <- comp3_speaker_pool_temp[which(comp3_speaker_pool_temp$Participant!=comp3_full_pool_late[i,2]),]
  listener_pool <- temp_listener_pool$LateResponse[which(temp_listener_pool$SpeakerLate==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp3_late_diffs_exp[i,] <- temp_pulls-comp3_full_pool_late[i]
}

comp4_early_diffs_exp <- matrix(NaN, nrow = nrow(comp4_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp4_full_pool_early)){
  temp_listener_pool <- comp4_speaker_pool_temp[which(comp4_speaker_pool_temp$Participant!=comp4_full_pool_early[i,2]),]
  listener_pool <- temp_listener_pool$EarlyResponse[which(temp_listener_pool$SpeakerEarly==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp4_early_diffs_exp[i,] <- temp_pulls-comp4_full_pool_early[i]
}

comp4_late_diffs_exp <- matrix(NaN, nrow = nrow(comp4_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp4_full_pool_late)){
  temp_listener_pool <- comp4_speaker_pool_temp[which(comp4_speaker_pool_temp$Participant!=comp4_full_pool_late[i,2]),]
  listener_pool <- temp_listener_pool$LateResponse[which(temp_listener_pool$SpeakerLate==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp4_late_diffs_exp[i,] <- temp_pulls-comp4_full_pool_late[i]
}

comp5_early_diffs_exp <- matrix(NaN, nrow = nrow(comp5_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp5_full_pool_early)){
  temp_listener_pool <- comp5_speaker_pool_temp[which(comp5_speaker_pool_temp$Participant!=comp5_full_pool_early[i,2]),]
  listener_pool <- temp_listener_pool$EarlyResponse[which(temp_listener_pool$SpeakerEarly==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp5_early_diffs_exp[i,] <- temp_pulls-comp5_full_pool_early[i]
}

comp5_late_diffs_exp <- matrix(NaN, nrow = nrow(comp5_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp5_full_pool_late)){
  temp_listener_pool <- comp5_speaker_pool_temp[which(comp5_speaker_pool_temp$Participant!=comp5_full_pool_late[i,2]),]
  listener_pool <- temp_listener_pool$LateResponse[which(temp_listener_pool$SpeakerLate==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp5_late_diffs_exp[i,] <- temp_pulls-comp5_full_pool_late[i]
}

comp6_early_diffs_exp <- matrix(NaN, nrow = nrow(comp6_full_pool_early),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp6_full_pool_early)){
  temp_listener_pool <- comp6_speaker_pool_temp[which(comp6_speaker_pool_temp$Participant!=comp6_full_pool_early[i,2]),]
  listener_pool <- temp_listener_pool$EarlyResponse[which(temp_listener_pool$SpeakerEarly==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp6_early_diffs_exp[i,] <- temp_pulls-comp6_full_pool_early[i]
}

comp6_late_diffs_exp <- matrix(NaN, nrow = nrow(comp6_full_pool_late),ncol = nPulls)
temp_pulls <- rep(NaN, nPulls)
for(i in 1:nrow(comp6_full_pool_late)){
  temp_listener_pool <- comp6_speaker_pool_temp[which(comp6_speaker_pool_temp$Participant!=comp6_full_pool_late[i,2]),]
  listener_pool <- temp_listener_pool$LateResponse[which(temp_listener_pool$SpeakerLate==1)]
  for(j in 1:nPulls){
    pull_index <- round(runif(1,min = 1, max = length(listener_pool)),0)
    temp_pulls[j] <- listener_pool[pull_index]
  }
  comp6_late_diffs_exp[i,] <- temp_pulls-comp6_full_pool_late[i]
}

par(mfrow=c(2,1))
all_exp_diff <- c(comp1_late_diffs_exp,comp2_late_diffs_exp,comp3_late_diffs_exp,
                  comp4_late_diffs_exp,comp5_late_diffs_exp,comp6_late_diffs_exp,
                  comp1_early_diffs_exp,comp2_early_diffs_exp,comp3_early_diffs_exp,
                  comp4_early_diffs_exp,comp5_early_diffs_exp,comp6_early_diffs_exp)
hist(c(all_exp_diff),xlim=c(-1,1))
hist(all_nov_diff,xlim=c(-1,1))

par(mfrow=c(3,2))

hist(comp1_early_diffs_exp)
hist(comp2_early_diffs_exp)
hist(comp3_early_diffs_exp)
hist(comp4_early_diffs_exp)
hist(comp5_early_diffs_exp)
hist(comp6_early_diffs_exp)

hist(comp1_late_diffs_exp)
hist(comp2_late_diffs_exp)
hist(comp3_late_diffs_exp)
hist(comp4_late_diffs_exp)
hist(comp5_late_diffs_exp)
hist(comp6_late_diffs_exp)
