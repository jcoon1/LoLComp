# clears workspace:  
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load(".RData")

## load necessary packages

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)
library(fMultivar)

## load bootstrapping function

source("BootstrapMean95.R")

## load experienced participant data

Experienced_csv <- read_csv("Experienced_csv.csv",col_names=TRUE)

# Verhulst logistic growth model, by participant (hierarchical) ----------------------------------------------------

nReps = 12


nExp_aud = length(which(Experienced_csv$Expert_aud==1))/nReps
nNov_aud = length(which(Experienced_csv$Expert_aud==0))/nReps

y_exp <- Experienced_csv$TVJ[which(Experienced_csv$Expert_aud==1)]
y_nov <- Experienced_csv$TVJ[which(Experienced_csv$Expert_aud==0)]
p_exp <- Experienced_csv$Prior[which(Experienced_csv$Expert_aud==1)]
p_nov <- Experienced_csv$Prior[which(Experienced_csv$Expert_aud==0)]

data <- list("nReps","y_exp","y_nov","p_exp","p_nov","nExp_aud","nNov_aud") # to be passed on to JAGS
# myinits <- list(
#   list("muGrand" = runif(1,0,1), "v" = 3,"prior_v" = 3))

# parameters to be monitored:	
parameters <- c("alphaGrand_exp","alphaGrand_nov","betaGrand_exp","betaGrand_nov",
                "alpha_sigma_exp","alpha_sigma_nov","beta_sigma_exp","beta_sigma_nov",
                "alpha_exp","beta_exp","alpha_nov","beta_nov")
# parameters <- c("alphaGrand_exp","alphaGrand_nov","betaGrand_exp","betaGrand_nov",
#                 "alphaDelta_exp","alphaDelta_nov","betaDelta_exp","betaDelta_nov",
#                 "alpha_exp","beta_exp","alpha_nov","beta_nov")

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
# samples <- jags(data, 
#                 #inits=myinits, 
#                 parameters.to.save = parameters,
#                 model.file="hierarchical_aud_byPar.txt", n.chains=1,n.burnin = 1000,n.iter=4000, DIC=T)

#### For running multiple chains in parallel. Involves different jags call and initialization####
#Way this works is that when running chains in parallel, R basically initializes multiple jags calls with
#a single chain, so you build an initialization function to create new initializations
#each time R initiates the separate, parallel chain
# inits.jagsParallel=function(){
#   return(list(list("muGrand" = runif(1,0,1), "v" = 3, "prior_v" = 3)))
# }

#Now use jags.parallel to run multiple chains much quicker, adjust chains in n.chains

samples =jags.parallel(data,
                       #inits=inits.jagsParallel(), 
                       parameters.to.save = parameters,
                       model.file="hierarchical_aud_byPar.txt",n.chains=4,n.thin=1,n.burnin=1000,n.iter=4000)

# End parallel running section

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

samples_hierarchical <- samples

alphaGrand_exp_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$alphaGrand_exp
alphaGrand_nov_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$alphaGrand_nov
betaGrand_exp_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$betaGrand_exp
betaGrand_nov_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$betaGrand_nov

alphaDelta_exp_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$alphaDelta_exp
alphaDelta_nov_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$alphaDelta_nov
betaDelta_exp_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$betaDelta_exp
betaDelta_nov_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$betaDelta_nov

alpha_exp_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$alpha_exp
alpha_nov_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$alpha_nov
beta_exp_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$beta_exp
beta_nov_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$beta_nov

sigmaAlpha_nov_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$alpha_sigma_nov
                                
sigmaAlpha_exp_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$alpha_sigma_exp

sigmaBeta_nov_hierarchical    <- samples_hierarchical$BUGSoutput$sims.list$beta_sigma_nov
                                   
sigmaBeta_exp_hierarchical <- samples_hierarchical$BUGSoutput$sims.list$beta_sigma_exp

chainlength = nrow(alphaGrand_exp_hierarchical)/4

DIC_aud_hierarchical = samples_hierarchical$BUGSoutput$DIC


# Calculations for inference table (95% CI's) -------------------------------------


lo = 12000*.05
hi = 12000*.95
inference_table <- matrix(NaN, nrow = 8, ncol = 3)

inference_table[1,1] <- mean(alphaGrand_exp_hierarchical)
tempsort <- sort(alphaGrand_exp_hierarchical)
inference_table[1,2] <- tempsort[lo]
inference_table[1,3] <- tempsort[hi]

inference_table[2,1] <- mean(alphaGrand_nov_hierarchical)
tempsort <- sort(alphaGrand_nov_hierarchical)
inference_table[2,2] <- tempsort[lo]
inference_table[2,3] <- tempsort[hi]

inference_table[3,1] <- mean(betaGrand_exp_hierarchical)
tempsort <- sort(betaGrand_exp_hierarchical)
inference_table[3,2] <- tempsort[lo]
inference_table[3,3] <- tempsort[hi]

inference_table[4,1] <- mean(betaGrand_nov_hierarchical)
tempsort <- sort(betaGrand_nov_hierarchical)
inference_table[4,2] <- tempsort[lo]
inference_table[4,3] <- tempsort[hi]

inference_table[5,1] <- mean(sigmaAlpha_exp_hierarchical)
tempsort <- sort(sigmaAlpha_exp_hierarchical)
inference_table[5,2] <- tempsort[lo]
inference_table[5,3] <- tempsort[hi]

inference_table[6,1] <- mean(sigmaAlpha_nov_hierarchical)
tempsort <- sort(sigmaAlpha_nov_hierarchical)
inference_table[6,2] <- tempsort[lo]
inference_table[6,3] <- tempsort[hi]

inference_table[7,1] <- mean(sigmaBeta_exp_hierarchical)
tempsort <- sort(sigmaBeta_exp_hierarchical)
inference_table[7,2] <- tempsort[lo]
inference_table[7,3] <- tempsort[hi]

inference_table[8,1] <- mean(sigmaBeta_nov_hierarchical)
tempsort <- sort(sigmaBeta_nov_hierarchical)
inference_table[8,2] <- tempsort[lo]
inference_table[8,3] <- tempsort[hi]

mean(sigmaAlpha_exp_hierarchical)

# Set up binned overlay for main plot -------------------------------------

exp_data <- matrix(NaN,nrow = (nReps*length(p_exp)),ncol = 4)
for(i in 1:length(p_exp)){
  # for(j in 1:nReps){
    long_ticker = i
      # (i-1)*nReps+j
    exp_data[long_ticker,1] <- p_exp[long_ticker]
    exp_data[long_ticker,2] <- y_exp[long_ticker]
    if(exp_data[long_ticker,2]==2){
      exp_data[long_ticker,2]<- 0
    # }
  }
}

exp_overlay <- matrix(NaN,nrow = 10,ncol = 4)
exp_overlay[,1] <- seq(.035,.935,by = .1)

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

nov_data <- matrix(NaN,nrow = (nReps*length(p_nov)),ncol = 4)
for(i in 1:length(p_nov)){
  # for(j in 1:nReps){
    long_ticker = i
      # (i-1)*nReps+j
    nov_data[long_ticker,1] <- p_nov[long_ticker]
    nov_data[long_ticker,2] <- y_nov[long_ticker]
    if(nov_data[long_ticker,2]==2){
      nov_data[long_ticker,2]<- 0
    # }
  }
}

nov_overlay <- matrix(NaN,nrow = 10,ncol = 4)
nov_overlay[,1] <- seq(.065,.965,by = .1)

bin_size_nov = rep(NaN,10)

for(i in 1:nrow(nov_overlay)){
  tol = .0001
  lo = (i*.1)-.1
  hi = i * .1
  
  templo <- nov_data
  templo[,1] <- templo[,1]+tol
  temp_bin1 <- nov_data[which(templo[,1]>lo),1:2]
  
  temphi <- temp_bin1
  temphi[,1] <- temphi[,1]+tol
  if(i<10){
    temp_bin2 <- temp_bin1[which(temphi[,1]<hi),2]
  }
  else{temp_bin2 <- temp_bin1[,2]}
  bin_size_nov[i] = length(temp_bin2)
  
  nov_overlay[i,2] <- mean(temp_bin2)
  
  nSamples = 10
  nIterations = 1000
  
  if(i==1){nSamples = 2}
  
  bootRange <- BootstrapMean95(temp_bin2,nSamples,nIterations)
  
  nov_overlay[i,3:4] <- bootRange
}


# Calculate inexperienced interpretation for comparison -------------------

Inexperienced_replication <- read_csv("Inexperienced_replication_csv.csv",col_names=TRUE)

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

# Alex plot for hierarchical by participant --------------------------------------------------

par(mar = c(4, 4, 2, 2))

par(mfrow=c(1,1))

hist(alphaGrand_exp_hierarchical-alphaGrand_nov_hierarchical)
hist(betaGrand_exp_hierarchical-betaGrand_nov_hierarchical)

plot_index <- floor(runif(4000, min=1, max=12001))

q = seq(0,1,by = .01)
q_odds <- q/(1-q)
q_logodds <- log(q_odds)

temp_hierarchical <- q

temp_hierarchical <- 1/(1 + exp(-betaGrand_exp_hierarchical[plot_index[1]]*(q_logodds-alphaGrand_exp_hierarchical[plot_index[1]])))
plot(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)),frame.plot = FALSE,
     axes = FALSE,xlab = '',ylab = '')

xlabels = c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
xat = c(.1,.2,.3,.4,.5,.6,.7,.8,.9)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=0,cex.axis = 1.5,lwd = 2)

lines(x = c(0,1),y=c(0,0),lwd =2)

lines(x = c(0,1),y = c(.5,.5),lty = 4)
lines(x = c(.5,.5),y = c(0,1),lty = 4)


for(i in 2:length(alphaGrand_exp_hierarchical)){
  temp_hierarchical <- 1/(1 + exp(-betaGrand_exp_hierarchical[plot_index[i]]*(q_logodds-alphaGrand_exp_hierarchical[plot_index[i]])))
  lines(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)))
}

# plot_index <- floor(runif(4000, min=1, max=12001))

for(i in 1:length(alphaGrand_exp_hierarchical)){
  temp_hierarchical <- 1/(1 + exp(-betaGrand_nov_hierarchical[plot_index[i]]*(q_logodds-alphaGrand_nov_hierarchical[plot_index[i]])))
  lines(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)))
}

points(x=exp_overlay[,1],y=exp_overlay[,2],col = "darkblue",pch = 16, cex = 2.2)

for(i in 1:nrow(exp_overlay)){
  if(exp_overlay[i,3]-exp_overlay[i,4]!=0){
    if(exp_overlay[i,3]==0){
      exp_overlay[i,3] <- exp_overlay[i,3]+.005
    }
    arrows(x0 = exp_overlay[i,1], y0 = exp_overlay[i,3],
           x1 = exp_overlay[i,1], y1 = exp_overlay[i,4],code = 3,
           col = "darkblue",angle = 90, length = .1,lty = 2, lwd = 2)
  }
}

points(x=nov_overlay[,1],y=nov_overlay[,2], col = "firebrick",pch = 15,cex = 2.2,lwd = 2)

for(i in 1:nrow(nov_overlay)){
  if(nov_overlay[i,3]-nov_overlay[i,4]!=0){
    if(nov_overlay[i,3]==0){
      nov_overlay[i,3] <- nov_overlay[i,3]+.005
    }
    arrows(x0 = nov_overlay[i,1], y0 = nov_overlay[i,3],
           x1 = nov_overlay[i,1], y1 = nov_overlay[i,4],code = 3,
           col = "firebrick",angle = 90, length = .1,lty = 2,lwd = 2)
  }
}

points(x = SONA_descriptives[1,1],y = c(.05),pch = 18,cex = 2.2,lwd =2)
arrows(x0 = SONA_descriptives[2,1], y0 = c(.05),
       x1 = SONA_descriptives[3,1], y1 = c(.05),code = 3,
       lwd = 2, angle = 90, length = .1,lty=1)

legend(.005,.98,legend = c('Experienced Aud (Empirical)','Naive Aud (Empirical)',
                           'Experienced Aud (Inferred)','Naive Aud (Inferred)',
                           "Naive Interpretation"),
       pch = c(16,15,NA,NA,18),lty = c(NA,NA,1,1,NA),lwd = c(NA,2,2,2,2),pt.cex = c(2,2,NA,NA,2),
       col = c('darkblue','firebrick',"darkblue","firebrick","black"),
       cex = 1.2,bty = 'n',
       y.intersp = .7)



title(main=" ", line=.5, cex.main=1.8)
title(ylab="Endorsement Rate", line=1.6, cex.lab=1.8)
title(xlab="Estimated Applicability", line=1.8, cex.lab=1.8)

# Joint BF, hierarchical Verhulst ------------------------------------------------------
prior_alpha <- rnorm(12000, mean = 0, sd = 2)
prior_beta <- abs(rnorm(12000, mean = 5, sd = 20))
prior_deltaAlpha_hierarchical <- prior_alpha-rnorm(12000, mean = 0, sd = 2)
prior_deltaBeta_hierarchical <- prior_beta-abs(rnorm(12000, mean = 5, sd = 20))
post_deltaAlpha_hierarchical <- alphaGrand_exp_hierarchical - alphaGrand_nov_hierarchical
post_deltaBeta_hierarchical <- betaGrand_exp_hierarchical - betaGrand_nov_hierarchical

prior_jointDelta_hierarchical <- density2d(x = prior_deltaAlpha_hierarchical, y = prior_deltaBeta_hierarchical, n=200)
post_jointDelta_hierarchical <- density2d(x = post_deltaAlpha_hierarchical, y = post_deltaBeta_hierarchical, n=120)

par(mar = c(4, 4, 2, 2))

contour(prior_jointDelta_hierarchical,xlim = c(-1,1),ylim=c(-10,10),xlab = "Shift",ylab = "Scale")
contour(post_jointDelta_hierarchical,add=TRUE,col="green")

abline(h=0,lty = 3)
abline(v=0, lty = 3)

point0_index_prior_hierarchical <- c(which(abs(prior_jointDelta_hierarchical$x) == min(abs(prior_jointDelta_hierarchical$x))),
                                 which(abs(prior_jointDelta_hierarchical$y) == min(abs(prior_jointDelta_hierarchical$y))))

prior_density_jointDelta_hierarchical <- prior_jointDelta_hierarchical$z[point0_index_prior_hierarchical[1],point0_index_prior_hierarchical[2]]

point0_index_post_hierarchical <- c(which(abs(post_jointDelta_hierarchical$x) == min(abs(post_jointDelta_hierarchical$x))),
                                     which(abs(post_jointDelta_hierarchical$y) == min(abs(post_jointDelta_hierarchical$y))))

post_density_jointDelta_hierarchical <- post_jointDelta_hierarchical$z[point0_index_post_hierarchical[1],point0_index_post_hierarchical[2]]

SavDick01_jointDelta_hierarchical <- post_density_jointDelta_hierarchical/prior_density_jointDelta_hierarchical

SavDick01_jointDelta_hierarchical


# Joint BF sanity check ---------------------------------------------------
jointprior <- matrix(NaN,nrow = length(prior_deltaAlpha_hierarchical),ncol = 2)
jointprior[,1] <- prior_deltaAlpha_hierarchical
jointprior[,2] <- prior_deltaBeta_hierarchical
jointpost <- cbind(post_deltaAlpha_hierarchical,post_deltaBeta_hierarchical)

radius_check = matrix(NaN,nrow=3,ncol=10)

radius_check[1,] <- c(15,10,1,.5,.4,.3,.2,.15,.12,.1)

for(j in 1:ncol(radius_check)){
  ticker_post = 0
  ticker_prior = 0
  for(i in 1:nrow(jointprior)){
    radius_comp_post = sqrt(((10*jointpost[i,1])^2) + (jointpost[i,2]^2))
    if(radius_check[1,j]>radius_comp_post){
      ticker_post = ticker_post+1
    }
    radius_comp_prior = sqrt(((10*jointprior[i,1])^2) + (jointprior[i,2]^2))
    if(radius_check[1,j]>radius_comp_prior){
      ticker_prior = ticker_prior+1
    }
  }
  radius_check[2,j] <- ticker_prior
  radius_check[3,j] <- ticker_post
}
# Bayes Factor Analysis (single variable) ------------------------------------------------

prior_density_deltaAlpha = dlogspline(0,logspline(prior_deltaAlpha_hierarchical))
post_density_deltaAlpha <- dlogspline(0,logspline(post_deltaAlpha_hierarchical))
SavDick01_deltaAlpha_hierarchical <- post_density_deltaAlpha/prior_density_deltaAlpha

prior_density_deltaBeta = dlogspline(0,logspline(prior_deltaBeta_hierarchical))
post_density_deltaBeta <- dlogspline(0,logspline(post_deltaBeta_hierarchical))
SavDick01_deltaBeta_hierarchical <- post_density_deltaBeta/prior_density_deltaBeta

deltaAlpha_sd  = sd(prior_deltaAlpha_hierarchical)
prior_prob_deltaAlpha_negative = pnorm(0,0,deltaAlpha_sd)
prior_odds_deltaAlpha_negative = prior_prob_deltaAlpha_negative/(1-prior_prob_deltaAlpha_negative)
posterior_prob_deltaAlpha_negative <- mean(post_deltaAlpha_hierarchical<0)
posterior_odds_deltaAlpha_negative <- posterior_prob_deltaAlpha_negative/(1-posterior_prob_deltaAlpha_negative)
BF10_deltaAlpha_negative <- posterior_odds_deltaAlpha_negative/prior_odds_deltaAlpha_negative

deltaBeta_sd  = sd(prior_deltaBeta_hierarchical)
prior_prob_deltaBeta_negative = pnorm(0,0,deltaBeta_sd)
prior_odds_deltaBeta_negative = prior_prob_deltaBeta_negative/(1-prior_prob_deltaBeta_negative)
posterior_prob_deltaBeta_negative <- mean(post_deltaBeta_hierarchical<0)
posterior_odds_deltaBeta_negative <- posterior_prob_deltaBeta_negative/(1-posterior_prob_deltaBeta_negative)
BF10_deltaBeta_negative <- posterior_odds_deltaBeta_negative/prior_odds_deltaBeta_negative

# # Linking function by participant -----------------------------------------
#  #Only uncomment this if you really need it (R freaks out)
# par(mfrow=c(5,4))
# par(mar = c(3, 4, 2, 2))
# 
# for(j in 1:20){
#   hi = j*nReps
#   lo= hi - (nReps-1)
#   temp_hierarchical <- q
# 
#   temp_hierarchical <- 1/(1 + exp(-beta_exp_hierarchical[1,j]*(q_logodds-alpha_exp_hierarchical[1,j])))
#   plot(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)),xlab = '',ylab = '',main = j)
#   lines(x = c(0,1),y = c(.5,.5),lty = 4)
#   lines(x = c(.5,.5),y = c(0,1),lty = 4)
#   points(x=p_exp[lo:hi],y=y_exp[lo:hi])
# 
#   for(i in 2:length(alpha_exp_hierarchical[,j])){
#     temp_hierarchical <- 1/(1 + exp(-beta_exp_hierarchical[i,j]*(q_logodds-alpha_exp_hierarchical[i,j])))
#     lines(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)))
#   }
# 
# }

# for(j in 21:nExp_aud){
#   hi = j*nReps
#   lo= hi - (nReps-1)
#   temp_hierarchical <- q
# 
#   temp_hierarchical <- 1/(1 + exp(-beta_exp_hierarchical[1,j]*(q_logodds-alpha_exp_hierarchical[1,j])))
#   plot(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)),xlab = '',ylab = '',main = j)
#   lines(x = c(0,1),y = c(.5,.5),lty = 4)
#   lines(x = c(.5,.5),y = c(0,1),lty = 4)
#   points(x=p_exp[lo:hi],y=y_exp[lo:hi])
# 
#   for(i in 2:length(alpha_exp_hierarchical[,j])){
#     temp_hierarchical <- 1/(1 + exp(-beta_exp_hierarchical[i,j]*(q_logodds-alpha_exp_hierarchical[i,j])))
#     lines(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)))
#   }
# 
# }
# 
# for(j in 1:20){
#   hi = j*nReps
#   lo= hi - (nReps-1)
#   temp_hierarchical <- q
# 
#   temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[1,j]*(q_logodds-alpha_nov_hierarchical[1,j])))
#   plot(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)),xlab = '',ylab = '',main = j)
#   lines(x = c(0,1),y = c(.5,.5),lty = 4)
#   lines(x = c(.5,.5),y = c(0,1),lty = 4)
#   points(x=p_nov[lo:hi],y=y_nov[lo:hi])
# 
#   for(i in 2:length(alpha_nov_hierarchical[,j])){
#     temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[i,j]*(q_logodds-alpha_nov_hierarchical[i,j])))
#     lines(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)))
#   }
# 
# }
# 
# for(j in 21:nNov_aud){
#   hi = j*nReps
#   lo= hi - (nReps-1)
#   temp_hierarchical <- q
# 
#   temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[1,j]*(q_logodds-alpha_nov_hierarchical[1,j])))
#   plot(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)),xlab = '',ylab = '',main = j)
#   lines(x = c(0,1),y = c(.5,.5),lty = 4)
#   lines(x = c(.5,.5),y = c(0,1),lty = 4)
#   points(x=p_nov[lo:hi],y=y_nov[lo:hi])
# 
#   for(i in 2:length(alpha_nov_hierarchical[,j])){
#     temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[i,j]*(q_logodds-alpha_nov_hierarchical[i,j])))
#     lines(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)))
#   }
# 
# }
# 
# 
# # Graph by type of participant --------------------------------------------
# 
# tight_exp_index <- c(1,8,11,14,20,23,25,28,30,38)
# tight_nov_index <- c(3,6,7,10,11,13,16,17,21,23,26,32)
# 
# par(mfrow=c(4,4))
# par(mar = c(3, 4, 2, 2))
# 
# for(j in 1:length(tight_exp_index)){
#   hi = tight_exp_index[j]*nReps
#   lo= hi - (nReps-1)
#   temp_hierarchical <- q
#   
#   temp_hierarchical <- 1/(1 + exp(-beta_exp_hierarchical[1,tight_exp_index[j]]*(q_logodds-alpha_exp_hierarchical[1,tight_exp_index[j]])))
#   plot(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)),xlab = '',ylab = '',main = tight_exp_index[j])
#   lines(x = c(0,1),y = c(.5,.5),lty = 4)
#   lines(x = c(.5,.5),y = c(0,1),lty = 4)
#   points(x=p_exp[lo:hi],y=y_exp[lo:hi])
#   
#   for(i in 2:length(alpha_exp_hierarchical[,tight_exp_index[j]])){
#     temp_hierarchical <- 1/(1 + exp(-beta_exp_hierarchical[i,tight_exp_index[j]]*(q_logodds-alpha_exp_hierarchical[i,tight_exp_index[j]])))
#     lines(q,temp_hierarchical,type = 'l',col=c(rgb(0,0,.5451,0.002)))
#   }
# }
# 
# for(j in 1:6){
#   hi = tight_nov_index[j]*nReps
#   lo= hi - (nReps-1)
#   temp_hierarchical <- q
#   
#   temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[1,tight_nov_index[j]]*(q_logodds-alpha_nov_hierarchical[1,tight_nov_index[j]])))
#   plot(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)),xlab = '',ylab = '',main = tight_nov_index[j])
#   lines(x = c(0,1),y = c(.5,.5),lty = 4)
#   lines(x = c(.5,.5),y = c(0,1),lty = 4)
#   points(x=p_nov[lo:hi],y=y_nov[lo:hi])
#   
#   for(i in 2:length(alpha_nov_hierarchical[,tight_nov_index[j]])){
#     temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[i,tight_nov_index[j]]*(q_logodds-alpha_nov_hierarchical[i,tight_nov_index[j]])))
#     lines(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)))
#   }
# }
# 
# for(j in 7:length(tight_nov_index)){
#   hi = tight_nov_index[j]*nReps
#   lo= hi - (nReps-1)
#   temp_hierarchical <- q
#   
#   temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[1,tight_nov_index[j]]*(q_logodds-alpha_nov_hierarchical[1,tight_nov_index[j]])))
#   plot(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)),xlab = '',ylab = '',main = tight_nov_index[j])
#   lines(x = c(0,1),y = c(.5,.5),lty = 4)
#   lines(x = c(.5,.5),y = c(0,1),lty = 4)
#   points(x=p_nov[lo:hi],y=y_nov[lo:hi])
#   
#   for(i in 2:length(alpha_nov_hierarchical[,tight_nov_index[j]])){
#     temp_hierarchical <- 1/(1 + exp(-beta_nov_hierarchical[i,tight_nov_index[j]]*(q_logodds-alpha_nov_hierarchical[i,tight_nov_index[j]])))
#     lines(q,temp_hierarchical,type = 'l',col=c(rgb(.698,.1333,.1333,0.002)))
#   }
# }

# Illustrative graph --------------------------------------------------

par(mfrow=c(1,1))

par(mar = c(3.5, 3, 1, 2))

temp_hierarchical <- q

plot(q,q,type = 'l',col=c(rgb(.698,.1333,.1333,0.0002)),frame.plot = FALSE,
     axes = FALSE,xlab = '',ylab = '',xlim=c(0,1),ylim=c(0,1))

xlabels = c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
xat = c(.1,.2,.3,.4,.5,.6,.7,.8,.9)

axis(side=1,pos=0,cex.axis = 1.4,lwd = 2,labels = xlabels,at = xat)
axis(side = 2, pos=0,cex.axis = 1.5,lwd = 2)

title(main=" ", line=.5, cex.main=1.8)
title(ylab="Endorsement Rate", line=1.6, cex.lab=1.8)
title(xlab="Estimated Applicability", line=1.8, cex.lab=1.8)

lines(x = c(0,1),y=c(0,0),lwd =2)


example_threshold <- .5
example_alpha <- log(example_threshold/(1-example_threshold))
example_beta <- 45
temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
lines(q,temp_hierarchical,type = 'l',lwd = 2,col="tomato",lty = 2)

example_threshold <- .5
example_alpha <- log(example_threshold/(1-example_threshold))
example_beta <- 1
temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
lines(q,temp_hierarchical,type = 'l',lwd = 2,col="tomato",lty = 2)

example_threshold <- .6
example_alpha <- log(example_threshold/(1-example_threshold))
example_beta <- 5
temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
lines(q,temp_hierarchical,type = 'l',lwd = 2,col="springgreen4",lty = 6)

example_threshold <- .4
example_alpha <- log(example_threshold/(1-example_threshold))
example_beta <- 5
temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
lines(q,temp_hierarchical,type = 'l',lwd = 2,col="springgreen4",lty = 6)

lines(x = c(0,1),y = c(.5,.5),lty = 4)
lines(x = c(.5,.5),y = c(0,1),lty = 4)

legend(.005,.98,legend = c(expression(paste(alpha, " varies; ", beta, " = ", 5)),
                           expression(paste(alpha, " = log-odds(", .5,"); ",beta, " varies"))),
       pch = c(NA,NA),lty = c(6,2),lwd = c(2),
       col = c("springgreen4","tomato"),
       cex = 1.5,bty = 'n',
       y.intersp = .7)

# legend(.005,.98,legend = expression(paste(alpha, " varies; ", beta, " = ", 5)),
#        pch = c(NA),lty = c(6),lwd = c(2),
#        col = c("springgreen4"),
#        cex = 1.5,bty = 'n')