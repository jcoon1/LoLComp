# clears workspace:  
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("splitsig.RData")

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)

front_end <- runif(10000,min = 0, max = 1)

for(i in 1:1000){}