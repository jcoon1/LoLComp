# clears workspace:  
rm(list=ls()) 

setwd('C:/Users/Jeff/Documents/R/LoLComp')

load(".RData")

library(tidyverse)
library(BayesFactor)
library(corrplot)
library(polspline)
library(R2jags)


# Load expert data --------------------------------------------------------

OnlineLoad <- read_csv("LoLCompResults_collection_forupload2.csv",col_names=TRUE)

nExperts = nrow(OnlineLoad)/60
nItems = 12

Compositions = matrix(NaN,nrow = 12,ncol = 8)

Compositions[1,1:5] <- c(40,62,24,15,12)
Compositions[2,1:5] <- c(59,83,89,46,142)
Compositions[3,1:5] <- c(85,25,81,65,11)
Compositions[4,1:5] <- c(133,92,140,46,130)
Compositions[5,1:5] <- c(29,22,98,47,78)
Compositions[6,1:5] <- c(54,44,21,73,66)

Compositions[7,1:5] <- c(68,36,84,124,11)
Compositions[8,1:5] <- c(85,143,61,7,66)
Compositions[9,1:5] <- c(33,71,130,46,76)
Compositions[10,1:5] <- c(87,62,112,49,76)
Compositions[11,1:5] <- c(118,36,61,73,100)
Compositions[12,1:5] <- c(16,62,86,65,11)

CompCheck <- function(Comp,CompositionTable) {
  for(i in 1:nrow(CompositionTable))
    if(CompositionTable[i,1]==Comp[1] && CompositionTable[i,2]==Comp[2] &&
       CompositionTable[i,3]==Comp[3] && CompositionTable[i,4]==Comp[4] &&
       CompositionTable[i,5]==Comp[5]){
      CompID = i
    }
  return(CompID)
}

for(i in 1:12){
  Borderline = rep(c(1,1,2,2,3,3),2)
  Compositions[i,6] <- i
  Compositions[i,7] <- round(i/2+.2)
  Compositions[i,8] <- Borderline[i]
}

tempComp = rep(NaN,5)

totalExpert = nItems*nExperts

Comp = rep(NaN,totalExpert)
CompType = rep(NaN,totalExpert)
EarlyResponse = rep(NaN,totalExpert)
LateResponse = rep(NaN,totalExpert)
SpeakerEarly = rep(NaN,totalExpert)
SpeakerLate = rep(NaN,totalExpert)
ParticipantID = rep(NaN,totalExpert)
EarlyTVJ = rep(NaN,totalExpert)
LateTVJ = rep(NaN,totalExpert)
EarlyPrior = rep(NaN,totalExpert)
LatePrior = rep(NaN,totalExpert)
GenPrior = rep(NaN,totalExpert)
EarlyShift = rep(NaN,totalExpert)
LateShift = rep(NaN,totalExpert)
CompBorderline = rep(NaN,totalExpert)
Skill = rep(NaN,totalExpert)
ChampFam_total = rep(NaN,totalExpert)
ChampFam = rep(NaN,totalExpert)


for(i in 1:totalExpert){
  EarlyContinue = (i*5)-1
  LateContinue = i*5
  
  EarlyResponse[i] <- OnlineLoad$speaker_excel_early[EarlyContinue]
  LateResponse[i] <- OnlineLoad$speaker_excel_late[LateContinue]
  
  tempComp <- c(OnlineLoad$blueTopIndex[EarlyContinue],OnlineLoad$blueJungleIndex[EarlyContinue],
                OnlineLoad$blueMidIndex[EarlyContinue],OnlineLoad$blueBotIndex[EarlyContinue],
                OnlineLoad$blueSupportIndex[EarlyContinue])
  
  CompID <- CompCheck(tempComp,Compositions)
  
  Comp[i] <- CompID
  CompType[i] <- Compositions[CompID,7]
  
  SpeakerEarly[i] <- OnlineLoad$speaker_early_index[LateContinue]
  SpeakerLate[i] <- OnlineLoad$speaker_late_index[LateContinue]
  
  ParticipantID[i] <- OnlineLoad$participant_id[LateContinue]
  
  TVJContinue = (i*5)-3
  EarlyTVJ[i] <- OnlineLoad$TVJearly[TVJContinue]
  LateTVJ[i] <- OnlineLoad$TVJlate[TVJContinue]
  
  PriorsContinue = (i*5)-2
  EarlyPrior[i] <- OnlineLoad$blue_excel_early[PriorsContinue]
  LatePrior[i] <- OnlineLoad$blue_excel_late[PriorsContinue]
  GenPrior[i] <- OnlineLoad$blue_excel_general[PriorsContinue]
  
  EarlyShift[i] <- EarlyResponse[i] - EarlyPrior[i]
  LateShift[i] <- LateResponse[i] - LatePrior[i]
  
  CompBorderline[i] <- Compositions[CompID,8]
  
  Skill[i] <- OnlineLoad$league_level[EarlyContinue]
  
  ChampFam[i] <- OnlineLoad$ChampFam[EarlyContinue]
  
  ChampFam_total[i] <- OnlineLoad$champfam_total[EarlyContinue]
  
}

OnlineData <- data.frame(Comp, CompType, EarlyResponse, LateResponse, SpeakerEarly,
                         SpeakerLate, ParticipantID, EarlyTVJ,LateTVJ,EarlyPrior,
                         LatePrior,GenPrior, EarlyShift, LateShift, CompBorderline,
                         Skill, ChampFam, ChampFam_total)

OnlineData$Participant <- rep(1:nExperts,each=12)
# Load full SONA data (with prior questions) ------------------------------

SONALoad <- read_csv("LoLCompResults_SONA_contaminant.csv",col_names=TRUE)

nNovices = nrow(SONALoad)/60

totalNovice = nItems*nNovices

Comp = rep(NaN,totalNovice)
CompType = rep(NaN,totalNovice)
EarlyResponse = rep(NaN,totalNovice)
LateResponse = rep(NaN,totalNovice)
SpeakerEarly = rep(NaN,totalNovice)
SpeakerLate = rep(NaN,totalNovice)
ParticipantID = rep(NaN,totalNovice)
EarlyTVJ = rep(NaN,totalNovice)
LateTVJ = rep(NaN,totalNovice)
EarlyPrior = rep(NaN,totalNovice)
LatePrior = rep(NaN,totalNovice)
GenPrior = rep(NaN,totalNovice)
EarlyShift = rep(NaN,totalNovice)
LateShift = rep(NaN,totalNovice)
CompBorderline = rep(NaN,totalNovice)

for(i in 1:totalNovice){
  EarlyContinue = (i*5)-1
  LateContinue = i*5
  EarlyResponse[i] <- SONALoad$speaker_excel_early[EarlyContinue]
  LateResponse[i] <- SONALoad$speaker_excel_late[LateContinue]
  
  tempComp <- c(SONALoad$blueTopIndex[EarlyContinue],SONALoad$blueJungleIndex[EarlyContinue],
                SONALoad$blueMidIndex[EarlyContinue],SONALoad$blueBotIndex[EarlyContinue],
                SONALoad$blueSupportIndex[EarlyContinue])
  
  CompID <- CompCheck(tempComp,Compositions)
  
  Comp[i] <- CompID
  CompType[i] <- Compositions[CompID,7]
  
  SpeakerEarly[i] <- SONALoad$speaker_early_index[LateContinue]
  SpeakerLate[i] <- SONALoad$speaker_late_index[LateContinue]
  
  ParticipantID[i] <- SONALoad$participant_id[LateContinue]
  
  TVJContinue = (i*5)-3
  EarlyTVJ[i] <- SONALoad$TVJearly[TVJContinue]
  LateTVJ[i] <- SONALoad$TVJlate[TVJContinue]
  
  PriorsContinue = (i*5)-2
  EarlyPrior[i] <- SONALoad$blue_excel_early[PriorsContinue]
  LatePrior[i] <- SONALoad$blue_excel_late[PriorsContinue]
  GenPrior[i] <- SONALoad$blue_excel_general[PriorsContinue]
  
  EarlyShift[i] <- EarlyResponse[i] - EarlyPrior[i]
  LateShift[i] <- LateResponse[i] - LatePrior
  
  CompBorderline[i] <- Compositions[CompID,8]
  
}

SONAData <- data.frame(Comp, CompType, EarlyResponse, LateResponse, SpeakerEarly,
                       SpeakerLate, ParticipantID, EarlyTVJ, LateTVJ, EarlyPrior,
                       LatePrior, GenPrior, EarlyShift, LateShift, CompBorderline)

# Start preliminary graphing ----------------------------------------------


## Histogram of all responses (not split by speaker yes/no)
par(mfrow=c(2,2))
hist(SONAData$EarlyResponse, ylim = c(0,150), main = "Early Novice",breaks = 10)
abline(v= mean(SONAData$EarlyResponse), col='#00CC00')
abline(v= mean(OnlineData$EarlyResponse), col='#7B68EE')
hist(OnlineData$EarlyResponse, ylim = c(0,150), main = "Early Expert",breaks = 10)
hist(SONAData$LateResponse, ylim = c(0,150), main = "Late Novice",breaks = 10)
abline(v= mean(SONAData$LateResponse), col='#00CC00')
abline(v= mean(OnlineData$LateResponse), col='#7B68EE')
hist(OnlineData$LateResponse, ylim = c(0,150), main = "Late Expert",breaks = 10)


## Histogram of all responses by composition (not split by speaker yes/no)
EarlyCompTypeTitles= c('Mega Early Novice','Borderline Early Novice','Not Early Novice',
                       'Mega Early Expert','Borderline Early Expert','Not Early Expert')

LateCompTypeTitles= c('Mega Late Novice','Borderline Late Novice','Not Late Novice',
                      'Mega Late Expert','Borderline Late Expert','Not Late Expert')

for(j in 1:2){
  par(mfrow=c(3,2))
  for(i in 1:3){
    if(j == 1){
      hist(SONAData$EarlyResponse[which(SONAData$CompType==i)],main = EarlyCompTypeTitles[i],xlab = "",
           ylim = c(0,30), xlim = c(0,1),breaks = 10
      )
      abline(v = mean(SONAData$EarlyResponse[which(SONAData$CompType==i)]),col = '#00CC00')
      abline(v = mean(OnlineData$EarlyResponse[which(OnlineData$CompType==i)]),col='#7B68EE')
      hist(OnlineData$EarlyResponse[which(OnlineData$CompType==i)],main = EarlyCompTypeTitles[i+3],xlab = "",
           ylim = c(0,30), xlim = c(0,1),breaks = 10
      )
    }
    if(j == 2){
      hist(SONAData$LateResponse[which(SONAData$CompType==i+3)],main = LateCompTypeTitles[i],xlab = "",
           ylim = c(0,30), xlim = c(0,1),breaks = 10
      )
      abline(v = mean(SONAData$LateResponse[which(SONAData$CompType==i+3)]),col='#00CC00')
      abline(v = mean(OnlineData$LateResponse[which(OnlineData$CompType==i+3)]),col = '#7B68EE')
      hist(OnlineData$LateResponse[which(OnlineData$CompType==i+3)],main = LateCompTypeTitles[i+3],xlab = "",
           ylim = c(0,30), xlim = c(0,1),breaks = 10
      )
      
    }
  }
}

## Histogram of all responses split by speaker yes/no

EarlySpeakerTitles = c('No Early Novice','Yes Early Novice','No Early Expert','Yes Early Expert')
LateSpeakerTitles = c('No Late Novice','Yes Late Novice','No Late Expert','Yes Late Expert')

for(j in 1:2){
  par(mfrow=c(2,2))
  for(i in 1:2){
    if(j==1){
      hist(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==i-1)],main = EarlySpeakerTitles[i],xlab = "",
           ylim = c(0,110), xlim = c(0,1),breaks = 10
      )
      abline(v= mean(SONAData$EarlyResponse[which(SONAData$SpeakerEarly==i-1)]),col='#00CC00')
      abline(v= mean(OnlineData$EarlyResponse[which(OnlineData$SpeakerEarly==i-1)]),col='#7B68EE')
      hist(OnlineData$EarlyResponse[which(OnlineData$SpeakerEarly==i-1)],main = EarlySpeakerTitles[i+2],xlab = "",
           ylim = c(0,110), xlim = c(0,1),breaks = 10
      )
    }
    if(j==2){
      hist(SONAData$LateResponse[which(SONAData$SpeakerLate==i-1)],main = LateSpeakerTitles[i],xlab = "",
           ylim = c(0,110), xlim = c(0,1),breaks = 10
      )
      abline(v= mean(SONAData$LateResponse[which(SONAData$SpeakerLate==i-1)]),col='#00CC00')
      abline(v= mean(OnlineData$LateResponse[which(OnlineData$SpeakerLate==i-1)]),col='#7B68EE')
      hist(OnlineData$LateResponse[which(OnlineData$SpeakerLate==i-1)],main = LateSpeakerTitles[i+2],xlab = "",
           ylim = c(0,110), xlim = c(0,1),breaks = 10
      )
    }
  }
}

## Histogram of responses on non-referenced trait

EarlySpeakerTitles = c('No Early Novice (Oppo)','Yes Early Novice (Oppo)','No Early Expert (Oppo)','Yes Early Expert (Oppo)')
LateSpeakerTitles = c('No Late Novice (Oppo)','Yes Late Novice (Oppo)','No Late Expert (Oppo)','Yes Late Expert (Oppo)')

for(j in 1:2){
  par(mfrow=c(2,2))
  for(i in 1:2){
    if(j==1){
      hist(SONAData$LateResponse[which(SONAData$SpeakerEarly==i-1)],main = EarlySpeakerTitles[i],xlab = "",
           ylim = c(0,80), xlim = c(0,1),breaks = 10
      )
      hist(OnlineData$LateResponse[which(OnlineData$SpeakerEarly==i-1)],main = EarlySpeakerTitles[i+2],xlab = "",
           ylim = c(0,80), xlim = c(0,1),breaks = 10
      )
    }
    if(j==2){
      hist(SONAData$EarlyResponse[which(SONAData$SpeakerLate==i-1)],main = LateSpeakerTitles[i],xlab = "",
           ylim = c(0,80), xlim = c(0,1),breaks = 10
      )
      hist(OnlineData$EarlyResponse[which(OnlineData$SpeakerLate==i-1)],main = LateSpeakerTitles[i+2],xlab = "",
           ylim = c(0,80), xlim = c(0,1),breaks = 10
      )
    }
  }
}

## Shift from prior to interpretation (not split by speaker yes/no)
par(mfrow=c(2,2))
hist(SONAData$EarlyShift, ylim = c(0,250), main = "Early Shift Novice", xlim = c(-1,1),breaks = 10)
hist(OnlineData$EarlyShift, ylim = c(0,250), main = "Early Shift Expert", xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8))
hist(SONAData$LateShift, ylim = c(0,250), main = "Late Shift Novice", xlim = c(-1,1),breaks = 10)
hist(OnlineData$LateShift, ylim = c(0,250), main = "Late Shift Expert", xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8))

## Shift from prior to interpretation split by speaker yes/no
EarlySpeakerTitles = c('No Early Shift Novice','Yes Early Shift Novice','No Early Shift Expert','Yes Early Shift Expert')
LateSpeakerTitles = c('No Late Shift Novice','Yes Late Shift Novice','No Late Shift Expert','Yes Late Shift Expert')

for(j in 1:2){
  par(mfrow=c(2,2))
  for(i in 1:2){
    if(j==1){
      hist(SONAData$EarlyShift[which(SONAData$SpeakerEarly==i-1)],main = EarlySpeakerTitles[i],xlab = "",
           ylim = c(0,150), xlim = c(-1,1),breaks = 10
      )
      abline(v = mean(SONAData$EarlyShift[which(SONAData$SpeakerEarly==i-1)]), col = '#00CC00')
      abline(v = mean(OnlineData$EarlyShift[which(OnlineData$SpeakerEarly==i-1)]), col = '#7B68EE')
      hist(OnlineData$EarlyShift[which(OnlineData$SpeakerEarly==i-1)],main = EarlySpeakerTitles[i+2],xlab = "",
           ylim = c(0,150), xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8)
      )
    }
    if(j==2){
      hist(SONAData$LateShift[which(SONAData$SpeakerLate==i-1)],main = LateSpeakerTitles[i],xlab = "",
           ylim = c(0,150), xlim = c(-1,1),breaks = 10
      )
      abline(v = mean(SONAData$EarlyShift[which(SONAData$SpeakerEarly==i-1)]), col = '#00CC00')
      abline(v = mean(OnlineData$LateShift[which(OnlineData$SpeakerLate==i-1)]), col = '#7B68EE')
      hist(OnlineData$LateShift[which(OnlineData$SpeakerLate==i-1)],main = LateSpeakerTitles[i+2],xlab = "",
           ylim = c(0,150), xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8)
      )
    }
  }
}

## Shift from prior to interpretation contingent on whether the participant agrees with speaker decision
EarlySpeakerTitles = c('Agree Early Shift Novice','Disagree Early Shift Novice','Agree Early Shift Expert','Disagree Early Shift Expert')
LateSpeakerTitles = c('Agree Late Shift Novice','Disagree Late Shift Novice','Agree Late Shift Expert','Disagree Late Shift Expert')

for(j in 1:2){
  par(mfrow=c(2,2))
  if(j==1){
    hist(SONAData$EarlyShift[which(SONAData$SpeakerEarly==SONAData$EarlyTVJ)],main = EarlySpeakerTitles[1],xlab = "",
         ylim = c(0,170), xlim = c(-1,1),breaks = 10
    )
    hist(OnlineData$EarlyShift[which(OnlineData$SpeakerEarly==OnlineData$EarlyTVJ)],main = EarlySpeakerTitles[3],xlab = "",
         ylim = c(0,170), xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8)
    )     
    hist(SONAData$EarlyShift[which(SONAData$SpeakerEarly!=SONAData$EarlyTVJ)],main = EarlySpeakerTitles[2],xlab = "",
         ylim = c(0,170), xlim = c(-1,1),breaks = 10
    )
    hist(OnlineData$EarlyShift[which(OnlineData$SpeakerEarly!=OnlineData$EarlyTVJ)],main = EarlySpeakerTitles[4],xlab = "",
         ylim = c(0,170), xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8)
    )
  }
  if(j==2){
    hist(SONAData$LateShift[which(SONAData$SpeakerLate==SONAData$LateTVJ)],main = LateSpeakerTitles[1],xlab = "",
         ylim = c(0,170), xlim = c(-1,1)
    )
    hist(OnlineData$LateShift[which(OnlineData$SpeakerLate==OnlineData$LateTVJ)],main = LateSpeakerTitles[3],xlab = "",
         ylim = c(0,170), xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8)
    )     
    hist(SONAData$LateShift[which(SONAData$SpeakerLate!=SONAData$LateTVJ)],main = LateSpeakerTitles[2],xlab = "",
         ylim = c(0,170), xlim = c(-1,1)
    )
    hist(OnlineData$LateShift[which(OnlineData$SpeakerLate!=OnlineData$LateTVJ)],main = LateSpeakerTitles[4],xlab = "",
         ylim = c(0,170), xlim = c(-1,1),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8)
    )
  }
}

## TVJ response rate
novice_EarlyTVJRate = mean(SONAData$EarlyTVJ)
novice_LateTVJRate = mean(SONAData$LateTVJ)
expert_EarlyTVJRate = mean(OnlineData$EarlyTVJ)
expert_LateTVJRate = mean(OnlineData$LateTVJ)

par(mfrow=c(2,2))

TVJBarGraph6 <- c(mean(SONAData$EarlyTVJ[which(SONAData$CompType==1)]),
                  mean(SONAData$EarlyTVJ[which(SONAData$CompType==2)]),
                  mean(SONAData$EarlyTVJ[which(SONAData$CompType==3)]))

barplot(TVJBarGraph6,main="Novice Early TVJ")

TVJBarGraph6 <- c(mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==1)]),
                  mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==2)]),
                  mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==3)]))

barplot(TVJBarGraph6,main="Expert Early TVJ")

TVJBarGraph6 <- c(mean(SONAData$LateTVJ[which(SONAData$CompType==4)]),
                  mean(SONAData$LateTVJ[which(SONAData$CompType==5)]),
                  mean(SONAData$LateTVJ[which(SONAData$CompType==6)]))

barplot(TVJBarGraph6,main="Novice Late TVJ")

TVJBarGraph6 <- c(mean(OnlineData$LateTVJ[which(OnlineData$CompType==4)]),
                  mean(OnlineData$LateTVJ[which(OnlineData$CompType==5)]),
                  mean(OnlineData$LateTVJ[which(OnlineData$CompType==6)]))

barplot(TVJBarGraph6,main="Expert Late TVJ")

## Responses split by comp type and expert speaker yes/no
EarlyCompTypeTitles= c('Mega Early Yes Novice','Borderline Early Yes Novice','Not Early Yes Novice',
                       'Mega Early Yes Expert','Borderline Early Yes Expert','Not Early Yes Expert')

LateCompTypeTitles= c('Mega Late Yes Novice','Borderline Late Yes Novice','Not Late Yes Novice',
                      'Mega Late Yes Expert','Borderline Late Yes Expert','Not Late Yes Expert')

for(j in 1:2){
  par(mfrow=c(3,2))
  for(i in 1:3){
    if(j == 1){
      hist(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==1)]==i)],main = EarlyCompTypeTitles[i],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
      abline(v=mean(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==1)]==i)]))
      hist(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]==i)],main = EarlyCompTypeTitles[i+3],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
      abline(v=mean(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]==i)]))
    }
    if(j == 2){
      hist(SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==1)]==i+3)],main = LateCompTypeTitles[i],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
      hist(OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==1)]==i+3)],main = LateCompTypeTitles[i+3],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
    }
  }
}

EarlyCompTypeTitles= c('Mega Early No Novice','Borderline Early No Novice','Not Early No Novice',
                       'Mega Early No Expert','Borderline Early No Expert','Not Early No Expert')

LateCompTypeTitles= c('Mega Late No Novice','Borderline Late No Novice','Not Late No Novice',
                      'Mega Late No Expert','Borderline Late No Expert','Not Late No Expert')

for(j in 1:2){
  par(mfrow=c(3,2))
  for(i in 1:3){
    if(j == 1){
      hist(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==0)]==i)],main = EarlyCompTypeTitles[i],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
      hist(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==0)]==i)],main = EarlyCompTypeTitles[i+3],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
    }
    if(j == 2){
      hist(SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==0)]==i+3)],main = LateCompTypeTitles[i],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
      hist(OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==0)]==i+3)],main = LateCompTypeTitles[i+3],xlab = "",
           ylim = c(0,15), xlim = c(0,1),breaks = 10
      )
    }
  }
}

## Combine early and late by composition type

histMegaYesNovice <- c(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==1)]==1)],
                       SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==1)]==4)])

histMegaYesExpert <- c(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]==1)],
                       OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==1)]==4)])

histBorderlineYesNovice <- c(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==1)]==2)],
                             SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==1)]==5)])

histBorderlineYesExpert <- c(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]==2)],
                             OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==1)]==5)])

histMegaNotYesNovice <- c(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==1)]==3)],
                          SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==1)]==6)])

histMegaNotYesExpert <- c(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]==3)],
                          OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==1)]==6)])

histMegaNoNovice <- c(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==0)]==1)],
                      SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==0)]==4)])

histMegaNoExpert <- c(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==0)]==1)],
                      OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==0)]==4)])

histBorderlineNoNovice <- c(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==0)]==2)],
                            SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==0)]==5)])

histBorderlineNoExpert <- c(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==0)]==2)],
                            OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==0)]==5)])

histMegaNotNoNovice <- c(SONAData$EarlyResponse[which(SONAData$CompType[which(SONAData$SpeakerEarly==0)]==3)],
                         SONAData$LateResponse[which(SONAData$CompType[which(SONAData$SpeakerLate==0)]==6)])

histMegaNotNoExpert <- c(OnlineData$EarlyResponse[which(OnlineData$CompType[which(OnlineData$SpeakerEarly==0)]==3)],
                         OnlineData$LateResponse[which(OnlineData$CompType[which(OnlineData$SpeakerLate==0)]==6)])

histBorderlineYes <- list(histMegaYesNovice,histMegaYesExpert,histBorderlineYesNovice,histBorderlineYesExpert,
                          histMegaNotYesNovice,histMegaNotYesExpert)

histBorderlineNo <- list(histMegaNoNovice,histMegaNoExpert,histBorderlineNoNovice,histBorderlineNoExpert,
                         histMegaNotNoNovice,histMegaNotNoExpert)

YesCompTypeTitles= c('Mega Yes Novice','Mega Yes Expert','Borderline Yes Novice','Borderline Yes Expert',
                     'Mega Not Yes Novice','Mega Not Yes Expert')

NoCompTypeTitles= c('Mega No Novice','Mega No Expert','Borderline No Novice','Borderline No Expert',
                    'Mega Not No Novice','Mega Not No Expert')

for(j in 1:2){
  par(mfrow=c(3,2))
  for(i in 1:6){
    if(j == 1){
      hist(histBorderlineYes[[i]],main = YesCompTypeTitles[i],xlab = "",
           ylim = c(0,30), xlim = c(0,1),breaks = 10
      )
      abline(v=mean(histBorderlineYes[[i]]))
    }
    if(j == 2){
      hist(histBorderlineNo[[i]],main = NoCompTypeTitles[i],xlab = "",
           ylim = c(0,30), xlim = c(0,1),breaks = 10
      )
      abline(v=mean(histBorderlineNo[[i]]))
    }
  }
}

histMegaNovicePrior <- c(SONAData$EarlyPrior[which(SONAData$CompType==1)],
                         SONAData$LatePrior[which(SONAData$CompType==4)])

histMegaExpertPrior <- c(OnlineData$EarlyPrior[which(OnlineData$CompType==1)],
                         OnlineData$LatePrior[which(OnlineData$CompType==4)])

histBorderlineNovicePrior <- c(SONAData$EarlyPrior[which(SONAData$CompType==2)],
                               SONAData$LatePrior[which(SONAData$CompType==5)])

histBorderlineExpertPrior <- c(OnlineData$EarlyPrior[which(OnlineData$CompType==2)],
                               OnlineData$LatePrior[which(OnlineData$CompType==5)])

histMegaNotNovicePrior <- c(SONAData$EarlyPrior[which(SONAData$CompType==3)],
                            SONAData$LatePrior[which(SONAData$CompType==6)])

histMegaNotExpertPrior <- c(OnlineData$EarlyPrior[which(OnlineData$CompType==3)],
                            OnlineData$LatePrior[which(OnlineData$CompType==6)])

histBorderlinePrior <- list(histMegaNovicePrior,histMegaExpertPrior,histBorderlineNovicePrior,histBorderlineExpertPrior,
                            histMegaNotNovicePrior,histMegaNotExpertPrior)

PriorCompTypeTitles= c('Mega Novice Prior','Mega Expert Prior','Borderline Novice Prior','Borderline Expert Prior',
                       'Mega Not Novice Prior','Mega Not Expert Prior')

## Priors, combine early and late by composition type
par(mfrow=c(3,2))
for(i in 1:6){
  hist(histBorderlinePrior[[i]],main = PriorCompTypeTitles[i],xlab = "",
       ylim = c(0,60), xlim = c(0,1),breaks = 10
  )
  abline(v=mean(histBorderlinePrior[[i]]))
}

ConditionSize = length(which(OnlineData$SpeakerEarly==1))+length(which(SONAData$SpeakerEarly==1))
NoviceConditionSize = length(which(SONAData$SpeakerEarly==1))
ExpertConditionSize = length(which(OnlineData$SpeakerEarly==1))
ExpertStart = NoviceConditionSize + 1

NoviceIDs <- rep(1:nNovices,each=6)
IDStart = 1+nNovices
IDEnd = nNovices+nExperts
ExpertIDs <- rep(IDStart:IDEnd,each=6)

## Look for participants doing weird stuff
SONAData$FromSpeakerEarly <- abs(SONAData$SpeakerEarly - SONAData$EarlyResponse)
SONAData$FromSpeakerLate <- abs(SONAData$SpeakerLate - SONAData$LateResponse)

SerialOffendersNoviceEarly <- SONAData$Participant[which(SONAData$FromSpeakerEarly>.8)]
SerialOffendersNoviceLate <- SONAData$Participant[which(SONAData$FromSpeakerLate>.8)]

# SerialOffendersNoviceEarly<- rbind(SerialOffendersNoviceEarly,which(SONAData$FromSpeakerEarly>.8)%%12)
# SerialOffendersNoviceLate<-rbind(SerialOffendersNoviceLate,which(SONAData$FromSpeakerLate>.8)%%12)


# Correlation Matrix of Compositions --------------------------------------


RawExpertCorMat_EarlyYes <- matrix(NaN,nrow=nExperts,ncol=6)
RawExpertCorMat_LateYes <- matrix(NaN,nrow=nExperts,ncol=6)
RawExpertCorMat_EarlyNo <- matrix(NaN,nrow=nExperts,ncol=6)
RawExpertCorMat_LateNo <- matrix(NaN,nrow=nExperts,ncol=6)
for(i in 1:nExperts){
  for(j in 1:nItems){
    ticker = ((i-1)*nItems)+j
    if(OnlineData$SpeakerEarly[ticker]==1){
      RawExpertCorMat_EarlyYes[i,OnlineData$CompType[ticker]] <- OnlineData$ProbitEarly[ticker]
    }
    if(OnlineData$SpeakerEarly[ticker]==0){
      RawExpertCorMat_EarlyNo[i,OnlineData$CompType[ticker]] <- OnlineData$ProbitEarly[ticker]
    }
    if(OnlineData$SpeakerLate[ticker]==1){
      RawExpertCorMat_LateYes[i,OnlineData$CompType[ticker]] <- OnlineData$ProbitLate[ticker]
    }
    if(OnlineData$SpeakerLate[ticker]==0){
      RawExpertCorMat_LateNo[i,OnlineData$CompType[ticker]] <- OnlineData$ProbitLate[ticker]
    }
  }
}

RawExpertCorMat_EarlyYes <- data.frame(RawExpertCorMat_EarlyYes)
RawExpertCorMat_LateYes <- data.frame(RawExpertCorMat_LateYes)
names(RawExpertCorMat_EarlyYes) <- c("mE","bE","nE","mL","bL","nL")
names(RawExpertCorMat_LateYes) <- c("mE","bE","nE","mL","bL","nL")

RawExpertCorMat_EarlyNo <- data.frame(RawExpertCorMat_EarlyNo)
RawExpertCorMat_LateNo <- data.frame(RawExpertCorMat_LateNo)
names(RawExpertCorMat_EarlyNo) <- c("mE","bE","nE","mL","bL","nL")
names(RawExpertCorMat_LateNo) <- c("mE","bE","nE","mL","bL","nL")

RawNoviceCorMat_EarlyYes <- matrix(NaN,nrow=nNovices,ncol=6)
RawNoviceCorMat_LateYes <- matrix(NaN,nrow=nNovices,ncol=6)
RawNoviceCorMat_EarlyNo <- matrix(NaN,nrow=nNovices,ncol=6)
RawNoviceCorMat_LateNo <- matrix(NaN,nrow=nNovices,ncol=6)
for(i in 1:nNovices){
  for(j in 1:nItems){
    ticker = ((i-1)*nItems)+j
    if(SONAData$SpeakerEarly[ticker]==1){
      RawNoviceCorMat_EarlyYes[i,SONAData$CompType[ticker]] <- SONAData$ProbitEarly[ticker]
    }
    if(SONAData$SpeakerEarly[ticker]==0){
      RawNoviceCorMat_EarlyNo[i,SONAData$CompType[ticker]] <- SONAData$ProbitEarly[ticker]
    }
    if(SONAData$SpeakerLate[ticker]==1){
      RawNoviceCorMat_LateYes[i,SONAData$CompType[ticker]] <- SONAData$ProbitLate[ticker]
    }
    if(SONAData$SpeakerLate[ticker]==0){
      RawNoviceCorMat_LateNo[i,SONAData$CompType[ticker]] <- SONAData$ProbitLate[ticker]
    }
  }
}

RawNoviceCorMat_EarlyYes <- data.frame(RawNoviceCorMat_EarlyYes)
RawNoviceCorMat_LateYes <- data.frame(RawNoviceCorMat_LateYes)
names(RawNoviceCorMat_EarlyYes) <- c("mE","bE","nE","mL","bL","nL")
names(RawNoviceCorMat_LateYes) <- c("mE","bE","nE","mL","bL","nL")

RawNoviceCorMat_EarlyNo <- data.frame(RawNoviceCorMat_EarlyNo)
RawNoviceCorMat_LateNo <- data.frame(RawNoviceCorMat_LateNo)
names(RawNoviceCorMat_EarlyNo) <- c("mE","bE","nE","mL","bL","nL")
names(RawNoviceCorMat_LateNo) <- c("mE","bE","nE","mL","bL","nL")

par(mfrow=c(2,2))
NoviceCorMat_EarlyYes <- cor(RawNoviceCorMat_EarlyYes)
NoviceCorMat_LateYes <- cor(RawNoviceCorMat_LateYes)
ExpertCorMat_EarlyYes <- cor(RawExpertCorMat_EarlyYes)
ExpertCorMat_LateYes <- cor(RawExpertCorMat_LateYes)
corrplot(NoviceCorMat_EarlyYes, method="number",type = 'upper',title = 'Novice Early Yes')
corrplot(ExpertCorMat_EarlyYes, method="number",type = 'upper', title = 'Expert Early Yes')
corrplot(NoviceCorMat_LateYes, method="number",type = 'upper',title = 'Novice Late Yes')
corrplot(ExpertCorMat_LateYes, method="number",type = 'upper', title = 'Expert Late Yes')

par(mfrow=c(2,2))
NoviceCorMat_EarlyNo <- cor(RawNoviceCorMat_EarlyNo)
NoviceCorMat_LateNo <- cor(RawNoviceCorMat_LateNo)
ExpertCorMat_EarlyNo <- cor(RawExpertCorMat_EarlyNo)
ExpertCorMat_LateNo <- cor(RawExpertCorMat_LateNo)
corrplot(NoviceCorMat_EarlyNo, method="number",type = 'upper',title = 'Novice Early No')
corrplot(ExpertCorMat_EarlyNo, method="number",type = 'upper', title = 'Expert Early No')
corrplot(NoviceCorMat_LateNo, method="number",type = 'upper',title = 'Novice Late No')
corrplot(ExpertCorMat_LateNo, method="number",type = 'upper', title = 'Expert Late No')

# Prep data for raw JAGS inputs -----------------------------------------------
nParticipants = nExperts+nNovices
nTrials = 6
nExpResponses = nExperts*nTrials
nNovResponses = nNovices*nTrials
nBeta = 7

XyesEarly <- matrix(data=rep(0,nParticipants*6*nBeta),nrow = nParticipants*6,ncol = nBeta)
YyesEarly <- rep(NaN,nParticipants*6)
zyesEarly <- matrix(data=rep(0,nParticipants*6*2),nrow = nParticipants*6,ncol = 2)

XnoEarly <- matrix(data=rep(0,nParticipants*6*nBeta),nrow = nParticipants*6,ncol = nBeta)
YnoEarly <- rep(NaN,nParticipants*6)
znoEarly <- matrix(data=rep(0,nParticipants*6*2),nrow = nParticipants*6,ncol = 2)

XyesLate <- matrix(data=rep(0,nParticipants*6*nBeta),nrow = nParticipants*6,ncol = nBeta)
YyesLate <- rep(NaN,nParticipants*6)
zyesLate <- matrix(data=rep(0,nParticipants*6*2),nrow = nParticipants*6,ncol = 2)

XnoLate <- matrix(data=rep(0,nParticipants*6*nBeta),nrow = nParticipants*6,ncol = nBeta)
YnoLate <- rep(NaN,nParticipants*6)
znoLate <- matrix(data=rep(0,nParticipants*6*2),nrow = nParticipants*6,ncol = 2)

XyesEarly[(nExpResponses+1):(nParticipants*6),1] <- 1
zyesEarly[(nExpResponses+1):(nParticipants*6),2] <- 1
zyesEarly[1:nExpResponses,1] <- 1

XnoEarly[(nExpResponses+1):(nParticipants*6),1] <- 1
znoEarly[(nExpResponses+1):(nParticipants*6),2] <- 1
znoEarly[1:nExpResponses,1] <- 1

XyesLate[(nExpResponses+1):(nParticipants*6),1] <- 1
zyesLate[(nExpResponses+1):(nParticipants*6),2] <- 1
zyesLate[1:nExpResponses,1] <- 1

XnoLate[(nExpResponses+1):(nParticipants*6),1] <- 1
znoLate[(nExpResponses+1):(nParticipants*6),2] <- 1
znoLate[1:nExpResponses,1] <- 1

tickeryesEarly = 1
tickernoEarly = 1
tickeryesLate = 1
tickernoLate = 1

for(i in 1:nrow(OnlineData)){
  
  
  if(OnlineData$SpeakerEarly[i]==1){
    
    YyesEarly[tickeryesEarly] <- OnlineData$EarlyResponse[i]
    
    XyesEarly[tickeryesEarly,(OnlineData$CompType[i]+1)] <- 1
    
    tickeryesEarly = tickeryesEarly+1
  }
  
  if(OnlineData$SpeakerEarly[i]==0){
    
    YnoEarly[tickernoEarly] <- OnlineData$EarlyResponse[i]
    
    XnoEarly[tickernoEarly,(OnlineData$CompType[i]+1)] <- 1
    
    tickernoEarly = tickernoEarly+1
  }
  
  if(OnlineData$SpeakerLate[i]==1){
    
    YyesLate[tickeryesLate] <- OnlineData$LateResponse[i]
    
    XyesLate[tickeryesLate,(OnlineData$CompType[i]+1)] <- 1
    
    tickeryesLate = tickeryesLate+1
  }
  
  if(OnlineData$SpeakerLate[i]==0){
    
    YnoLate[tickernoLate] <- OnlineData$LateResponse[i]
    
    XnoLate[tickernoLate,(OnlineData$CompType[i]+1)] <- 1
    
    tickernoLate = tickernoLate+1
  }
}

for(i in 1:nrow(SONAData)){
  
  
  if(SONAData$SpeakerEarly[i]==1){
    
    YyesEarly[tickeryesEarly] <- SONAData$EarlyResponse[i]
    
    #XyesEarly[tickeryesEarly,(SONAData$CompType[i]+1)] <- 1
    
    tickeryesEarly = tickeryesEarly+1
  }
  
  if(SONAData$SpeakerEarly[i]==0){
    
    YnoEarly[tickernoEarly] <- SONAData$EarlyResponse[i]
    
    #XnoEarly[tickernoEarly,(SONAData$CompType[i]+1)] <- 1
    
    tickernoEarly = tickernoEarly+1
  }
  
  if(SONAData$SpeakerLate[i]==1){
    
    YyesLate[tickeryesLate] <- SONAData$LateResponse[i]
    
    #XyesLate[tickeryesLate,(SONAData$CompType[i]+1)] <- 1
    
    tickeryesLate = tickeryesLate+1
  }
  
  if(SONAData$SpeakerLate[i]==0){
    
    YnoLate[tickernoLate] <- SONAData$LateResponse[i]
    
    #XnoLate[tickernoLate,(SONAData$CompType[i]+1)] <- 1
    
    tickernoLate = tickernoLate+1
  }
}


# MVN Precision Matrix (for raw JAGS betas) ----------------------------------------
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Run raw linear regression in JAGS -------------------------------------------
nSamples = 4000

rawBetayesEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
rawDIC <- rep(NaN,4)

X <- XyesEarly
Y<- YyesEarly
z <- zyesEarly

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

rawBetayesEarlysamples <- samples$BUGSoutput$sims.list$BETA

rawsigma_yesEarlysamples <- samples$BUGSoutput$sims.list$sigma

rawtau_yesEarlysamples <- samples$BUGSoutput$sims.list$tau

rawsummaryyesEarly <- samples$BUGSoutput$summary

rawDIC[1] <- samples$BUGSoutput$DIC

##No Early Regression##
rawBetanoEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoEarly
Y<- YnoEarly
z <- znoEarly

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

rawBetanoEarlysamples <- samples$BUGSoutput$sims.list$BETA

rawsigma_noEarlysamples <- samples$BUGSoutput$sims.list$sigma

rawtau_noEarlysamples <- samples$BUGSoutput$sims.list$tau

rawsummarynoEarly <- samples$BUGSoutput$summary

rawDIC[2] <- samples$BUGSoutput$DIC

##Yes Late Regression##

rawBetayesLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XyesLate
Y<- YyesLate
z <- zyesLate

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

rawBetayesLatesamples <- samples$BUGSoutput$sims.list$BETA

rawsigma_yesLatesamples <- samples$BUGSoutput$sims.list$sigma

rawtau_yesLatesamples <- samples$BUGSoutput$sims.list$tau

rawsummaryyesLate <- samples$BUGSoutput$summary

rawDIC[3] <- samples$BUGSoutput$DIC

##No Late Regression##
rawBetanoLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoLate
Y<- YnoLate
z <- znoLate

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

rawBetanoLatesamples <- samples$BUGSoutput$sims.list$BETA

rawsigma_noLatesamples <- samples$BUGSoutput$sims.list$sigma

rawtau_noLatesamples <- samples$BUGSoutput$sims.list$tau

rawsummarynoLate <- samples$BUGSoutput$summary

rawDIC[4] <- samples$BUGSoutput$DIC

# Run split sigma linear regression in JAGS -------------------------------------------
nSamples = 4000

splitsig_BetayesEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
splitsig_DIC <- rep(NaN,4)

X <- XyesEarly
Y<- YyesEarly
z <- zyesEarly

data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta","z") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_nov","sigma_exp","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_jags_splitsig.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

splitsig_BetayesEarlysamples <- samples$BUGSoutput$sims.list$BETA

sigmanov_yesEarlysamples <- samples$BUGSoutput$sims.list$sigma_nov

sigmaexp_yesEarlysamples <- samples$BUGSoutput$sims.list$sigma_exp

tau_yesEarlysamples <- samples$BUGSoutput$sims.list$tau

splitsig_summaryyesEarly <- samples$BUGSoutput$summary

splitsig_DIC[1] <- samples$BUGSoutput$DIC

##No Early Regression##
splitsig_BetanoEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoEarly
Y<- YnoEarly
z <- znoEarly

data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta","z") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_nov","sigma_exp","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_jags_splitsig.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

splitsig_BetanoEarlysamples <- samples$BUGSoutput$sims.list$BETA

sigmanov_noEarlysamples <- samples$BUGSoutput$sims.list$sigma_nov

sigmaexp_noEarlysamples <- samples$BUGSoutput$sims.list$sigma_exp

tau_noEarlysamples <- samples$BUGSoutput$sims.list$tau

splitsig_summarynoEarly <- samples$BUGSoutput$summary

splitsig_DIC[2] <- samples$BUGSoutput$DIC

##Yes Late Regression##

splitsig_BetayesLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XyesLate
Y<- YyesLate
z <- zyesLate

data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta","z") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_nov","sigma_exp","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_jags_splitsig.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

splitsig_BetayesLatesamples <- samples$BUGSoutput$sims.list$BETA

sigmanov_yesLatesamples <- samples$BUGSoutput$sims.list$sigma_nov

sigmaexp_yesLatesamples <- samples$BUGSoutput$sims.list$sigma_exp

tau_yesLatesamples <- samples$BUGSoutput$sims.list$tau

splitsig_summaryyesLate <- samples$BUGSoutput$summary

splitsig_DIC[3] <- samples$BUGSoutput$DIC

##No Late Regression##
splitsig_BetanoLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoLate
Y<- YnoLate
z <- znoLate

data <- list("nParticipants","nTrials","X", "Y","bPrec","nBeta","z") # to be passed on to JAGS
# parameters to be monitored:	
parameters <- c("BETA","tau","sigma_nov","sigma_exp","alpha","Ypred")

#Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
myinits <- list(list("BETA" = runif(7,-.5,.5)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, #inits=myinits, 
                parameters.to.save=parameters,
                model.file="LoLComp_LinearReg_jags_splitsig.txt", n.chains=3,n.burnin = 1000, n.iter=3000, DIC=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

splitsig_BetanoLatesamples <- samples$BUGSoutput$sims.list$BETA

sigmanov_noLatesamples <- samples$BUGSoutput$sims.list$sigma_nov

sigmaexp_noLatesamples <- samples$BUGSoutput$sims.list$sigma_exp

tau_noLatesamples <- samples$BUGSoutput$sims.list$tau

splitsig_summarynoLate <- samples$BUGSoutput$summary

splitsig_DIC[4] <- samples$BUGSoutput$DIC


# Analyze sigmas and taus from standard model -----------------------------
RhoTest_yesEarly <- rep(NaN,length(rawsigma_yesEarlysamples))
RhoTest_noEarly <- matrix(NaN,length(rawsigma_noEarlysamples))
RhoTest_yesLate <- matrix(NaN,length(rawsigma_yesLatesamples))
RhoTest_noLate <- matrix(NaN,length(rawsigma_noLatesamples))

for(i in 1:length(rawsigma_yesEarlysamples)){
  RhoTest_yesEarly[i] <- rawtau_yesEarlysamples[i]^2/
    (rawtau_yesEarlysamples[i]^2 + rawsigma_yesEarlysamples[i]^2)
  
  RhoTest_yesLate[i] <- rawtau_yesLatesamples[i]^2/
    (rawtau_yesLatesamples[i]^2 + rawsigma_yesLatesamples[i]^2)
}

ind975 = length(rawsigma_yesEarlysamples)*.975
ind075 = length(rawsigma_yesEarlysamples)*.075

sortRhoTest_yesEarly <- sort(RhoTest_yesEarly)
sortRhoTest_yesLate <- sort(RhoTest_yesLate)

mean(RhoTest_yesEarly)
sortRhoTest_yesEarly[ind975]
sortRhoTest_yesEarly[ind075]

mean(RhoTest_yesLate)
sortRhoTest_yesLate[ind975]
sortRhoTest_yesLate[ind075]

# Compare expert and novice sigmas ----------------------------------------

SigmaTest_yesEarly <- matrix(NaN,nrow = length(sigmanov_yesEarlysamples),ncol = 2)
SigmaTest_noEarly <- matrix(NaN,nrow = length(sigmanov_noEarlysamples),ncol = 2)
SigmaTest_yesLate <- matrix(NaN,nrow = length(sigmanov_yesLatesamples),ncol = 2)
SigmaTest_noLate <- matrix(NaN,nrow = length(sigmanov_noLatesamples),ncol = 2)

SigmaDiff_yesEarly <- rep(NaN,length(sigmanov_yesEarlysamples))
SigmaDiff_noEarly <- rep(NaN,length(sigmanov_yesEarlysamples))
SigmaDiff_yesLate <- rep(NaN,length(sigmanov_yesEarlysamples))
SigmaDiff_noLate <- rep(NaN,length(sigmanov_yesEarlysamples))

for(i in 1:length(sigmanov_noLatesamples)){
  SigmaTest_yesEarly[i,1] <- tau_yesEarlysamples[i]^2/
    (tau_yesEarlysamples[i]^2 + sigmaexp_yesEarlysamples[i]^2)
  SigmaTest_yesEarly[i,2] <- tau_yesEarlysamples[i]^2/
    (tau_yesEarlysamples[i]^2 + sigmanov_yesEarlysamples[i]^2)
  
  SigmaDiff_yesEarly[i] <- sigmaexp_yesEarlysamples[i]-sigmanov_yesEarlysamples[i]
  SigmaDiff_yesLate[i] <- sigmaexp_yesLatesamples[i]-sigmanov_yesLatesamples[i]
  
  SigmaTest_yesLate[i,1] <- tau_yesLatesamples[i]^2/
    (tau_yesLatesamples[i]^2 + sigmaexp_yesLatesamples[i]^2)
  SigmaTest_yesLate[i,2] <- tau_yesEarlysamples[i]^2/
    (tau_yesLatesamples[i]^2 + sigmanov_yesLatesamples[i]^2)
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

sortSigmaTest_yesEarly <- matrix(NaN,nrow = length(sigmanov_yesEarlysamples),ncol = 2)
sortSigmaTest_noEarly <- matrix(NaN,nrow = length(sigmanov_noEarlysamples),ncol = 2)
sortSigmaTest_yesLate <- matrix(NaN,nrow = length(sigmanov_yesLatesamples),ncol = 2)
sortSigmaTest_noLate <- matrix(NaN,nrow = length(sigmanov_noLatesamples),ncol = 2)

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
posterior_density_sigDiff_noEarly <- dlogspline(0,logspline(SigmaDiff_noEarly))
posterior_density_sigDiff_noLate <- dlogspline(0,logspline(SigmaDiff_noLate))

SavDick <- data.frame(matrix(NaN,nrow=6,ncol=8))
names(SavDick) <- c("Yes Early","12B Yes Early","No Early","12B No Early","Yes Late","12B Yes Late",
                    "No Late","12B No Late")

SavDick_sig <- c(posterior_density_sigDiff_yesEarly/prior_density_SigDiff,
                 posterior_density_sigDiff_yesLate/prior_density_SigDiff)

# Plot beta inferences from raw -------------------------------------------
rawPlotValues_EarlyYes <- matrix(NaN, nrow = 7,ncol = 3)
rawPlotValues_EarlyNo <- matrix(NaN, nrow = 7,ncol = 3)
rawPlotValues_LateYes <- matrix(NaN, nrow = 7,ncol = 3)
rawPlotValues_LateNo <- matrix(NaN, nrow = 7,ncol = 3)

loconf = nrow(rawBetayesEarlysamples) * .025
hiconf = nrow(rawBetayesEarlysamples) * .975


tempsort <- sort(rawBetayesEarlysamples[,1])
rawPlotValues_EarlyYes[1,1] <- mean(tempsort[1500:1001])
rawPlotValues_EarlyYes[1,2] <- tempsort[loconf]
rawPlotValues_EarlyYes[1,3] <- tempsort[hiconf]

tempsort <- sort(rawBetanoEarlysamples[,1])
rawPlotValues_EarlyNo[1,1] <- mean(tempsort[1500:1001])
rawPlotValues_EarlyNo[1,2] <- tempsort[loconf]
rawPlotValues_EarlyNo[1,3] <- tempsort[hiconf]

tempsort <- sort(rawBetayesLatesamples[,1])
rawPlotValues_LateYes[1,1] <- mean(tempsort[1500:1001])
rawPlotValues_LateYes[1,2] <- tempsort[loconf]
rawPlotValues_LateYes[1,3] <- tempsort[hiconf]

tempsort <- sort(rawBetanoLatesamples[,1])
rawPlotValues_LateNo[1,1] <- mean(tempsort[1500:1001])
rawPlotValues_LateNo[1,2] <- tempsort[loconf]
rawPlotValues_LateNo[1,3] <- tempsort[hiconf]

for(i in 2:7){
  
  tempsort <- sort(rawBetayesEarlysamples[,i])
  rawPlotValues_EarlyYes[i,1] <- mean(tempsort[1500:1001])
  rawPlotValues_EarlyYes[i,2] <- tempsort[loconf]
  rawPlotValues_EarlyYes[i,3] <- tempsort[hiconf]
  
  tempsort <- sort(rawBetanoEarlysamples[,i])
  rawPlotValues_EarlyNo[i,1] <- mean(tempsort[1500:1001])
  rawPlotValues_EarlyNo[i,2] <- tempsort[loconf]
  rawPlotValues_EarlyNo[i,3] <- tempsort[hiconf]
  
  tempsort <- sort(rawBetayesLatesamples[,i])
  rawPlotValues_LateYes[i,1] <- mean(tempsort[1500:1001])
  rawPlotValues_LateYes[i,2] <- tempsort[loconf]
  rawPlotValues_LateYes[i,3] <- tempsort[hiconf]
  
  tempsort <- sort(rawBetanoLatesamples[,i])
  rawPlotValues_LateNo[i,1] <- mean(tempsort[1500:1001])
  rawPlotValues_LateNo[i,2] <- tempsort[loconf]
  rawPlotValues_LateNo[i,3] <- tempsort[hiconf]
}

plot(1:6,rawPlotValues_EarlyYes[2:7,1],ylim = c(.15,.9),main = "Yes Early",xlab = 'beta',ylab = 'Interpretation')
abline(h = rawPlotValues_EarlyYes[1,1],col = '#DC143C')
abline(h= rawPlotValues_EarlyYes[1,2],col = '#DC143C',lty = 2)
abline(h= rawPlotValues_EarlyYes[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = rawPlotValues_EarlyYes[i+1,2],
         x1 = i, y1 = rawPlotValues_EarlyYes[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}
legend(3.5,.45,legend = c('Novice Interpretation','Expert Interpretation'),pch = c(NA,1),lty = c(1,NA),
       col = c('#DC143C','black'),cex = .8)

plot(1:6,rawPlotValues_LateYes[2:7,1],ylim = c(.15,.9),main = "Yes Late",xlab = 'beta',ylab = 'Interpretation')
abline(h = rawPlotValues_LateYes[1,1],col = '#DC143C')
abline(h= rawPlotValues_LateYes[1,2],col = '#DC143C',lty = 2)
abline(h= rawPlotValues_LateYes[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = rawPlotValues_LateYes[i+1,2],
         x1 = i, y1 = rawPlotValues_LateYes[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}


plot(1:6,rawPlotValues_EarlyNo[2:7,1],ylim = c(.15,.9),main = "No Early",xlab = 'beta',ylab = 'Interpretation')
abline(h = rawPlotValues_EarlyNo[1,1],col = '#DC143C')
abline(h= rawPlotValues_EarlyNo[1,2],col = '#DC143C',lty = 2)
abline(h= rawPlotValues_EarlyNo[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = rawPlotValues_EarlyNo[i+1,2],
         x1 = i, y1 = rawPlotValues_EarlyNo[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}

plot(1:6,rawPlotValues_LateNo[2:7,1],ylim = c(.15,.9),main = "No Late",xlab = 'beta',ylab = 'Interpretation')
abline(h = rawPlotValues_LateNo[1,1],col = '#DC143C')
abline(h= rawPlotValues_LateNo[1,2],col = '#DC143C',lty = 2)
abline(h= rawPlotValues_LateNo[1,3],col = '#DC143C',lty = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = rawPlotValues_LateNo[i+1,2],
         x1 = i, y1 = rawPlotValues_LateNo[i+1,3],code = 3,
         col = 'black',angle = 90, length = .1,lty = 3)
}

# Test for expected order in "Yes" conditions -----------------------------

OrderTestEarly <- matrix(NaN,nrow=nrow(rawBetayesEarlysamples),ncol = 7)

OrderTestEarly[,1:6] <- rawBetayesEarlysamples[,2:7]

##Early hypothesis: 1 > 2 > 3, 6 > 5 > 4
for(i in 1:nrow(rawBetayesEarlysamples)){
  firstcheck <- OrderTestEarly[i,1] > OrderTestEarly[i,2] && 
    OrderTestEarly[i,2] > OrderTestEarly[i,3]
  secondcheck <- OrderTestEarly[i,6] > OrderTestEarly[i,5] && 
    OrderTestEarly[i,5] > OrderTestEarly[i,4]
  OrderTestEarly[i,7] <- firstcheck && secondcheck
}

OrderTestLate <- matrix(NaN,nrow=nrow(rawBetayesLatesamples),ncol = 7)

OrderTestLate[,1:6] <- rawBetayesLatesamples[,2:7]

##Late hypothesis: 3 > 2 > 1, 4 > 5 > 6
for(i in 1:nrow(rawBetayesLatesamples)){
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

# Prep data for overfull JAGS inputs -----------------------------------------------
nParticipants = nExperts+nNovices
nBeta = 12
nTrials = 6
nExpResponses = nExperts*nTrials
nNovResponses = nNovices*nTrials

XyesEarly_overfull <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YyesEarly <- rep(NaN,nParticipants*nTrials)

XnoEarly_overfull <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YnoEarly <- rep(NaN,nParticipants*nTrials)

XyesLate_overfull <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YyesLate <- rep(NaN,nParticipants*nTrials)

XnoLate_overfull <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YnoLate <- rep(NaN,nParticipants*nTrials)

tickeryesEarly = 1
tickernoEarly = 1
tickeryesLate = 1
tickernoLate = 1

for(i in 1:nrow(OnlineData)){
  
  
  if(OnlineData$SpeakerEarly[i]==1){
    
    YyesEarly[tickeryesEarly] <- OnlineData$EarlyResponse[i]
    
    XyesEarly_overfull[tickeryesEarly,(OnlineData$CompType[i])] <- 1
    
    tickeryesEarly = tickeryesEarly+1
  }
  
  if(OnlineData$SpeakerEarly[i]==0){
    
    YnoEarly[tickernoEarly] <- OnlineData$EarlyResponse[i]
    
    XnoEarly_overfull[tickernoEarly,(OnlineData$CompType[i])] <- 1
    
    tickernoEarly = tickernoEarly+1
  }
  
  if(OnlineData$SpeakerLate[i]==1){
    
    YyesLate[tickeryesLate] <- OnlineData$LateResponse[i]
    
    XyesLate_overfull[tickeryesLate,(OnlineData$CompType[i])] <- 1
    
    tickeryesLate = tickeryesLate+1
  }
  
  if(OnlineData$SpeakerLate[i]==0){
    
    YnoLate[tickernoLate] <- OnlineData$LateResponse[i]
    
    XnoLate_overfull[tickernoLate,(OnlineData$CompType[i])] <- 1
    
    tickernoLate = tickernoLate+1
  }
}

for(i in 1:nrow(SONAData)){
  
  
  if(SONAData$SpeakerEarly[i]==1){
    
    YyesEarly[tickeryesEarly] <- SONAData$EarlyResponse[i]
    
    XyesEarly_overfull[tickeryesEarly,(SONAData$CompType[i]+6)] <- 1
    
    tickeryesEarly = tickeryesEarly+1
  }
  
  if(SONAData$SpeakerEarly[i]==0){
    
    YnoEarly[tickernoEarly] <- SONAData$EarlyResponse[i]
    
    XnoEarly_overfull[tickernoEarly,(SONAData$CompType[i]+6)] <- 1
    
    tickernoEarly = tickernoEarly+1
  }
  
  if(SONAData$SpeakerLate[i]==1){
    
    YyesLate[tickeryesLate] <- SONAData$LateResponse[i]
    
    XyesLate_overfull[tickeryesLate,(SONAData$CompType[i]+6)] <- 1
    
    tickeryesLate = tickeryesLate+1
  }
  
  if(SONAData$SpeakerLate[i]==0){
    
    YnoLate[tickernoLate] <- SONAData$LateResponse[i]
    
    XnoLate_overfull[tickernoLate,(SONAData$CompType[i]+6)] <- 1
    
    tickernoLate = tickernoLate+1
  }
}


# MVN Precision Matrix (for overfull JAGS betas) ----------------------------------------
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Run overfull linear regression in JAGS -------------------------------------------
nSamples = 4000

overfullBetayesEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
overfullDIC <- rep(NaN,4)

X <- XyesEarly_overfull
Y<- YyesEarly

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

overfullBetayesEarlysamples <- samples$BUGSoutput$sims.list$BETA

overfullsummaryyesEarly <- samples$BUGSoutput$summary

overfullDIC[1] <- samples$BUGSoutput$DIC

##No Early Regression##
overfullBetanoEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoEarly_overfull
Y<- YnoEarly

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

overfullBetanoEarlysamples <- samples$BUGSoutput$sims.list$BETA

overfullsummarynoEarly <- samples$BUGSoutput$summary

overfullDIC[2] <- samples$BUGSoutput$DIC

##Yes Late Regression##

overfullBetayesLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XyesLate_overfull
Y<- YyesLate

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

overfullBetayesLatesamples <- samples$BUGSoutput$sims.list$BETA

overfullsummaryyesLate <- samples$BUGSoutput$summary

overfullDIC[3] <- samples$BUGSoutput$DIC

##No Late Regression##
overfullBetanoLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoLate_overfull
Y<- YnoLate

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

overfullBetanoLatesamples <- samples$BUGSoutput$sims.list$BETA

overfullsummarynoLate <- samples$BUGSoutput$summary

overfullDIC[4] <- samples$BUGSoutput$DIC



# Compare DIC's (between main and overfull) -----------------------------------------------------------

par(mfrow=c(1,1))

plot(1:4,rawDIC,ylim = c(-700,-300),xlab = 'Condition',ylab = 'DIC', main = "Main Model and Overfull DIC Comparison")

points(1:4,overfullDIC,col = '#DC143C',pch=3)

legend(3, -200, legend=c("Main Model", "Overfull"),
       col=c("black", "red"), pch=c(1,3), cex=0.8)

text_rawDIC <- matrix(NaN,nrow = 4, ncol =2)
text_rawDIC[,1] <- 1:4
text_rawDIC[,2] <- round(rawDIC,0)
text_rawDIC <- data.frame(text_rawDIC)
names(text_rawDIC) <- c("Condition","DIC")

text(DIC~Condition, labels=DIC,data=text_rawDIC, cex=0.9, font=1,pos=1)

text_overfullDIC <- matrix(NaN,nrow = 4, ncol =2)
text_overfullDIC[,1] <- 1:4
text_overfullDIC[,2] <- round(overfullDIC,0)
text_overfullDIC <- data.frame(text_overfullDIC)
names(text_overfullDIC) <- c("Condition","DIC")

text(DIC~Condition, labels=DIC,data=text_overfullDIC, cex=0.9, font=1,pos=3,col='#DC143C')


# Compare alphas (between main and overfull) ------------------------------

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

par(mfrow=c(2,2))
b <- min(c(rawsummaryyesEarly[8:56,5],rawsummaryyesEarly[57:105,5])) - .001
e <- max(c(rawsummaryyesEarly[8:56,5],rawsummaryyesEarly[57:105,5])) + .001
ax <- pretty(b:e,n=12)
start = ax[1]
shift = ax[2]-ax[1]
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_rawsummaryyesEarly_Expert <- hist(rawsummaryyesEarly[8:56,5],
                                           plot=FALSE,breaks = ax)
hg_alpha_rawsummaryyesEarly_Novice <- hist(rawsummaryyesEarly[57:105,5],plot = FALSE,breaks = ax)
plot(hg_alpha_rawsummaryyesEarly_Expert,main ='Yes Early Alpha (Main Model)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_rawsummaryyesEarly_Novice,add = TRUE,col = c2)

b <- min(c(overfullsummaryyesEarly[13:61,5],overfullsummaryyesEarly[62:110,5])) - .001
e <- max(c(overfullsummaryyesEarly[13:61,5],overfullsummaryyesEarly[62:110,5])) + .001
ax <- pretty(b:e,n=12)
start = ax[1]
shift = ax[2]-ax[1]
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_overfullsummaryyesEarly_Expert <- hist(overfullsummaryyesEarly[13:61,5],
                                                plot=FALSE,breaks = ax)
hg_alpha_overfullsummaryyesEarly_Novice <- hist(overfullsummaryyesEarly[62:110,5],plot = FALSE,breaks = ax)
plot(hg_alpha_overfullsummaryyesEarly_Expert,main ='Yes Early Alpha (12 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_overfullsummaryyesEarly_Novice,add = TRUE,col = c2)

b <- min(c(rawsummarynoEarly[8:56,5],rawsummarynoEarly[57:105,5])) - .001
e <- max(c(rawsummarynoEarly[8:56,5],rawsummarynoEarly[57:105,5])) + .001
ax <- pretty(b:e,n=12)
ax[1] = -.35
shift = .06
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_rawsummarynoEarly_Expert <- hist(rawsummarynoEarly[8:56,5],
                                          plot=FALSE,breaks = ax)
hg_alpha_rawsummarynoEarly_Novice <- hist(rawsummarynoEarly[57:105,5],plot = FALSE,breaks = ax)
plot(hg_alpha_rawsummarynoEarly_Expert,main ='No Early Alpha (Main Model)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_rawsummarynoEarly_Novice,add = TRUE,col = c2)

b <- min(c(overfullsummarynoEarly[13:61,5],overfullsummarynoEarly[62:110,5])) - .001
e <- max(c(overfullsummarynoEarly[13:61,5],overfullsummarynoEarly[62:110,5])) + .001
ax <- pretty(b:e,n=12)
start = ax[1]
shift = ax[2]-ax[1]
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_overfullsummarynoEarly_Expert <- hist(overfullsummarynoEarly[13:61,5],
                                               plot=FALSE,breaks = ax)
hg_alpha_overfullsummarynoEarly_Novice <- hist(overfullsummarynoEarly[62:110,5],plot = FALSE,breaks = ax)
plot(hg_alpha_overfullsummarynoEarly_Expert,main ='No Early Alpha (12 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_overfullsummarynoEarly_Novice,add = TRUE,col = c2)

b <- min(c(rawsummaryyesLate[8:56,5],rawsummaryyesLate[57:105,5])) - .001
e <- max(c(rawsummaryyesLate[8:56,5],rawsummaryyesLate[57:105,5])) + .001
ax <- pretty(b:e,n=12)
ax[1] = -.35
shift = .05
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_rawsummaryyesLate_Expert <- hist(rawsummaryyesLate[8:56,5],
                                          plot=FALSE,breaks = ax)
hg_alpha_rawsummaryyesLate_Novice <- hist(rawsummaryyesLate[57:105,5],plot = FALSE,breaks = ax)
plot(hg_alpha_rawsummaryyesLate_Expert,main ='Yes Late Alpha (Main Model)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_rawsummaryyesLate_Novice,add = TRUE,col = c2)

b <- min(c(overfullsummaryyesLate[13:61,5],overfullsummaryyesLate[62:110,5])) - .001
e <- max(c(overfullsummaryyesLate[13:61,5],overfullsummaryyesLate[62:110,5])) + .001
ax <- pretty(b:e,n=12)
ax[1] =-.35
shift = .05
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_overfullsummaryyesLate_Expert <- hist(overfullsummaryyesLate[13:61,5],
                                               plot=FALSE,breaks = ax)
hg_alpha_overfullsummaryyesLate_Novice <- hist(overfullsummaryyesLate[62:110,5],plot = FALSE,breaks = ax)
plot(hg_alpha_overfullsummaryyesLate_Expert,main ='Yes Late Alpha (12 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_overfullsummaryyesLate_Novice,add = TRUE,col = c2)

b <- min(c(rawsummarynoLate[8:56,5],rawsummarynoLate[57:105,5])) - .001
e <- max(c(rawsummarynoLate[8:56,5],rawsummarynoLate[57:105,5])) + .001
ax <- pretty(b:e,n=12)
ax[1] = 0
shift = .07
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_rawsummarynoLate_Expert <- hist(rawsummarynoLate[8:56,5],
                                         plot=FALSE,breaks = ax)
hg_alpha_rawsummarynoLate_Novice <- hist(rawsummarynoLate[57:105,5],plot = FALSE,breaks = ax)
plot(hg_alpha_rawsummarynoLate_Expert,main ='No Late Alpha (Main Model)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_rawsummarynoLate_Novice,add = TRUE,col = c2)

b <- min(c(overfullsummarynoLate[13:61,5],overfullsummarynoLate[62:110,5])) - .001
e <- max(c(overfullsummarynoLate[13:61,5],overfullsummarynoLate[62:110,5])) + .001
ax <- pretty(b:e,n=12)
ax[1] = -.1
shift = .07
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_overfullsummarynoLate_Expert <- hist(overfullsummarynoLate[13:61,5],
                                              plot=FALSE,breaks = ax)
hg_alpha_overfullsummarynoLate_Novice <- hist(overfullsummarynoLate[62:110,5],plot = FALSE,breaks = ax)
plot(hg_alpha_overfullsummarynoLate_Expert,main ='No Late Alpha (12 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_overfullsummarynoLate_Novice,add = TRUE,col = c2)


# Expected order test for experts in "Yes" condition, overfull model ----------------------

OrderTestEarly_expertoverfull <- matrix(NaN,nrow=nrow(overfullBetayesEarlysamples),ncol = 7)

OrderTestEarly_expertoverfull[,1:6] <- overfullBetayesEarlysamples[,1:6]

##Early hypothesis: 1 > 2 > 3, 6 > 5 > 4
for(i in 1:nrow(overfullBetayesEarlysamples)){
  firstcheck <- OrderTestEarly_expertoverfull[i,1] > OrderTestEarly_expertoverfull[i,2] && 
    OrderTestEarly_expertoverfull[i,2] > OrderTestEarly_expertoverfull[i,3]
  secondcheck <- OrderTestEarly_expertoverfull[i,6] > OrderTestEarly_expertoverfull[i,5] && 
    OrderTestEarly_expertoverfull[i,5] > OrderTestEarly_expertoverfull[i,4]
  OrderTestEarly_expertoverfull[i,7] <- firstcheck && secondcheck
}

OrderTestLate_expertoverfull <- matrix(NaN,nrow=nrow(overfullBetayesLatesamples),ncol = 7)

OrderTestLate_expertoverfull[,1:6] <- overfullBetayesLatesamples[,1:6]

##Late hypothesis: 3 > 2 > 1, 4 > 5 > 6
for(i in 1:nrow(overfullBetayesLatesamples)){
  firstcheck <- OrderTestLate_expertoverfull[i,3] > OrderTestLate_expertoverfull[i,2] && 
    OrderTestLate_expertoverfull[i,2] > OrderTestLate_expertoverfull[i,1]
  secondcheck <- OrderTestLate_expertoverfull[i,4] > OrderTestLate_expertoverfull[i,5] && 
    OrderTestLate_expertoverfull[i,5] > OrderTestLate_expertoverfull[i,6]
  OrderTestLate_expertoverfull[i,7] <- firstcheck && secondcheck
}

EarlyOrderProb_expertoverfull <- mean(OrderTestEarly_expertoverfull[,7])
LateOrderProb_expertoverfull <- mean(OrderTestLate_expertoverfull[,7])

EarlyOrderOdds_expertoverfull <- EarlyOrderProb_expertoverfull/(1-EarlyOrderProb_expertoverfull)
LateOrderOdds_expertoverfull <- LateOrderProb_expertoverfull/(1-LateOrderProb_expertoverfull)

PriorOdds_expertoverfull <- (1/36)/(1-(1/36))

BFEarlyOrder_expertoverfull <- EarlyOrderOdds_expertoverfull/PriorOdds_expertoverfull
BFLateOrder_expertoverfull <- LateOrderOdds_expertoverfull/PriorOdds_expertoverfull

# Expected order test for novices in "Yes" condition, overfull model ----------------------

OrderTestEarly_noviceoverfull <- matrix(NaN,nrow=nrow(overfullBetayesEarlysamples),ncol = 7)

OrderTestEarly_noviceoverfull[,1:6] <- overfullBetayesEarlysamples[,7:12]

##Early hypothesis: 1 > 2 > 3, 6 > 5 > 4
for(i in 1:nrow(overfullBetayesEarlysamples)){
  firstcheck <- OrderTestEarly_noviceoverfull[i,1] > OrderTestEarly_noviceoverfull[i,2] && 
    OrderTestEarly_noviceoverfull[i,2] > OrderTestEarly_noviceoverfull[i,3]
  secondcheck <- OrderTestEarly_noviceoverfull[i,6] > OrderTestEarly_noviceoverfull[i,5] && 
    OrderTestEarly_noviceoverfull[i,5] > OrderTestEarly_noviceoverfull[i,4]
  OrderTestEarly_noviceoverfull[i,7] <- firstcheck && secondcheck
}

OrderTestLate_noviceoverfull <- matrix(NaN,nrow=nrow(overfullBetayesLatesamples),ncol = 7)

OrderTestLate_noviceoverfull[,1:6] <- overfullBetayesLatesamples[,7:12]

##Late hypothesis: 3 > 2 > 1, 4 > 5 > 6
for(i in 1:nrow(overfullBetayesLatesamples)){
  firstcheck <- OrderTestLate_noviceoverfull[i,3] > OrderTestLate_noviceoverfull[i,2] && 
    OrderTestLate_noviceoverfull[i,2] > OrderTestLate_noviceoverfull[i,1]
  secondcheck <- OrderTestLate_noviceoverfull[i,4] > OrderTestLate_noviceoverfull[i,5] && 
    OrderTestLate_noviceoverfull[i,5] > OrderTestLate_noviceoverfull[i,6]
  OrderTestLate_noviceoverfull[i,7] <- firstcheck && secondcheck
}

EarlyOrderProb_noviceoverfull <- mean(OrderTestEarly_noviceoverfull[,7])
LateOrderProb_noviceoverfull <- mean(OrderTestLate_noviceoverfull[,7])

EarlyOrderOdds_noviceoverfull <- EarlyOrderProb_noviceoverfull/(1-EarlyOrderProb_noviceoverfull)
LateOrderOdds_noviceoverfull <- LateOrderProb_noviceoverfull/(1-LateOrderProb_noviceoverfull)

PriorOdds_noviceoverfull <- (1/36)/(1-(1/36))

BFEarlyOrder_noviceoverfull <- EarlyOrderOdds_noviceoverfull/PriorOdds_noviceoverfull
BFLateOrder_noviceoverfull <- LateOrderOdds_noviceoverfull/PriorOdds_noviceoverfull

# Prep data for sparse JAGS inputs -----------------------------------------------
nParticipants = nExperts+nNovices
nBeta = 2
nTrials = 6
nExpResponses = nExperts*nTrials
nNovResponses = nNovices*nTrials

XyesEarly_sparse <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YyesEarly <- rep(NaN,nParticipants*nTrials)

XnoEarly_sparse <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YnoEarly <- rep(NaN,nParticipants*nTrials)

XyesLate_sparse <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YyesLate <- rep(NaN,nParticipants*nTrials)

XnoLate_sparse <- matrix(data=rep(0,nParticipants*nTrials*nBeta),nrow = nParticipants*nTrials,ncol = nBeta)
YnoLate <- rep(NaN,nParticipants*nTrials)

tickeryesEarly = 1
tickernoEarly = 1
tickeryesLate = 1
tickernoLate = 1

XyesEarly_sparse[1:nExpResponses,1] <- 1
XyesEarly_sparse[(nExpResponses+1):(nParticipants*nTrials),2] <- 1

for(i in 1:nrow(OnlineData)){
  
  
  if(OnlineData$SpeakerEarly[i]==1){
    
    YyesEarly[tickeryesEarly] <- OnlineData$EarlyResponse[i]
    
    tickeryesEarly = tickeryesEarly+1
  }
  
  if(OnlineData$SpeakerEarly[i]==0){
    
    YnoEarly[tickernoEarly] <- OnlineData$EarlyResponse[i]
    
    tickernoEarly = tickernoEarly+1
  }
  
  if(OnlineData$SpeakerLate[i]==1){
    
    YyesLate[tickeryesLate] <- OnlineData$LateResponse[i]
    
    tickeryesLate = tickeryesLate+1
  }
  
  if(OnlineData$SpeakerLate[i]==0){
    
    YnoLate[tickernoLate] <- OnlineData$LateResponse[i]
    
    tickernoLate = tickernoLate+1
  }
}

for(i in 1:nrow(SONAData)){
  
  
  if(SONAData$SpeakerEarly[i]==1){
    
    YyesEarly[tickeryesEarly] <- SONAData$EarlyResponse[i]
    
    tickeryesEarly = tickeryesEarly+1
  }
  
  if(SONAData$SpeakerEarly[i]==0){
    
    YnoEarly[tickernoEarly] <- SONAData$EarlyResponse[i]
    
    tickernoEarly = tickernoEarly+1
  }
  
  if(SONAData$SpeakerLate[i]==1){
    
    YyesLate[tickeryesLate] <- SONAData$LateResponse[i]
    
    tickeryesLate = tickeryesLate+1
  }
  
  if(SONAData$SpeakerLate[i]==0){
    
    YnoLate[tickernoLate] <- SONAData$LateResponse[i]
    
    tickernoLate = tickernoLate+1
  }
}


# MVN Precision Matrix (for sparse JAGS betas) ----------------------------------------
bPrec <- matrix(0,nrow = nBeta, ncol = nBeta)

for(i in 1:nBeta){
  bPrec[i,i] <- 16
}

# Run sparse linear regression in JAGS -------------------------------------------
nSamples = 4000

sparseBetayesEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))
sparseDIC <- rep(NaN,4)

X <- XyesEarly_sparse
Y<- YyesEarly

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

sparseBetayesEarlysamples <- samples$BUGSoutput$sims.list$BETA

sparsesummaryyesEarly <- samples$BUGSoutput$summary

sparseDIC[1] <- samples$BUGSoutput$DIC

##No Early Regression##
sparseBetanoEarlysamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoEarly_sparse
Y<- YnoEarly

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

sparseBetanoEarlysamples <- samples$BUGSoutput$sims.list$BETA

sparsesummarynoEarly <- samples$BUGSoutput$summary

sparseDIC[2] <- samples$BUGSoutput$DIC

##Yes Late Regression##

sparseBetayesLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XyesLate_sparse
Y<- YyesLate

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

sparseBetayesLatesamples <- samples$BUGSoutput$sims.list$BETA

sparsesummaryyesLate <- samples$BUGSoutput$summary

sparseDIC[3] <- samples$BUGSoutput$DIC

##No Late Regression##
sparseBetanoLatesamples <- array(rep(NaN,nSamples*nBeta),c(nSamples,nBeta))

X <- XnoLate_sparse
Y<- YnoLate

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

sparseBetanoLatesamples <- samples$BUGSoutput$sims.list$BETA

sparsesummarynoLate <- samples$BUGSoutput$summary

sparseDIC[4] <- samples$BUGSoutput$DIC

# Plot beta inferences from sparse -------------------------------------------
sparsePlotValues_EarlyYes <- matrix(NaN, nrow = 2,ncol = 3)
sparsePlotValues_EarlyNo <- matrix(NaN, nrow = 2,ncol = 3)
sparsePlotValues_LateYes <- matrix(NaN, nrow = 2,ncol = 3)
sparsePlotValues_LateNo <- matrix(NaN, nrow = 2,ncol = 3)

loconf = nrow(sparseBetayesEarlysamples) * .025
hiconf = nrow(sparseBetayesEarlysamples) * .975


tempsort <- sort(sparseBetayesEarlysamples[,1])
sparsePlotValues_EarlyYes[1,1] <- mean(tempsort[1500:1501])
sparsePlotValues_EarlyYes[1,2] <- tempsort[loconf]
sparsePlotValues_EarlyYes[1,3] <- tempsort[hiconf]

tempsort <- sort(sparseBetanoEarlysamples[,1])
sparsePlotValues_EarlyNo[1,1] <- mean(tempsort[1500:1501])
sparsePlotValues_EarlyNo[1,2] <- tempsort[loconf]
sparsePlotValues_EarlyNo[1,3] <- tempsort[hiconf]

tempsort <- sort(sparseBetayesLatesamples[,1])
sparsePlotValues_LateYes[1,1] <- mean(tempsort[1500:1501])
sparsePlotValues_LateYes[1,2] <- tempsort[loconf]
sparsePlotValues_LateYes[1,3] <- tempsort[hiconf]

tempsort <- sort(sparseBetanoLatesamples[,1])
sparsePlotValues_LateNo[1,1] <- mean(tempsort[1500:1501])
sparsePlotValues_LateNo[1,2] <- tempsort[loconf]
sparsePlotValues_LateNo[1,3] <- tempsort[hiconf]

tempsort <- sort(sparseBetayesEarlysamples[,2])
sparsePlotValues_EarlyYes[2,1] <- mean(tempsort[1500:1501])
sparsePlotValues_EarlyYes[2,2] <- tempsort[loconf]
sparsePlotValues_EarlyYes[2,3] <- tempsort[hiconf]

tempsort <- sort(sparseBetanoEarlysamples[,2])
sparsePlotValues_EarlyNo[2,1] <- mean(tempsort[1500:1501])
sparsePlotValues_EarlyNo[2,2] <- tempsort[loconf]
sparsePlotValues_EarlyNo[2,3] <- tempsort[hiconf]

tempsort <- sort(sparseBetayesLatesamples[,2])
sparsePlotValues_LateYes[2,1] <- mean(tempsort[1500:1501])
sparsePlotValues_LateYes[2,2] <- tempsort[loconf]
sparsePlotValues_LateYes[2,3] <- tempsort[hiconf]

tempsort <- sort(sparseBetanoLatesamples[,2])
sparsePlotValues_LateNo[2,1] <- mean(tempsort[1500:1501])
sparsePlotValues_LateNo[2,2] <- tempsort[loconf]
sparsePlotValues_LateNo[2,3] <- tempsort[hiconf]


plot(1,sparsePlotValues_EarlyYes[1,1],xlim = c(.5,2.5),ylim = c(0,1),main = "Yes Early (Only 2 Betas)",xlab = 'beta',ylab = 'Interpretation')

arrows(x0 = 1, y0 = sparsePlotValues_EarlyYes[1,2],
       x1 = 1, y1 = sparsePlotValues_EarlyYes[1,3],code = 3,
       col = 'black',angle = 90, length = .1)

points(2,sparsePlotValues_EarlyYes[2,1],col = '#00CC00',pch = 2)

arrows(x0 = 2, y0 = sparsePlotValues_EarlyYes[2,2],
       x1 = 2, y1 = sparsePlotValues_EarlyYes[2,3],code = 3,
       col = '#00CC00',angle = 90, length = .1,lty = 3)

plot(1,sparsePlotValues_EarlyNo[1,1],xlim = c(.5,2.5),ylim = c(0,1),main = "No Early (Only 2 Betas)",xlab = 'beta',ylab = 'Interpretation')

arrows(x0 = 1, y0 = sparsePlotValues_EarlyNo[1,2],
       x1 = 1, y1 = sparsePlotValues_EarlyNo[1,3],code = 3,
       col = 'black',angle = 90, length = .1)

points(2,sparsePlotValues_EarlyNo[2,1],col = '#00CC00',pch = 2)


arrows(x0 = 2, y0 = sparsePlotValues_EarlyNo[2,2],
       x1 = 2, y1 = sparsePlotValues_EarlyNo[2,3],code = 3,
       col = '#00CC00',angle = 90, length = .1,lty = 3)

plot(1,sparsePlotValues_LateYes[1,1],xlim = c(.5,2.5),ylim = c(0,1),main = "Yes Late (Only 2 Betas)",xlab = 'beta',ylab = 'Interpretation')

arrows(x0 = 1, y0 = sparsePlotValues_LateYes[1,2],
       x1 = 1, y1 = sparsePlotValues_LateYes[1,3],code = 3,
       col = 'black',angle = 90, length = .1)

points(2,sparsePlotValues_LateYes[2,1],col = '#00CC00',pch = 2)

arrows(x0 = 2, y0 = sparsePlotValues_LateYes[2,2],
       x1 = 2, y1 = sparsePlotValues_LateYes[2,3],code = 3,
       col = '#00CC00',angle = 90, length = .1,lty = 3)

plot(1,sparsePlotValues_LateNo[1,1],xlim = c(.5,2.5),ylim = c(0,1),main = "No Late (Only 2 Betas)",xlab = 'beta',ylab = 'Interpretation')

arrows(x0 = 1, y0 = sparsePlotValues_LateNo[1,2],
       x1 = 1, y1 = sparsePlotValues_LateNo[1,3],code = 3,
       col = 'black',angle = 90, length = .1)

points(2,sparsePlotValues_LateNo[2,1],col = '#00CC00',pch = 2)


arrows(x0 = 2, y0 = sparsePlotValues_LateNo[2,2],
       x1 = 2, y1 = sparsePlotValues_LateNo[2,3],code = 3,
       col = '#00CC00',angle = 90, length = .1,lty = 3)

# Compare DIC's (between main and sparse) -----------------------------------------------------------

par(mfrow=c(1,1))

plot(1:4,rawDIC,ylim = c(-650,-200),xlab = 'Condition',ylab = 'DIC', main = "Main Model and Sparse DIC Comparison")

points(1:4,sparseDIC,col = '#DC143C',pch=3)

legend(3, -200, legend=c("Main Model", "Sparse"),
       col=c("black", "red"), pch=c(1,3), cex=0.8)

text_rawDIC <- matrix(NaN,nrow = 4, ncol =2)
text_rawDIC[,1] <- 1:4
text_rawDIC[,2] <- round(rawDIC,0)
text_rawDIC <- data.frame(text_rawDIC)
names(text_rawDIC) <- c("Condition","DIC")

text(DIC~Condition, labels=DIC,data=text_rawDIC, cex=0.9, font=1,pos=1)

text_sparseDIC <- matrix(NaN,nrow = 4, ncol =2)
text_sparseDIC[,1] <- 1:4
text_sparseDIC[,2] <- round(sparseDIC,0)
text_sparseDIC <- data.frame(text_sparseDIC)
names(text_sparseDIC) <- c("Condition","DIC")

text(DIC~Condition, labels=DIC,data=text_sparseDIC, cex=0.9, font=1,pos=3,col='#DC143C')


# Compare alphas (between main and sparse) ------------------------------

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

par(mfrow=c(2,2))
# 
# b <- min(c(rawsummaryyesEarly[8:56,5],rawsummaryyesEarly[57:105,5])) - .001
# e <- max(c(rawsummaryyesEarly[8:56,5],rawsummaryyesEarly[57:105,5])) + .001
# ax <- pretty(b:e,n=12)
# start = ax[1]
# shift = ax[2]-ax[1]
# for(i in 2:14){
#   ax[i]<-ax[i-1]+shift
# }
# hg_alpha_rawsummaryyesEarly_Expert <- hist(rawsummaryyesEarly[8:56,5],
#                                            plot=FALSE,breaks = ax)
# hg_alpha_rawsummaryyesEarly_Novice <- hist(rawsummaryyesEarly[57:105,5],plot = FALSE,breaks = ax)
# plot(hg_alpha_rawsummaryyesEarly_Expert,main ='Yes Early Alpha (Main Model)',xlim=c(-1,1),col=c1,
#      xlab = 'alpha estimate')
# plot(hg_alpha_rawsummaryyesEarly_Novice,add = TRUE,col = c2)


b <- min(c(sparsesummaryyesEarly[3:51,5],sparsesummaryyesEarly[52:100,5])) - .001
e <- max(c(sparsesummaryyesEarly[3:51,5],sparsesummaryyesEarly[52:100,5])) + .001
ax <- pretty(b:e,n=12)
start = ax[1]
shift = ax[2]-ax[1]
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_sparsesummaryyesEarly_Expert <- hist(sparsesummaryyesEarly[3:51,5],
                                              plot=FALSE,breaks=ax)
hg_alpha_sparsesummaryyesEarly_Novice <- hist(sparsesummaryyesEarly[52:100,5],plot = FALSE,breaks=ax)
plot(hg_alpha_sparsesummaryyesEarly_Expert,main ='Yes Early Alpha (Only 2 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_sparsesummaryyesEarly_Novice,add = TRUE,col = c2)

# b <- min(c(rawsummarynoEarly[8:56,5],rawsummarynoEarly[57:105,5])) - .001
# e <- max(c(rawsummarynoEarly[8:56,5],rawsummarynoEarly[57:105,5])) + .001
# ax <- pretty(b:e,n=12)
# ax[1] = -.35
# shift = .06
# for(i in 2:14){
#   ax[i]<-ax[i-1]+shift
# }
# hg_alpha_rawsummarynoEarly_Expert <- hist(rawsummarynoEarly[8:56,5],
#                                           plot=FALSE,breaks = ax)
# hg_alpha_rawsummarynoEarly_Novice <- hist(rawsummarynoEarly[57:105,5],plot = FALSE,breaks = ax)
# plot(hg_alpha_rawsummarynoEarly_Expert,main ='No Early Alpha (Main Model)',xlim=c(-1,1),col=c1,
#      xlab = 'alpha estimate')
# plot(hg_alpha_rawsummarynoEarly_Novice,add = TRUE,col = c2)

b <- min(c(sparsesummarynoEarly[3:51,5],sparsesummarynoEarly[52:100,5])) - .001
e <- max(c(sparsesummarynoEarly[3:51,5],sparsesummarynoEarly[52:100,5])) + .001
ax <- pretty(b:e,n=12)
start = 0
shift = .08
ax[1] <- start
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_sparsesummarynoEarly_Expert <- hist(sparsesummarynoEarly[3:51,5],
                                             plot=FALSE,breaks = ax)
hg_alpha_sparsesummarynoEarly_Novice <- hist(sparsesummarynoEarly[52:100,5],plot = FALSE,breaks = ax)
plot(hg_alpha_sparsesummarynoEarly_Expert,main ='No Early Alpha (Only 2 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_sparsesummarynoEarly_Novice,add = TRUE,col = c2)

# b <- min(c(rawsummaryyesLate[8:56,5],rawsummaryyesLate[57:105,5])) - .001
# e <- max(c(rawsummaryyesLate[8:56,5],rawsummaryyesLate[57:105,5])) + .001
# ax <- pretty(b:e,n=12)
# start = ax[1]
# shift = .06
# for(i in 2:14){
#   ax[i]<-ax[i-1]+shift
# }
# hg_alpha_rawsummaryyesLate_Expert <- hist(rawsummaryyesLate[8:56,5],
#                                           plot=FALSE,breaks = ax)
# hg_alpha_rawsummaryyesLate_Novice <- hist(rawsummaryyesLate[57:105,5],plot = FALSE,breaks = ax)
# plot(hg_alpha_rawsummaryyesLate_Expert,main ='Yes Late Alpha (Main Model)',xlim=c(-1,1),col=c1,
#      xlab = 'alpha estimate')
# plot(hg_alpha_rawsummaryyesLate_Novice,add = TRUE,col = c2)

b <- min(c(sparsesummaryyesLate[3:51,5],sparsesummaryyesLate[52:100,5])) - .001
e <- max(c(sparsesummaryyesLate[3:51,5],sparsesummaryyesLate[52:100,5])) + .001
ax <- pretty(b:e,n=12)
ax[1] = .25
shift = .06
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_sparsesummaryyesLate_Expert <- hist(sparsesummaryyesLate[3:51,5],
                                             plot=FALSE,breaks = ax)
hg_alpha_sparsesummaryyesLate_Novice <- hist(sparsesummaryyesLate[52:100,5],plot = FALSE,breaks = ax)
plot(hg_alpha_sparsesummaryyesLate_Expert,main ='Yes Late Alpha (Only 2 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_sparsesummaryyesLate_Novice,add = TRUE,col = c2)

# b <- min(c(rawsummarynoLate[3:51,5],rawsummarynoLate[52:100,5])) - .001
# e <- max(c(sparsesummarynoLate[3:51,5],sparsesummarynoLate[52:100,5])) + .001
# ax <- pretty(b:e,n=12)
# ax[1] = -.4
# shift = .07
# for(i in 2:14){
#   ax[i]<-ax[i-1]+shift
# }
# hg_alpha_rawsummarynoLate_Expert <- hist(rawsummarynoLate[8:56,5],
#                                          plot=FALSE,breaks = ax)
# hg_alpha_rawsummarynoLate_Novice <- hist(rawsummarynoLate[57:105,5],plot = FALSE,breaks = ax)
# plot(hg_alpha_rawsummarynoLate_Expert,main ='No Late Alpha (Main Model)',xlim=c(-1,1),col=c1,
#      xlab = 'alpha estimate')
# plot(hg_alpha_rawsummarynoLate_Novice,add = TRUE,col = c2)

b <- min(c(sparsesummarynoLate[3:51,5],sparsesummarynoLate[52:100,5])) - .001
e <- max(c(sparsesummarynoLate[3:51,5],sparsesummarynoLate[52:100,5])) + .001
ax <- pretty(b:e,n=12)
ax[1] = 0
shift = .07
for(i in 2:14){
  ax[i]<-ax[i-1]+shift
}
hg_alpha_sparsesummarynoLate_Expert <- hist(sparsesummarynoLate[3:51,5],
                                            plot=FALSE,breaks = ax)
hg_alpha_sparsesummarynoLate_Novice <- hist(sparsesummarynoLate[52:100,5],plot = FALSE,breaks = ax)
plot(hg_alpha_sparsesummarynoLate_Expert,main ='No Late Alpha (Only 2 Betas)',xlim=c(-1,1),col=c1,
     xlab = 'alpha estimate')
plot(hg_alpha_sparsesummarynoLate_Novice,add = TRUE,col = c2)


# Test model predictions w/ epsilon--------------------------------------------------
nBeta = 7

postpred_lo = nBeta+1
postpred_mid = postpred_lo+nExpResponses-1
postpred_hi = postpred_mid+nNovResponses

par(mfrow=c(2,2))
plot(YyesEarly[1:nExpResponses],rawsummaryyesEarly[postpred_lo:postpred_mid,1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'Yes Early',col = '#7B68EE')
points(YyesEarly[(nExpResponses+1):length(YyesEarly)],rawsummaryyesEarly[(postpred_mid+1):postpred_hi,1],
       col = '#DC143C')
abline(a=0,b=1)

legend(.6, .2, legend=c("Expert", "Novice"),
       col=c("blue", "red"), pch=c(1,1), cex=0.8)

plot(YyesLate[1:nExpResponses],rawsummaryyesLate[postpred_lo:postpred_mid,1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'Yes Late',col = '#7B68EE')
points(YyesLate[(nExpResponses+1):length(YyesLate)],rawsummaryyesLate[(postpred_mid+1):postpred_hi,1],
       col = '#DC143C')
abline(a=0,b=1)

plot(YnoEarly[1:nExpResponses],rawsummarynoEarly[postpred_lo:postpred_mid,1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'No Early',col = '#7B68EE')
points(YnoEarly[(nExpResponses+1):length(YnoEarly)],rawsummarynoEarly[(postpred_mid+1):postpred_hi,1],
       col = '#DC143C')
abline(a=0,b=1)

plot(YnoLate[1:nExpResponses],rawsummarynoLate[postpred_lo:postpred_mid,1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'No Late',col = '#7B68EE')
points(YnoLate[(nExpResponses+1):length(YnoLate)],rawsummarynoLate[(postpred_mid+1):postpred_hi,1],
       col = '#DC143C')
abline(a=0,b=1)


# Test model predictions w/o epsilon--------------------------------------------------
betas_yesEarly <- c(rawsummaryyesEarly[1:7,1])
betas_noEarly <- c(rawsummarynoEarly[1:7,1])
betas_yesLate <- c(rawsummaryyesLate[1:7,1])
betas_noLate <- c(rawsummarynoLate[1:7,1])

testTotal = nExpResponses+nNovResponses

raw_ModelPredictions_yesEarly <- matrix(NaN,nrow = testTotal,ncol = 4)
raw_ModelPredictions_noEarly <- matrix(NaN,nrow = testTotal,ncol = 4)
raw_ModelPredictions_yesLate <- matrix(NaN,nrow = testTotal,ncol = 4)
raw_ModelPredictions_noLate <- matrix(NaN,nrow = testTotal,ncol = 4)

raw_ModelPredictions_yesEarly[1:nExpResponses,1] <- OnlineData$Participant[which(OnlineData$SpeakerEarly==1)]
raw_ModelPredictions_yesEarly[1:nExpResponses,2] <- OnlineData$CompType[which(OnlineData$SpeakerEarly==1)]
raw_ModelPredictions_yesEarly[1:nExpResponses,3] <- OnlineData$EarlyResponse[which(OnlineData$SpeakerEarly==1)]
for(i in 1:nExpResponses){
  raw_ModelPredictions_yesEarly[i,4] <- betas_yesEarly[1+raw_ModelPredictions_yesEarly[i,2]]+
    rawsummaryyesEarly[nBeta+testTotal+raw_ModelPredictions_yesEarly[i,1],1]
}

raw_ModelPredictions_yesEarly[(nExpResponses+1):testTotal,1] <- SONAData$Participant[which(SONAData$SpeakerEarly==1)]
raw_ModelPredictions_yesEarly[(nExpResponses+1):testTotal,2] <- SONAData$CompType[which(SONAData$SpeakerEarly==1)]
raw_ModelPredictions_yesEarly[(nExpResponses+1):testTotal,3] <- SONAData$EarlyResponse[which(SONAData$SpeakerEarly==1)]
for(i in (nExpResponses+1):testTotal){
  raw_ModelPredictions_yesEarly[i,4] <- betas_yesEarly[1]+
    rawsummaryyesEarly[nBeta+testTotal+raw_ModelPredictions_yesEarly[i,1]+nExperts,1]
}
par(mfrow=c(2,2))
plot(raw_ModelPredictions_yesEarly[1:nExpResponses,3],raw_ModelPredictions_yesEarly[1:nExpResponses,4],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'Yes Early',col = '#7B68EE')
points(raw_ModelPredictions_yesEarly[(nExpResponses+1):length(YyesEarly),3],raw_ModelPredictions_yesEarly[(nExpResponses+1):length(YyesEarly),4],
       col = '#DC143C')
abline(a=0,b=1)

legend(.6, .2, legend=c("Expert", "Novice"),
       col=c("blue", "red"), pch=c(1,1), cex=0.8)

raw_ModelPredictions_noEarly[1:nExpResponses,1] <- OnlineData$Participant[which(OnlineData$SpeakerEarly==0)]
raw_ModelPredictions_noEarly[1:nExpResponses,2] <- OnlineData$CompType[which(OnlineData$SpeakerEarly==0)]
raw_ModelPredictions_noEarly[1:nExpResponses,3] <- OnlineData$EarlyResponse[which(OnlineData$SpeakerEarly==0)]
for(i in 1:nExpResponses){
  raw_ModelPredictions_noEarly[i,4] <- betas_noEarly[1+raw_ModelPredictions_noEarly[i,2]]+
    rawsummarynoEarly[nBeta+testTotal+raw_ModelPredictions_noEarly[i,1],1]
}

raw_ModelPredictions_noEarly[(nExpResponses+1):testTotal,1] <- SONAData$Participant[which(SONAData$SpeakerEarly==0)]
raw_ModelPredictions_noEarly[(nExpResponses+1):testTotal,2] <- SONAData$CompType[which(SONAData$SpeakerEarly==0)]
raw_ModelPredictions_noEarly[(nExpResponses+1):testTotal,3] <- SONAData$EarlyResponse[which(SONAData$SpeakerEarly==0)]
for(i in (nExpResponses+1):testTotal){
  raw_ModelPredictions_noEarly[i,4] <- betas_noEarly[1]+
    rawsummarynoEarly[nBeta+testTotal+raw_ModelPredictions_noEarly[i,1]+nExperts,1]
}

plot(raw_ModelPredictions_noEarly[1:nExpResponses,3],raw_ModelPredictions_noEarly[1:nExpResponses,4],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'No Early',col = '#7B68EE')
points(raw_ModelPredictions_noEarly[(nExpResponses+1):length(YnoEarly),3],raw_ModelPredictions_noEarly[(nExpResponses+1):length(YnoEarly),4],
       col = '#DC143C')
abline(a=0,b=1)

raw_ModelPredictions_yesLate[1:nExpResponses,1] <- OnlineData$Participant[which(OnlineData$SpeakerLate==1)]
raw_ModelPredictions_yesLate[1:nExpResponses,2] <- OnlineData$CompType[which(OnlineData$SpeakerLate==1)]
raw_ModelPredictions_yesLate[1:nExpResponses,3] <- OnlineData$LateResponse[which(OnlineData$SpeakerLate==1)]
for(i in 1:nExpResponses){
  raw_ModelPredictions_yesLate[i,4] <- betas_yesLate[1+raw_ModelPredictions_yesLate[i,2]]+
    rawsummaryyesLate[nBeta+testTotal+raw_ModelPredictions_yesLate[i,1],1]
}

raw_ModelPredictions_yesLate[(nExpResponses+1):testTotal,1] <- SONAData$Participant[which(SONAData$SpeakerLate==1)]
raw_ModelPredictions_yesLate[(nExpResponses+1):testTotal,2] <- SONAData$CompType[which(SONAData$SpeakerLate==1)]
raw_ModelPredictions_yesLate[(nExpResponses+1):testTotal,3] <- SONAData$LateResponse[which(SONAData$SpeakerLate==1)]
for(i in (nExpResponses+1):testTotal){
  raw_ModelPredictions_yesLate[i,4] <- betas_yesLate[1]+
    rawsummaryyesLate[nBeta+testTotal+raw_ModelPredictions_yesLate[i,1]+nExperts,1]
}

plot(raw_ModelPredictions_yesLate[1:nExpResponses,3],raw_ModelPredictions_yesLate[1:nExpResponses,4],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'Yes Late',col = '#7B68EE')
points(raw_ModelPredictions_yesLate[(nExpResponses+1):length(YyesLate),3],raw_ModelPredictions_yesLate[(nExpResponses+1):length(YyesLate),4],
       col = '#DC143C')
abline(a=0,b=1)

raw_ModelPredictions_noLate[1:nExpResponses,1] <- OnlineData$Participant[which(OnlineData$SpeakerLate==0)]
raw_ModelPredictions_noLate[1:nExpResponses,2] <- OnlineData$CompType[which(OnlineData$SpeakerLate==0)]
raw_ModelPredictions_noLate[1:nExpResponses,3] <- OnlineData$LateResponse[which(OnlineData$SpeakerLate==0)]
for(i in 1:nExpResponses){
  raw_ModelPredictions_noLate[i,4] <- betas_noLate[1+raw_ModelPredictions_noLate[i,2]]+
    rawsummarynoLate[nBeta+testTotal+raw_ModelPredictions_noLate[i,1],1]
}

raw_ModelPredictions_noLate[(nExpResponses+1):testTotal,1] <- SONAData$Participant[which(SONAData$SpeakerLate==0)]
raw_ModelPredictions_noLate[(nExpResponses+1):testTotal,2] <- SONAData$CompType[which(SONAData$SpeakerLate==0)]
raw_ModelPredictions_noLate[(nExpResponses+1):testTotal,3] <- SONAData$LateResponse[which(SONAData$SpeakerLate==0)]
for(i in (nExpResponses+1):testTotal){
  raw_ModelPredictions_noLate[i,4] <- betas_noLate[1]+
    rawsummarynoLate[nBeta+testTotal+raw_ModelPredictions_noLate[i,1]+nExperts,1]
}

plot(raw_ModelPredictions_noLate[1:nExpResponses,3],raw_ModelPredictions_noLate[1:nExpResponses,4],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'No Late',col = '#7B68EE')
points(raw_ModelPredictions_noLate[(nExpResponses+1):length(YnoLate),3],raw_ModelPredictions_noLate[(nExpResponses+1):length(YnoLate),4],
       col = '#DC143C')
abline(a=0,b=1)

# Posterior predictive accuracy by comp type ------------------------------
yesEarly_postpredbycomp <- matrix(NaN,nrow=4,ncol=6)
noEarly_postpredbycomp <- matrix(NaN,nrow=4,ncol=6)
yesLate_postpredbycomp <- matrix(NaN,nrow=4,ncol=6)
noLate_postpredbycomp <- matrix(NaN,nrow=4,ncol=6)

for(i in 1:6){
  tempindex = which(XyesEarly[,i+1]==1)
  tempindex_yPred = tempindex+nBeta
  yesEarly_postpredbycomp[1,i] = mean(YyesEarly[tempindex])
  yesEarly_postpredbycomp[2,i] = sd(YyesEarly[tempindex])
  yesEarly_postpredbycomp[3,i] = mean(rawsummaryyesEarly[tempindex_yPred,5])
  yesEarly_postpredbycomp[4,i] = sd(rawsummaryyesEarly[tempindex_yPred,5])
}

for(i in 1:6){
  tempindex = which(XnoEarly[,i+1]==1)
  tempindex_yPred = tempindex+nBeta
  noEarly_postpredbycomp[1,i] = mean(YnoEarly[tempindex])
  noEarly_postpredbycomp[2,i] = sd(YnoEarly[tempindex])
  noEarly_postpredbycomp[3,i] = mean(rawsummarynoEarly[tempindex_yPred,5])
  noEarly_postpredbycomp[4,i] = sd(rawsummarynoEarly[tempindex_yPred,5])
}

for(i in 1:6){
  tempindex = which(XyesLate[,i+1]==1)
  tempindex_yPred = tempindex+nBeta
  yesLate_postpredbycomp[1,i] = mean(YyesLate[tempindex])
  yesLate_postpredbycomp[2,i] = sd(YyesLate[tempindex])
  yesLate_postpredbycomp[3,i] = mean(rawsummaryyesLate[tempindex_yPred,5])
  yesLate_postpredbycomp[4,i] = sd(rawsummaryyesLate[tempindex_yPred,5])
}

for(i in 1:6){
  tempindex = which(XnoLate[,i+1]==1)
  tempindex_yPred = tempindex+nBeta
  noLate_postpredbycomp[1,i] = mean(YnoLate[tempindex])
  noLate_postpredbycomp[2,i] = sd(YnoLate[tempindex])
  noLate_postpredbycomp[3,i] = mean(rawsummarynoLate[tempindex_yPred,5])
  noLate_postpredbycomp[4,i] = sd(rawsummarynoLate[tempindex_yPred,5])
}

ConfInt <- function(Dev,n){
  ConfInt = 1.96*(Dev/sqrt(n))
  return(ConfInt)
}

par(mfrow = c(1,1))

plot(yesEarly_postpredbycomp[1,1:6],yesEarly_postpredbycomp[3,1:6],xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)
points(noEarly_postpredbycomp[1,1:6],noEarly_postpredbycomp[3,1:6],col = '#DC143C')
points(yesLate_postpredbycomp[1,1:6],yesLate_postpredbycomp[3,1:6],col = '#00CC00')
points(noLate_postpredbycomp[1,1:6],noLate_postpredbycomp[3,1:6],col = '#7B68EE')
legend(.8, .2, legend=c("Yes Early", "No Early","Yes Late","No Late"),
       col=c("black", "red",'#00CC00','#7B68EE'), pch=c(1,1,1,1), cex=0.8)

# Test model predictions w/ epsilon,by comp type--------------------------------------------------
nBeta = 7

postpred_lo = nBeta+1
postpred_mid = postpred_lo+nExpResponses-1
postpred_hi = postpred_mid+nNovResponses

colors = c('black','#DC143C','#7B68EE','#00CC00','brown','purple')
par(mfrow=c(2,2))
plot(YyesEarly[which(XyesEarly[,1+1]==1)],rawsummaryyesEarly[(which(XyesEarly[,1+1]==1)+nBeta),1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'Yes Early',col = 'black')
for(i in 2:6){
  points(YyesEarly[which(XyesEarly[,i+1]==1)],rawsummaryyesEarly[(which(XyesEarly[,1+1]==1)+nBeta),1],
         col = colors[i])
}
abline(a=0,b=1)

legend(.6, .3, legend=c("Comp 1", "Comp 2","Comp 3","Comp 4","Comp 5","Comp 6"),
       col=colors, pch=rep(1,6), cex=0.8)

plot(YnoEarly[which(XnoEarly[,1+1]==1)],rawsummarynoEarly[(which(XnoEarly[,1+1]==1)+nBeta),1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'No Early',col = 'black')
for(i in 2:6){
  points(YnoEarly[which(XnoEarly[,i+1]==1)],rawsummarynoEarly[(which(XnoEarly[,1+1]==1)+nBeta),1],
         col = colors[i])
}
abline(a=0,b=1)

legend(.6, .3, legend=c("Comp 1", "Comp 2","Comp 3","Comp 4","Comp 5","Comp 6"),
       col=colors, pch=rep(1,6), cex=0.8)


plot(YyesLate[which(XyesLate[,1+1]==1)],rawsummaryyesLate[(which(XyesLate[,1+1]==1)+nBeta),1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'Yes Late',col = 'black')
for(i in 2:6){
  points(YyesLate[which(XyesLate[,i+1]==1)],rawsummaryyesLate[(which(XyesLate[,1+1]==1)+nBeta),1],
         col = colors[i])
}
abline(a=0,b=1)

legend(.6, .3, legend=c("Comp 1", "Comp 2","Comp 3","Comp 4","Comp 5","Comp 6"),
       col=colors, pch=rep(1,6), cex=0.8)

plot(YnoLate[which(XnoLate[,1+1]==1)],rawsummarynoLate[(which(XnoLate[,1+1]==1)+nBeta),1],xlim=c(0,1),ylim=c(0,1),
     xlab = 'Response',ylab='Posterior Predictive',main = 'No Late',col = 'black')
for(i in 2:6){
  points(YnoLate[which(XnoLate[,i+1]==1)],rawsummarynoLate[(which(XnoLate[,1+1]==1)+nBeta),1],
         col = colors[i])
}
abline(a=0,b=1)

legend(.6, .3, legend=c("Comp 1", "Comp 2","Comp 3","Comp 4","Comp 5","Comp 6"),
       col=colors, pch=rep(1,6), cex=0.8)

# Posterior predictive accuracy by comp type ------------------------------
yesEarly_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
noEarly_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
yesLate_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)
noLate_bycomp_expert <- matrix(NaN,nrow=5,ncol=6)

for(i in 1:6){
  tempindex <- which(XyesEarly[1:nExpResponses,(i+1)]==1)
  yesEarly_bycomp_expert[1,i] = mean(YyesEarly[tempindex])
  yesEarly_bycomp_expert[2,i] = sd(YyesEarly[tempindex])
  yesEarly_bycomp_expert[3,i] = mean(YyesEarly[(nExpResponses+1):testTotal])
  yesEarly_bycomp_expert[4,i] = sd(YyesEarly[(nExpResponses+1):testTotal])
  tempindex_novice <- SONAData[which(SONAData$CompType==i),]
  yesEarly_bycomp_expert[5,i] = 
    mean(tempindex_novice$EarlyResponse[which(tempindex_novice$SpeakerEarly==1)])
}

for(i in 1:6){
  tempindex = which(XnoEarly[1:nExpResponses,i+1]==1)
  noEarly_bycomp_expert[1,i] = mean(YnoEarly[tempindex])
  noEarly_bycomp_expert[2,i] = sd(YnoEarly[tempindex])
  noEarly_bycomp_expert[3,i] = mean(YnoEarly[(nExpResponses+1):testTotal])
  noEarly_bycomp_expert[4,i] = sd(YnoEarly[(nExpResponses+1):testTotal])
  tempindex_novice <- SONAData[which(SONAData$CompType==i),]
  noEarly_bycomp_expert[5,i] = 
    mean(tempindex_novice$EarlyResponse[which(tempindex_novice$SpeakerEarly==0)])
}

for(i in 1:6){
  tempindex = which(XyesLate[1:nExpResponses,i+1]==1)
  yesLate_bycomp_expert[1,i] = mean(YyesLate[tempindex])
  yesLate_bycomp_expert[2,i] = sd(YyesLate[tempindex])
  yesLate_bycomp_expert[3,i] = mean(YyesLate[(nExpResponses+1):testTotal])
  yesLate_bycomp_expert[4,i] = sd(YyesLate[(nExpResponses+1):testTotal])
  tempindex_novice <- SONAData[which(SONAData$CompType==i),]
  yesLate_bycomp_expert[5,i] = 
    mean(tempindex_novice$LateResponse[which(tempindex_novice$SpeakerLate==1)])
}

for(i in 1:6){
  tempindex = which(XnoLate[1:nExpResponses,i+1]==1)
  noLate_bycomp_expert[1,i] = mean(YnoLate[tempindex])
  noLate_bycomp_expert[2,i] = sd(YnoLate[tempindex])
  noLate_bycomp_expert[3,i] = mean(YnoLate[(nExpResponses+1):testTotal])
  noLate_bycomp_expert[4,i] = sd(YnoLate[(nExpResponses+1):testTotal])
  tempindex_novice <- SONAData[which(SONAData$CompType==i),]
  noLate_bycomp_expert[5,i] = 
    mean(tempindex_novice$LateResponse[which(tempindex_novice$SpeakerLate==0)])
}

ConfInt <- function(Dev,n){
  ConfInt = 1.96*(Dev/sqrt(n))
  return(ConfInt)
}

# Plot beta inferences from raw w/ actual data -------------------------------------------
source("PlotFunc.R")
par(mfrow=c(1,1))
par(mar = c(3.5, 3, 1, 2))
plot1 = MainPlot(rawPlotValues_EarlyYes[1,1:3],yesEarly_bycomp_expert[3,1],yesEarly_bycomp_expert[5,1:6],
                 rawPlotValues_EarlyYes[2:7,1:3],yesEarly_bycomp_expert[1,1:6],'',
                 if_legend = 1,leftmost = 1)
par(mar = c(3.5, 3, 1, 2))
plot2 = MainPlot(rawPlotValues_LateYes[1,1:3],yesLate_bycomp_expert[3,1],yesLate_bycomp_expert[5,1:6],
                 rawPlotValues_LateYes[2:7,1:3],yesLate_bycomp_expert[1,1:6],'',
                 if_legend = 0, leftmost = 0)

MainPlot(rawPlotValues_EarlyNo[1,1:3],noEarly_bycomp_expert[3,1],noEarly_bycomp_expert[5,1:6],
         rawPlotValues_EarlyNo[2:7,1:3],noEarly_bycomp_expert[1,1:6],'Interpretations of Early Game Generalizations (No)',
         if_legend = 1, leftmost = 1)

MainPlot(rawPlotValues_LateNo[1,1:3],noLate_bycomp_expert[3,1],noLate_bycomp_expert[5,1:6],
         rawPlotValues_LateNo[2:7,1:3],noLate_bycomp_expert[1,1:6],'Interpretations of Late Game Generalizations (No)',
         if_legend = 0, leftmost = 0)

# Handmade directional BF tests -------------------------------------------------------
BF_betacomparisons_yesEarly <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(rawBetayesEarlysamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(rawBetayesEarlysamples),ncol = 3)

post[,1] <- rawBetayesEarlysamples[,1]

for(i in 1:6){
  post[,2] <- rawBetayesEarlysamples[,i+1]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]>post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_betacomparisons_yesEarly[i] <- post_odds/alternate_odds
}

BF_betacomparisons_noEarly <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(rawBetanoEarlysamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(rawBetanoEarlysamples),ncol = 3)

post[,1] <- rawBetanoEarlysamples[,1]

for(i in 1:6){
  post[,2] <- rawBetanoEarlysamples[,i+1]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]<post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob) 
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_betacomparisons_noEarly[i] <- post_odds/alternate_odds
}

BF_betacomparisons_yesLate_higherExp <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(rawBetayesLatesamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(rawBetayesLatesamples),ncol = 3)

post[,1] <- rawBetayesLatesamples[,1]

for(i in 1:6){
  post[,2] <- rawBetayesLatesamples[,i+1]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]<post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_betacomparisons_yesLate_higherExp[i] <- post_odds/alternate_odds
}

BF_betacomparisons_yesLate <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(rawBetayesLatesamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(rawBetayesLatesamples),ncol = 3)

post[,1] <- rawBetayesLatesamples[,1]

for(i in 1:6){
  post[,2] <- rawBetayesLatesamples[,i+1]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]>post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_betacomparisons_yesLate[i] <- post_odds/alternate_odds
}

BF_betacomparisons_noLate <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(rawBetanoLatesamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(rawBetanoLatesamples),ncol = 3)

post[,1] <- rawBetanoLatesamples[,1]

for(i in 1:6){
  post[,2] <- rawBetanoLatesamples[,i+1]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]<post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_betacomparisons_noLate[i] <- post_odds/alternate_odds
}

par(mfrow = c(3,2))

for(i in 1:6){
  title = paste0("Yes Early: B0 vs. B",i,sep = NULL)
  plot(density(rawBetayesEarlysamples[,1]),main = title,xlim=c(.45,.75))
  lines(density(rawBetayesEarlysamples[,i+1]),col='#DC143C')
}

for(i in 1:6){
  title = paste0("No Early: B0 vs. B",i,sep = NULL)
  plot(density(rawBetanoEarlysamples[,1]),main = title,xlim=c(.27,.57))
  lines(density(rawBetanoEarlysamples[,i+1]),col='#DC143C')
}

for(i in 1:6){
  title = paste0("Yes Late: B0 vs. B",i,sep = NULL)
  plot(density(rawBetayesLatesamples[,1]),main = title,xlim=c(.45,.75))
  lines(density(rawBetayesLatesamples[,i+1]),col='#DC143C')
}

for(i in 1:6){
  title = paste0("No Late: B0 vs. B",i,sep = NULL)
  plot(density(rawBetanoLatesamples[,1]),main = title,xlim=c(.27,.57))
  lines(density(rawBetanoLatesamples[,i+1]),col='#DC143C')
}

# Handmade directional BF tests, 12 Beta -------------------------------------------------------
BF_12B_betacomparisons_yesEarly <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(overfullBetayesEarlysamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(overfullBetayesEarlysamples),ncol = 3)



for(i in 1:6){
  post[,1] <- overfullBetayesEarlysamples[,i+6]
  post[,2] <- overfullBetayesEarlysamples[,i]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]>post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_12B_betacomparisons_yesEarly[i] <- post_odds/alternate_odds
}

BF_12B_betacomparisons_noEarly <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(overfullBetanoEarlysamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(overfullBetanoEarlysamples),ncol = 3)

for(i in 1:6){
  post[,1] <- overfullBetanoEarlysamples[,i+6]
  post[,2] <- overfullBetanoEarlysamples[,i]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]<post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_12B_betacomparisons_noEarly[i] <- post_odds/alternate_odds
}

BF_12B_betacomparisons_yesLate_higherExp <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(overfullBetayesLatesamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(overfullBetayesLatesamples),ncol = 3)

for(i in 1:6){
  post[,1] <- overfullBetayesLatesamples[,i+6]
  post[,2] <- overfullBetayesLatesamples[,i]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]<post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_12B_betacomparisons_yesLate_higherExp[i] <- post_odds/alternate_odds
}

BF_12B_betacomparisons_yesLate <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(overfullBetayesLatesamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(overfullBetayesLatesamples),ncol = 3)

for(i in 1:6){
  post[,1] <- overfullBetayesLatesamples[,i+6]
  post[,2] <- overfullBetayesLatesamples[,i]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]>post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_12B_betacomparisons_yesLate[i] <- post_odds/alternate_odds
}

BF_12B_betacomparisons_noLate <- rep(NaN,6)

alternate <- matrix(NaN,nrow = nrow(overfullBetanoLatesamples),ncol = 3)
post <- matrix(NaN,nrow = nrow(overfullBetanoLatesamples),ncol = 3)

for(i in 1:6){
  post[,1] <- overfullBetanoLatesamples[,i+6]
  post[,2] <- overfullBetanoLatesamples[,i]
  for(j in 1:nrow(post)){
    post[j,3] <- post[j,1]<post[j,2]
  }
  post_prob <- mean(post[,3])
  if(post_prob==1){
    post_prob = (nSamples-1)/nSamples
  }
  if(post_prob==0){
    post_prob = 1/nSamples
  }
  post_odds <- post_prob/(1-post_prob)
  
  alternate_prob <- 1 - post_prob
  alternate_odds <- alternate_prob/(1-alternate_prob)
  
  BF_12B_betacomparisons_noLate[i] <- post_odds/alternate_odds
}

for(i in 1:6){
  title = paste0("12 Beta Yes Early: B",i, "vs. B",i+6,sep = NULL)
  plot(density(overfullBetayesEarlysamples[,1]),main = title,xlim=c(.45,.75),col = '#DC143C')
  lines(density(overfullBetayesEarlysamples[,i+1]))
}

for(i in 1:6){
  title = paste0("12 Beta No Early: B0 vs. B",i,sep = NULL)
  plot(density(overfullBetanoEarlysamples[,1]),main = title,xlim=c(.27,.57),col='#DC143C')
  lines(density(overfullBetanoEarlysamples[,i+1]))
}

for(i in 1:6){
  title = paste0("12 Beta Yes Late: B0 vs. B",i,sep = NULL)
  plot(density(overfullBetayesLatesamples[,1]),main = title,xlim=c(.45,.75),col='#DC143C')
  lines(density(overfullBetayesLatesamples[,i+1]))
}

for(i in 1:6){
  title = paste0("12 Beta No Late: B0 vs. B",i,sep = NULL)
  plot(density(overfullBetanoLatesamples[,1]),main = title,xlim=c(.27,.57),col='#DC143C')
  lines(density(overfullBetanoLatesamples[,i+1]))
}

BF_betacomparisons_results <- data.frame(matrix(NaN,nrow = 6,ncol = 10))
names(BF_betacomparisons_results) <- c("Yes Early","12 Beta Yes Early","No Early","12 Beta No Early",
                                       "Yes Late","12 Beta Yes Late","Yes Late Higher?","12 Beta Yes Late Higher?",
                                       "No Late","12 Beta No Late")
BF_betacomparisons_results[,1] <- round(BF_betacomparisons_yesEarly,3)
BF_betacomparisons_results[,2] <- round(BF_12B_betacomparisons_yesEarly,3)
BF_betacomparisons_results[,3] <- round(BF_betacomparisons_noEarly,3)
BF_betacomparisons_results[,4] <- round(BF_12B_betacomparisons_noEarly,3)
BF_betacomparisons_results[,5] <- round(BF_betacomparisons_yesLate,3)
BF_betacomparisons_results[,6] <- round(BF_12B_betacomparisons_yesLate,3)
BF_betacomparisons_results[,7] <- round(BF_12B_betacomparisons_yesLate_higherExp,3)
BF_betacomparisons_results[,8] <- round(BF_12B_betacomparisons_yesLate_higherExp,3)
BF_betacomparisons_results[,9] <- round(BF_betacomparisons_noLate,3)
BF_betacomparisons_results[,10] <- round(BF_12B_betacomparisons_noLate,3)

# Correlation matrix of beta pulls ----------------------------------------

cordata_rawBetayesEarlysamples <- data.frame(rawBetayesEarlysamples)
cordata_rawBetanoEarlysamples <- data.frame(rawBetanoEarlysamples)
cordata_rawBetayesLatesamples <- data.frame(rawBetayesLatesamples)
cordata_rawBetanoLatesamples <- data.frame(rawBetanoLatesamples)
names(cordata_rawBetayesEarlysamples) <- c('B0','B1','B2','B3','B4','B5','B6')
names(cordata_rawBetanoEarlysamples) <- c('B0','B1','B2','B3','B4','B5','B6')
names(cordata_rawBetayesLatesamples) <- c('B0','B1','B2','B3','B4','B5','B6')
names(cordata_rawBetanoLatesamples) <- c('B0','B1','B2','B3','B4','B5','B6')

cormat_rawBetayesEarlysamples <- cor(rawBetayesEarlysamples)
cormat_rawBetanoEarlysamples <- cor(rawBetanoEarlysamples)
cormat_rawBetayesLatesamples <- cor(rawBetayesLatesamples)
cormat_rawBetanoLatesamples <- cor(rawBetanoLatesamples)

par(mfrow=c(2,2))
corrplot(cormat_rawBetayesEarlysamples, method="number",type = 'upper',title = 'Early Yes')
corrplot(cormat_rawBetanoEarlysamples, method="number",type = 'upper', title = 'Early No')
corrplot(cormat_rawBetayesLatesamples, method="number",type = 'upper',title = 'Late Yes')
corrplot(cormat_rawBetanoLatesamples, method="number",type = 'upper', title = 'Late No')

# Hypothesis test using Savage-Dickey -------------------------------------

prior_density = dnorm(0,0,.3537)
posterior_density_yesEarly <- rep(NaN,6)
posterior_density_noEarly <- rep(NaN,6)
posterior_density_yesLate <- rep(NaN,6)
posterior_density_noLate <- rep(NaN,6)

posterior_12B_density_yesEarly <- rep(NaN,6)
posterior_12B_density_noEarly <- rep(NaN,6)
posterior_12B_density_yesLate <- rep(NaN,6)
posterior_12B_density_noLate <- rep(NaN,6)

for(i in 1:6){
  tempdif <- rawBetayesEarlysamples[,1]-rawBetayesEarlysamples[,i+1]
  posterior_density_yesEarly[i] <- dlogspline(0,logspline(tempdif))
  tempdif <- overfullBetayesEarlysamples[,i+6]-overfullBetayesEarlysamples[,i]
  posterior_12B_density_yesEarly[i] <- dlogspline(0,logspline(tempdif))
  
  tempdif <- rawBetanoEarlysamples[,1]-rawBetanoEarlysamples[,i+1]
  posterior_density_noEarly[i] <- dlogspline(0,logspline(tempdif))
  tempdif <- overfullBetanoEarlysamples[,i+6]-overfullBetanoEarlysamples[,i]
  posterior_12B_density_noEarly[i] <- dlogspline(0,logspline(tempdif))
  
  tempdif <- rawBetayesLatesamples[,1]-rawBetayesLatesamples[,i+1]
  posterior_density_yesLate[i] <- dlogspline(0,logspline(tempdif))
  tempdif <- overfullBetayesLatesamples[,i+6]-overfullBetayesLatesamples[,i]
  posterior_12B_density_yesLate[i] <- dlogspline(0,logspline(tempdif))
  
  tempdif <- rawBetanoLatesamples[,1]-rawBetanoLatesamples[,i+1]
  posterior_density_noLate[i] <- dlogspline(0,logspline(tempdif))
  tempdif <- overfullBetanoLatesamples[,i+6]-overfullBetanoLatesamples[,i]
  posterior_12B_density_noLate[i] <- dlogspline(0,logspline(tempdif))
}

SavDick <- data.frame(matrix(NaN,nrow=6,ncol=8))
names(SavDick) <- c("Yes Early","12B Yes Early","No Early","12B No Early","Yes Late","12B Yes Late",
                    "No Late","12B No Late")

SavDick[,1] <- posterior_density_yesEarly/prior_density
SavDick[,2] <- posterior_12B_density_yesEarly/prior_density
SavDick[,3] <- posterior_density_noEarly/prior_density
SavDick[,4] <- posterior_12B_density_noEarly/prior_density
SavDick[,5] <- posterior_density_yesLate/prior_density
SavDick[,6] <- posterior_12B_density_yesLate/prior_density
SavDick[,7] <- posterior_density_noLate/prior_density
SavDick[,8] <- posterior_12B_density_noLate/prior_density

# Plot beta inferences from overfull -------------------------------------------
par(mfrow = c(2,2))

overfullPlotValues_EarlyYes_Expert <- matrix(NaN, nrow = 6,ncol = 3)
overfullPlotValues_EarlyNo_Expert <- matrix(NaN, nrow = 6,ncol = 3)
overfullPlotValues_LateYes_Expert <- matrix(NaN, nrow = 6,ncol = 3)
overfullPlotValues_LateNo_Expert <- matrix(NaN, nrow = 6,ncol = 3)

overfullPlotValues_EarlyYes_Novice <- matrix(NaN, nrow = 6,ncol = 3)
overfullPlotValues_EarlyNo_Novice <- matrix(NaN, nrow = 6,ncol = 3)
overfullPlotValues_LateYes_Novice <- matrix(NaN, nrow = 6,ncol = 3)
overfullPlotValues_LateNo_Novice <- matrix(NaN, nrow = 6,ncol = 3)

loconf = nrow(overfullBetayesEarlysamples) * .025
hiconf = nrow(overfullBetayesEarlysamples) * .975


for(i in 1:6){
  
  tempsort <- sort(overfullBetayesEarlysamples[,i])
  overfullPlotValues_EarlyYes_Expert[i,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_EarlyYes_Expert[i,2] <- tempsort[loconf]
  overfullPlotValues_EarlyYes_Expert[i,3] <- tempsort[hiconf]
  
  tempsort <- sort(overfullBetanoEarlysamples[,i])
  overfullPlotValues_EarlyNo_Expert[i,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_EarlyNo_Expert[i,2] <- tempsort[loconf]
  overfullPlotValues_EarlyNo_Expert[i,3] <- tempsort[hiconf]
  
  tempsort <- sort(overfullBetayesLatesamples[,i])
  overfullPlotValues_LateYes_Expert[i,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_LateYes_Expert[i,2] <- tempsort[loconf]
  overfullPlotValues_LateYes_Expert[i,3] <- tempsort[hiconf]
  
  tempsort <- sort(overfullBetanoLatesamples[,i])
  overfullPlotValues_LateNo_Expert[i,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_LateNo_Expert[i,2] <- tempsort[loconf]
  overfullPlotValues_LateNo_Expert[i,3] <- tempsort[hiconf]
}

for(i in 7:12){
  
  tempsort <- sort(overfullBetayesEarlysamples[,i])
  overfullPlotValues_EarlyYes_Novice[i-6,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_EarlyYes_Novice[i-6,2] <- tempsort[loconf]
  overfullPlotValues_EarlyYes_Novice[i-6,3] <- tempsort[hiconf]
  
  tempsort <- sort(overfullBetanoEarlysamples[,i])
  overfullPlotValues_EarlyNo_Novice[i-6,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_EarlyNo_Novice[i-6,2] <- tempsort[loconf]
  overfullPlotValues_EarlyNo_Novice[i-6,3] <- tempsort[hiconf]
  
  tempsort <- sort(overfullBetayesLatesamples[,i])
  overfullPlotValues_LateYes_Novice[i-6,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_LateYes_Novice[i-6,2] <- tempsort[loconf]
  overfullPlotValues_LateYes_Novice[i-6,3] <- tempsort[hiconf]
  
  tempsort <- sort(overfullBetanoLatesamples[,i])
  overfullPlotValues_LateNo_Novice[i-6,1] <- mean(tempsort[1500:1001])
  overfullPlotValues_LateNo_Novice[i-6,2] <- tempsort[loconf]
  overfullPlotValues_LateNo_Novice[i-6,3] <- tempsort[hiconf]
}

plot(1:6,overfullPlotValues_EarlyYes_Expert[1:6,1],ylim = c(0,1),main = "Yes Early (All 12 Betas)",xlab = 'beta',ylab = 'Interpretation')

for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_EarlyYes_Expert[i,2],
         x1 = i, y1 = overfullPlotValues_EarlyYes_Expert[i,3],code = 3,
         col = 'black',angle = 90, length = .1)
}

points(1:6,overfullPlotValues_EarlyYes_Novice[1:6,1],col = '#00CC00',pch = 2)

abline(h= yesEarly_bycomp_expert[3,1],col='#00CC00')

for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_EarlyYes_Novice[i,2],
         x1 = i, y1 = overfullPlotValues_EarlyYes_Novice[i,3],code = 3,
         col = '#00CC00',angle = 90, length = .1,lty = 3)
}

legend(3.5,.35,legend = c("Novice","Expert","Novice Mean (Empirical)"),pch = c(2,1,NA),lty = c(NA,NA,1),
       col = c('#00CC00','black'),cex = .8)

plot(1:6,overfullPlotValues_LateYes_Expert[1:6,1],ylim = c(0,1),main = "Yes Late (All 12 Betas)",xlab = 'beta',ylab = 'Interpretation')
for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_LateYes_Expert[i,2],
         x1 = i, y1 = overfullPlotValues_LateYes_Expert[i,3],code = 3,
         col = 'black',angle = 90, length = .1)
}

points(1:6,overfullPlotValues_LateYes_Novice[1:6,1],col = '#00CC00',pch = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_LateYes_Novice[i,2],
         x1 = i, y1 = overfullPlotValues_LateYes_Novice[i,3],code = 3,
         col = '#00CC00',angle = 90, length = .1,lty = 3)
}

abline(h= yesLate_bycomp_expert[3,1],col='#00CC00')

plot(1:6,overfullPlotValues_EarlyNo_Expert[1:6,1],ylim = c(0,1),main = "No Early (All 12 Betas)",xlab = 'beta',ylab = 'Interpretation')
for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_EarlyNo_Expert[i,2],
         x1 = i, y1 = overfullPlotValues_EarlyNo_Expert[i,3],code = 3,
         col = 'black',angle = 90, length = .1)
}

points(1:6,overfullPlotValues_EarlyNo_Novice[1:6,1],col = '#00CC00',pch = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_EarlyNo_Novice[i,2],
         x1 = i, y1 = overfullPlotValues_EarlyNo_Novice[i,3],code = 3,
         col = '#00CC00',angle = 90, length = .1,lty = 3)
}

abline(h= noEarly_bycomp_expert[3,1],col='#00CC00')

plot(1:6,overfullPlotValues_LateNo_Expert[1:6,1],ylim = c(0,1),main = "No Late (All 12 Betas)",xlab = 'beta',ylab = 'Interpretation')
for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_LateNo_Expert[i,2],
         x1 = i, y1 = overfullPlotValues_LateNo_Expert[i,3],code = 3,
         col = 'black',angle = 90, length = .1)
}

points(1:6,overfullPlotValues_LateNo_Novice[1:6,1],col = '#00CC00',pch = 2)

for(i in 1:6){
  arrows(x0 = i, y0 = overfullPlotValues_LateNo_Novice[i,2],
         x1 = i, y1 = overfullPlotValues_LateNo_Novice[i,3],code = 3,
         col = '#00CC00',angle = 90, length = .1,lty = 3)
  
  abline(h= noLate_bycomp_expert[3,1],col='#00CC00')
}


# Graph alphas by skill level ---------------------------------------------

Skill_options <- c("Silver","Gold","Platinum","Diamond","Master","Challenger")
ChampFam_total_options <- c("0-30","31-60","61-90","91-120","121-150")
for(i in 1:nrow(OnlineData)){
  OnlineData$Skill_Index[i] <- which(Skill_options==OnlineData$Skill[i])
  OnlineData$ChampFam_total_index[i] <- which(ChampFam_total_options==OnlineData$ChampFam_total[i])
}

Skill_Alpha_yesEarly <- matrix(NaN,nrow = nExperts,ncol = 4)
nYpred = nExperts*nTrials+nNovices*nTrials

for(i in 1:nExperts){
  ticker = i*nTrials*2
  Skill_Alpha_yesEarly[i,1] <- OnlineData$Participant[ticker]
  Skill_Alpha_yesEarly[i,2] <- rawsummaryyesEarly[i+nYpred+7,1]
  Skill_Alpha_yesEarly[i,3] <- OnlineData$Skill_Index[ticker]
  Skill_Alpha_yesEarly[i,4] <- OnlineData$ChampFam_total_index[ticker]
}


par(mfrow = c(2,2))
plot(Skill_Alpha_yesEarly[,3],Skill_Alpha_yesEarly[,2])
plot(Skill_Alpha_yesEarly[,4],Skill_Alpha_yesEarly[,2])

Skill_Alpha_noEarly <- matrix(NaN,nrow = nExperts,ncol = 4)

for(i in 1:nExperts){
  ticker = i*nTrials*2
  Skill_Alpha_noEarly[i,1] <- OnlineData$Participant[ticker]
  Skill_Alpha_noEarly[i,2] <- rawsummarynoEarly[i+nYpred+7,1]
  Skill_Alpha_noEarly[i,3] <- OnlineData$Skill_Index[ticker]
  Skill_Alpha_noEarly[i,4] <- OnlineData$ChampFam_total_index[ticker]
}

plot(Skill_Alpha_noEarly[,3],Skill_Alpha_noEarly[,2])
plot(Skill_Alpha_noEarly[,4],Skill_Alpha_noEarly[,2])

Skill_Alpha_yesLate <- matrix(NaN,nrow = nExperts,ncol = 4)

for(i in 1:nExperts){
  ticker = i*nTrials*2
  Skill_Alpha_yesLate[i,1] <- OnlineData$Participant[ticker]
  Skill_Alpha_yesLate[i,2] <- rawsummaryyesLate[i+nYpred+7,1]
  Skill_Alpha_yesLate[i,3] <- OnlineData$Skill_Index[ticker]
  Skill_Alpha_yesLate[i,4] <- OnlineData$ChampFam_total_index[ticker]
}

plot(Skill_Alpha_yesLate[,3],Skill_Alpha_yesLate[,2])
plot(Skill_Alpha_yesLate[,4],Skill_Alpha_yesLate[,2])

Skill_Alpha_noLate <- matrix(NaN,nrow = nExperts,ncol = 4)

for(i in 1:nExperts){
  ticker = i*nTrials*2
  Skill_Alpha_noLate[i,1] <- OnlineData$Participant[ticker]
  Skill_Alpha_noLate[i,2] <- rawsummarynoLate[i+nYpred+7,1]
  Skill_Alpha_noLate[i,3] <- OnlineData$Skill_Index[ticker]
  Skill_Alpha_noLate[i,4] <- OnlineData$ChampFam_total_index[ticker]
}

plot(Skill_Alpha_noLate[,3],Skill_Alpha_noLate[,2])
plot(Skill_Alpha_noLate[,4],Skill_Alpha_noLate[,2])


# Calculate and graph epsilon inferences ----------------------------------
nYpred = nExperts*nTrials+nNovices*nTrials
Epsilons_yesEarly <- matrix(NaN,nrow = length(YyesEarly),ncol = 2)
ticker = 0
for(i in 1:nExperts){
  for(j in 1:nTrials){
    ticker = ((i-1)*nTrials)+ j
    Beta_Index = which(XyesEarly[ticker,]==1)
    B_N = rawsummaryyesEarly[Beta_Index,1]
    Alpha_index = i+7+nYpred
    tempAlpha = rawsummaryyesEarly[Alpha_index,1]
    tempYpred = rawsummaryyesEarly[7+ticker,1]
    Epsilons_yesEarly[ticker,1] = tempYpred-B_N-tempAlpha
    Epsilons_yesEarly[ticker,2] = tempYpred
  }
}

B_0 = rawsummaryyesEarly[1,1]
for(i in 1:nNovices){
  for(j in 1:nTrials){
    ticker2 = ((i-1)*nTrials)+j + ticker
    Alpha_index = i+7+nYpred+nExperts
    tempAlpha = rawsummaryyesEarly[Alpha_index,1]
    tempYpred = rawsummaryyesEarly[7+ticker2,1]
    Epsilons_yesEarly[ticker2,1] = tempYpred-B_0-tempAlpha
    Epsilons_yesEarly[ticker2,2] = tempYpred
  }
}
par(mfrow = c(2,2))
plot(Epsilons_yesEarly[1:nExpResponses,2],Epsilons_yesEarly[1:nExpResponses,1],col = 'blue',xlim = c(0,1),ylim = c(-.01,.01))
points(Epsilons_yesEarly[(nExpResponses+1):length(YyesEarly),2],Epsilons_yesEarly[(nExpResponses+1):length(YyesEarly),1],
       col = 'red')
abline(h=0)

EpsilonGraph_yesEarly <- rep(NaN,nNovices+nExperts)
EpsilonGraph_yesLate <- rep(NaN,nNovices+nExperts)
for(i in 1:nExperts){
  eTickerhi = i*6
  eTickerlo = (i*6)-5
  EpsilonGraph_yesEarly[i] <- mean(Epsilons_yesEarly[eTickerlo:eTickerhi,1])
  EpsilonGraph_yesLate[i] <- mean(Epsilons_yesLate[eTickerlo:eTickerhi,1])
}

for(i in (nExperts+1):(nNovices+nExperts)){
  eTickerhi = i*6
  eTickerlo = (i*6)-5
  EpsilonGraph_yesEarly[i] <- mean(Epsilons_yesEarly[eTickerlo:eTickerhi,1])
  EpsilonGraph_yesLate[i] <- mean(Epsilons_yesLate[eTickerlo:eTickerhi,1])
}

par(mfrow = c(2,2))
hist(EpsilonGraph_yesEarly[1:nExperts])
hist(EpsilonGraph_yesEarly[(nExperts+1):(nNovices+nExperts)])

ExpertEpsilonMean_yesEarly = mean(EpsilonGraph_yesEarly[1:nExperts])
ExpertEpsilonSD_yesEarly = sd(EpsilonGraph_yesEarly[1:nExperts])
NoviceEpsilonMean_yesEarly = 
  mean(EpsilonGraph_yesEarly[(nExperts+1):(nNovices+nExperts)])
NoviceEpsilonSD_yesEarly = 
  sd(EpsilonGraph_yesEarly[(nExperts+1):(nNovices+nExperts)])

hist(EpsilonGraph_yesLate[1:nExperts])
hist(EpsilonGraph_yesLate[(nExperts+1):(nNovices+nExperts)])

ExpertEpsilonMean_yesLate = mean(EpsilonGraph_yesLate[1:nExperts])
ExpertEpsilonSD_yesLate = sd(EpsilonGraph_yesLate[1:nExperts])
NoviceEpsilonMean_yesLate = 
  mean(EpsilonGraph_yesLate[(nExperts+1):(nNovices+nExperts)])
NoviceEpsilonSD_yesLate = 
  sd(EpsilonGraph_yesLate[(nExperts+1):(nNovices+nExperts)])

ExpertEpsilonMean_all = mean(c(EpsilonGraph_yesEarly[1:nExperts],
                               EpsilonGraph_yesLate[1:nExperts]))
ExpertEpsilonSD_all = sd(c(EpsilonGraph_yesEarly[1:nExperts],
                           EpsilonGraph_yesLate[1:nExperts]))
NoviceEpsilonMean_all = 
  mean(c(EpsilonGraph_yesEarly[(nExperts+1):(nNovices+nExperts)],
         EpsilonGraph_yesLate[(nExperts+1):(nNovices+nExperts)]))
NoviceEpsilonSD_all = 
  sd(c(EpsilonGraph_yesEarly[(nExperts+1):(nNovices+nExperts)],
       EpsilonGraph_yesLate[(nExperts+1):(nNovices+nExperts)]))

Epsilons_noEarly <- matrix(NaN,nrow = length(YnoEarly),ncol = 2)
ticker = 0
for(i in 1:nExperts){
  for(j in 1:nTrials){
    ticker = ((i-1)*nTrials)+ j
    Beta_Index = which(XnoEarly[ticker,]==1)
    B_N = rawsummarynoEarly[Beta_Index,1]
    Alpha_index = i+7+nYpred
    tempAlpha = rawsummarynoEarly[Alpha_index,1]
    tempYpred = rawsummarynoEarly[7+ticker,1]
    Epsilons_noEarly[ticker,1] = tempYpred-B_N-tempAlpha
    Epsilons_noEarly[ticker,2] = tempYpred
  }
}

B_0 = rawsummarynoEarly[1,1]
for(i in 1:nNovices){
  for(j in 1:nTrials){
    ticker2 = ((i-1)*nTrials)+j + ticker
    Alpha_index = i+7+nYpred+nExperts
    tempAlpha = rawsummarynoEarly[Alpha_index,1]
    tempYpred = rawsummarynoEarly[7+ticker2,1]
    Epsilons_noEarly[ticker2,1] = tempYpred-B_0-tempAlpha
    Epsilons_noEarly[ticker2,2] = tempYpred
  }
}

plot(Epsilons_noEarly[1:nExpResponses,2],Epsilons_noEarly[1:nExpResponses,1],col = 'blue',xlim = c(0,1),ylim = c(-.01,.01))
points(Epsilons_noEarly[(nExpResponses+1):length(YnoEarly),2],Epsilons_noEarly[(nExpResponses+1):length(YnoEarly),1],
       col = 'red')
abline(h=0)

Epsilons_yesLate <- matrix(NaN,nrow = length(YyesLate),ncol = 2)
ticker = 0
for(i in 1:nExperts){
  for(j in 1:nTrials){
    ticker = ((i-1)*nTrials)+ j
    Beta_Index = which(XyesLate[ticker,]==1)
    B_N = rawsummaryyesLate[Beta_Index,1]
    Alpha_index = i+7+nYpred
    tempAlpha = rawsummaryyesLate[Alpha_index,1]
    tempYpred = rawsummaryyesLate[7+ticker,1]
    Epsilons_yesLate[ticker,1] = tempYpred-B_N-tempAlpha
    Epsilons_yesLate[ticker,2] = tempYpred
  }
}

B_0 = rawsummaryyesLate[1,1]
for(i in 1:nNovices){
  for(j in 1:nTrials){
    ticker2 = ((i-1)*nTrials)+j + ticker
    Alpha_index = i+7+nYpred+nExperts
    tempAlpha = rawsummaryyesLate[Alpha_index,1]
    tempYpred = rawsummaryyesLate[7+ticker2,1]
    Epsilons_yesLate[ticker2,1] = tempYpred-B_0-tempAlpha
    Epsilons_yesLate[ticker2,2] = tempYpred
  }
}

plot(Epsilons_yesLate[1:nExpResponses,2],Epsilons_yesLate[1:nExpResponses,1],col = 'blue',xlim = c(0,1),ylim = c(-.01,.01))
points(Epsilons_yesLate[(nExpResponses+1):length(YyesLate),2],Epsilons_yesLate[(nExpResponses+1):length(YyesLate),1],
       col = 'red')
abline(h=0)

hist(Epsilons_yesLate[1:nExpResponses,1])
hist(Epsilons_yesLate[(nExpResponses+1):length(YyesLate),1])

Epsilons_noLate <- matrix(NaN,nrow = length(YnoLate),ncol = 2)
ticker = 0
for(i in 1:nExperts){
  for(j in 1:nTrials){
    ticker = ((i-1)*nTrials)+ j
    Beta_Index = which(XnoLate[ticker,]==1)
    B_N = rawsummarynoLate[Beta_Index,1]
    Alpha_index = i+7+nYpred
    tempAlpha = rawsummarynoLate[Alpha_index,1]
    tempYpred = rawsummarynoLate[7+ticker,1]
    Epsilons_noLate[ticker,1] = tempYpred-B_N-tempAlpha
    Epsilons_noLate[ticker,2] = tempYpred
  }
}

B_0 = rawsummarynoLate[1,1]
for(i in 1:nNovices){
  for(j in 1:nTrials){
    ticker2 = ((i-1)*nTrials)+j + ticker
    Alpha_index = i+7+nYpred+nExperts
    tempAlpha = rawsummarynoLate[Alpha_index,1]
    tempYpred = rawsummarynoLate[7+ticker2,1]
    Epsilons_noLate[ticker2,1] = tempYpred-B_0-tempAlpha
    Epsilons_noLate[ticker2,2] = tempYpred
  }
}

plot(Epsilons_noLate[1:nExpResponses,2],Epsilons_noLate[1:nExpResponses,1],col = 'blue',xlim = c(0,1),ylim = c(-.01,.01))
points(Epsilons_noLate[(nExpResponses+1):length(YnoLate),2],Epsilons_noLate[(nExpResponses+1):length(YnoLate),1],
       col = 'red')
abline(h=0)


# Expert TVJ performance by comp type -------------------------------------

TVJ_Early_bycomp <- c(mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==1)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==2)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==3)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==4)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==5)]),
                      mean(OnlineData$EarlyTVJ[which(OnlineData$CompType==6)]))

SE_Binary = function(n,p.est){
  answer = sqrt(p.est*(1-p.est)/n)
  return(answer)
}

plot(1:6,TVJ_Early_bycomp,main = "Expert Endorsement Rate by Comp Type (Early)",xlab = 'Comp Type',ylab = '% "True" Responses',
     ylim = c(0,1))
abline(h=.5,lty = 2)

for(i in 1:6){
  SE = SE_Binary(nExpResponses/3,TVJ_Early_bycomp[i])
  arrows(x0 = i, y0 = TVJ_Early_bycomp[i]+(SE*1.96),
         x1 = i, y1 = TVJ_Early_bycomp[i]-(SE*1.96),code = 3,
         angle = 90, length = .1,lty = 3)
}

TVJ_Late_bycomp <- c(mean(OnlineData$LateTVJ[which(OnlineData$CompType==1)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==2)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==3)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==4)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==5)]),
                     mean(OnlineData$LateTVJ[which(OnlineData$CompType==6)]))

plot(1:6,TVJ_Late_bycomp,main = "Expert Endorsement Rate by Comp Type (Late)",xlab = 'Comp Type',ylab = '% "True" Responses',
     ylim = c(0,1))
abline(h=.5,lty = 2)
for(i in 1:6){
  SE = SE_Binary(nExpResponses/3,TVJ_Late_bycomp[i])
  arrows(x0 = i, y0 = TVJ_Late_bycomp[i]+(SE*1.96),
         x1 = i, y1 = TVJ_Late_bycomp[i]-(SE*1.96),code = 3,
         angle = 90, length = .1,lty = 3)
}


# Correlation between expert priors and interpretation? -------------------

Prior_Int_yes <- matrix(NaN,nrow = (nrow(OnlineData)/2),ncol = 4)
Prior_Int_no <- matrix(NaN,nrow = (nrow(OnlineData)/2),ncol = 4)
Earlyticker_yes = 1
Lateticker_yes = 1
Earlyticker_no = 1
Lateticker_no = 1

for(i in 1:totalExpert){
  if(OnlineData$SpeakerEarly[i]==1){
    Prior_Int_yes[Earlyticker_yes,1] <- OnlineData$EarlyPrior[i]
    Prior_Int_yes[Earlyticker_yes,2] <- OnlineData$EarlyResponse[i]
    Earlyticker_yes = Earlyticker_yes+1
  }
  if(OnlineData$SpeakerEarly[i]==0){
    Prior_Int_no[Earlyticker_no,1] <- OnlineData$EarlyPrior[i]
    Prior_Int_no[Earlyticker_no,2] <- OnlineData$EarlyResponse[i]
    Earlyticker_no = Earlyticker_no+1
  }
  if(OnlineData$SpeakerLate[i]==1){
    Prior_Int_yes[Lateticker_yes,3] <- OnlineData$LatePrior[i]
    Prior_Int_yes[Lateticker_yes,4] <- OnlineData$LateResponse[i]
    Lateticker_yes = Lateticker_yes+1
  }
  if(OnlineData$SpeakerLate[i]==0){
    Prior_Int_no[Lateticker_no,3] <- OnlineData$LatePrior[i]
    Prior_Int_no[Lateticker_no,4] <- OnlineData$LateResponse[i]
    Lateticker_no = Lateticker_no+1
  }
}

Prior_Int_yes <- data.frame(Prior_Int_yes)
Prior_Int_no <- data.frame(Prior_Int_no)

names(Prior_Int_yes) <- c("Early Prior","Early Int","Late Prior","Late Int")
names(Prior_Int_no) <- c("Early Prior","Early Int","Late Prior","Late Int")
par(mfrow = c(1,1))
cormatPrior_Int_yes <- cor(Prior_Int_yes)
cormatPrior_Int_no <- cor(Prior_Int_no)
corrplot(cormatPrior_Int_yes, method="number",type = 'upper',title = 'Prior->Interpretation Yes')
corrplot(cormatPrior_Int_no, method="number",type = 'upper', title = 'Prior->Interpretation No')

cormatPrior_Int_yes_pval <- rcorr(as.matrix(Prior_Int_yes))
cormatPrior_Int_no_pval <- rcorr(as.matrix(Prior_Int_no))


# Make tables -------------------------------------------------------------

Table_contents_early <- data.frame(matrix(NaN,nrow= 5,ncol = 7))
names(Table_contents_early) <- c('1','2','3','4','5','6','Novice')

Table_contents_early[1,1:7] <- c('Strong Early','Middling Early','Weak Early','Strong Late','Middling Late','Weak Late','-')

Table_contents_early[2,1:7] <- round(c(yesEarly_bycomp_expert[1,1:6],yesEarly_bycomp_expert[3,1]),3)

Table_contents_early[3,1:7] <- round(c(rawPlotValues_EarlyYes[2:7,1],rawPlotValues_EarlyYes[1,1]),3)
3.293
for(i in 2:7){
  Table_contents_early[4,i-1] <- round(sd(rawBetayesEarlysamples[,i]),3)
}

Table_contents_early[4,7] <- round(sd(rawBetayesEarlysamples[,1]),3)

Table_contents_early[5,1:6] <- round(SavDick$`Yes Early`[1:6],3)

Table_contents_early[5,7] <- '-'

Table_rownames <- c('Comp. Type','Mean (Empirical)','Mean (Inferred)','SD','BF')
row.names(Table_contents_early) <- Table_rownames

kbl_Table_early <- kbl(Table_contents_early,booktabs = T,"latex")%>%
  kable_styling(latex_options = "striped")

Table_contents_late <- data.frame(matrix(NaN,nrow= 5,ncol = 7))
names(Table_contents_late) <- c('1','2','3','4','5','6','Novice')

Table_contents_late[1,1:7] <- c('Strong Early','Middling Early','Weak Early','Strong Late','Middling Late','Weak Late','-')

Table_contents_late[2,1:7] <- round(c(yesLate_bycomp_expert[1,1:6],yesLate_bycomp_expert[3,1]),3)

Table_contents_late[3,1:7] <- round(c(rawPlotValues_LateYes[2:7,1],rawPlotValues_LateYes[1,1]),3)

for(i in 2:7){
  Table_contents_late[4,i-1] <- round(sd(rawBetayesLatesamples[,i]),3)
}

Table_contents_late[4,7] <- round(sd(rawBetayesLatesamples[,1]),3)

Table_contents_late[5,1:6] <- round(SavDick$`Yes Late`[1:6],3)

Table_contents_late[5,7] <- '-'

Table_rownames <- c('Comp. Type','Mean (Empirical)','Mean (Inferred)','SD','BF')
row.names(Table_contents_late) <- Table_rownames

# Try out trial-by-trial descriptive adequacy plot ------------------------
# par(mfrow = c(7,7),
#           oma = c(5,4,0,0) + 0.1,
#           mar = c(0,0,1,1) + 0.1)
# 
# for(k in 1:7){
#   for(j in 1:7){
#   participant_ticker = ((k-1)*7) + j
#   ypred_lo = ((participant_ticker*6)-5)+nBeta
#   ypred_hi = (participant_ticker*6)+nBeta
#   y_lo = (participant_ticker*6)-5
#   y_hi = (participant_ticker*6)
#   tempplot = YyesEarly[y_lo:y_hi]
#   trial_index = OnlineData$CompType[y_lo:y_hi]
#   tempplot_hi = rawsummaryyesEarly[ypred_lo:ypred_hi,7]
#   tempplot_lo = rawsummaryyesEarly[ypred_lo:ypred_hi,3]
#   tempplot <- tempplot[trial_index]
#   tempplot_hi <- tempplot_hi[trial_index]
#   tempplot_lo <- tempplot_lo[trial_index]
#   plot(1:6,tempplot,ylim = c(0,1),xaxt='n',yaxt='n')
#   for(i in 1:6){
#     arrows(x0 = i, y0 = tempplot_hi[i],
#            x1 = i, tempplot_lo[i],code = 3,
#            col = 'black',angle = 90, length = .1,lty = 3)
#   }
#   }
# }


# Illustrative graph for comp types ---------------------------------------

par(mfrow=c(1,1))
EarlyGame_Ill = c(.70,.63,.55,.55,.63,.70)
LateGame_Ill = c(.55,.63,.70,.70,.63,.55)
par(mar = c(3, 3, 1, 2))

plot(1:6,TVJ_Early_bycomp,ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'black',bty = 'n',xlim = c(.75,6),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 23)

points(1:6,TVJ_Late_bycomp,col='black',pch = 15,cex = 2)

title(ylab="Rate of Endorsement", line=1.5, cex.lab=1.8)
title(xlab="Composition Type", line=1.6, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.5,lwd = 2)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,6.25),y=c(0,0),lwd =2)
abline(v=0)

for(i in 1:5){
  lines(x = c(i,i+1),y = c(TVJ_Early_bycomp[i],TVJ_Early_bycomp[i+1]),lty = 5,lwd = 2)
  lines(x = c(i,i+1),y = c(TVJ_Late_bycomp[i],TVJ_Late_bycomp[i+1]),lty = 3,lwd = 2)
}

legend(3.75,.15,legend = c('Early-Game Generalization','Late-Game Generalization'),
      pch = c(23,15),pt.cex = c(2.2,2),lwd =c(2,2),lty = c(5,3),
      col = c('black','black'),
      cex = 1.3,bty = 'n')


# Alex version ------------------------------------------------------------

par(mfrow=c(1,2))
EarlyGame_Ill = c(.70,.63,.55,.55,.63,.70)
LateGame_Ill = c(.55,.63,.70,.70,.63,.55)
par(mar = c(3, 4, 2, 2))

plot(1:6,TVJ_Early_bycomp,ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'black',bty = 'n',xlim = c(.75,6),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 23)

title(main="Early Game ", line=.5, cex.main=1.8)
title(ylab="Endorsement Rate", line=2.1, cex.lab=1.8)
title(xlab="Composition Type", line=1.6, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.5,lwd = 2)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,6.25),y=c(0,0),lwd =2)
abline(v=0)

plot(1:6,TVJ_Late_bycomp,ylim = c(0,1),main = '', xlab = ' ',ylab = '',
     col = 'black',bty = 'n',xlim = c(.75,6),frame.plot = FALSE,
     axes = FALSE,cex = 2.2,lwd =2,pch = 23)

title(main="Late Game", line=.5, cex.main=1.8)
title(xlab="Composition Type", line=1.6, cex.lab=1.8)

axis(side=1,pos=0,cex.axis = 1.5,lwd = 2)
axis(side = 2, pos=.75,cex.axis = 1.5,lwd = 2)

lines(x = c(.75,6.25),y=c(0,0),lwd =2)
abline(v=0)

# Run all data through NormCDF JAGS model ---------------------------------
## Assumes all expert generic statements are produced by the same process

Allxy <- matrix(NaN, nrow = nExpResponses,ncol = 6)

for(i in 1:nExpResponses){
  Allxy[i,1] <- OnlineData$Participant[i]
  Allxy[i,2] <- OnlineData$EarlyPrior[i]
  Allxy[i,3] <- OnlineData$EarlyTVJ[i]
  Allxy[i,4] <- OnlineData$LatePrior[i]
  Allxy[i,5] <- OnlineData$LateTVJ[i]
}

Allxy <- data.frame(Allxy)

names(Allxy) <- c('participant','early_x','early_y','late_x','late_y')

x <- c(Allxy$early_x, Allxy$late_x)
y <- c(Allxy$early_y, Allxy$late_y)
n <- length(x)

for(i in 1:n){
  if(y[i] == TRUE){
    y[i] <- 1
  }
  if(y[i] == FALSE){
    y[i] <- 0
  }
}

data <- list("x","y", "n") # to be passed on to JAGS
myinits <- list(
  list("mu" = runif(1,0,1), "sigma" = .5))

# parameters to be monitored:	
parameters <- c("mu", "sigma")

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, inits=myinits, parameters,
                model.file="NormCDF.txt", n.chains=1,n.burnin = 1000, n.iter=4000,  n.thin=2, DIC=T)

#### For running multiple chains in parallel. Involves different jags call and initialization####
#Ways this works is that when running chains in parallel, R basically initializes multiple jags calls with
#a single chain, so you build an initialization function to create new initializations
#each time R initiates the separate, parallel chain
inits.jagsParallel=function(){
  return(list(list("mu" = runif(1,0,1),"sigma" = .5)))
}

#Now use jags.parallel to run multiple chains much quicker, adjust chains in n.chains

samples =jags.parallel(data,inits=inits.jagsParallel(), parameters.to.save = parameters,
                       model.file="NormCDF.txt",n.chains=4,n.iter=4000,n.thin = 2)

# End parallel running section

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

AllMu    <- samples$BUGSoutput$sims.list$mu
AllSigma <- samples$BUGSoutput$sims.list$sigma

chainlength = length(AllMu)/4

hist(AllMu)
hist(AllSigma)
plot(1:chainlength,AllMu[1:chainlength],type = 'l', main = 'All Mu Chains')
colors <- c('blue','red','green')
for(i in 1:3){
  low = (i*chainlength)+1
  hi = (i+1)*chainlength
  lines(1:chainlength,AllMu[low:hi],col = colors[i])
}

plot(1:chainlength,AllSigma[1:chainlength],type = 'l', main = 'All Sigma Chains')
colors <- c('blue','red','green')
for(i in 1:3){
  low = (i*chainlength)+1
  hi = (i+1)*chainlength
  lines(1:chainlength,AllSigma[low:hi],col = colors[i])
}

AllStats <- matrix(NaN,nrow = 2, ncol = 3)
AllStats <- data.frame(AllStats)
names(AllStats) <- c('mode','Conf_5','Conf_95')
row.names(AllStats) <- c('mu','sigma')

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

AllStats$mode[1] <- getmode(AllMu)
AllStats$mode[2] <- getmode(AllSigma)

AllMuSort <- sort(AllMu)
AllSigmaSort <- sort(AllSigma)
AllStats$Conf_5[1] <- AllMuSort[200]
AllStats$Conf_95[1] <- AllMuSort[3800]
AllStats$Conf_5[2] <- AllSigmaSort[200]
AllStats$Conf_95[2] <- AllSigmaSort[3800]

q = seq(0,1,by = .01)

plot(q,pnorm(q,AllStats$mode[1],AllStats$mode[2]), type = 'l',main = 'All Responses',
     ylab = 'Prob of True Response',xlab = 'Estimated Prevalence')
lines(q,pnorm(q,AllStats$mode[1],AllStats$Conf_5[2]), col = 'red', lty = 2)
lines(q,pnorm(q,AllStats$mode[1],AllStats$Conf_95[2]), col = 'red', lty = 2)
lines(q,pnorm(q,AllStats$Conf_5[1],AllStats$mode[2]), col = 'blue',lty = 3)
lines(q,pnorm(q,AllStats$Conf_95[1],AllStats$mode[2]), col = 'blue', lty = 3)
legend(.8, .3, col = c('black','red','blue'),
       legend = c('mode','5-95% Sigma','5-95% Mu'),lty = c(1,2,3))
abline(v = .5,lty = 4)
abline(h = .5, lty = 4)


# NormCDF assuming late and early have different processes ----------------
## Assumes early and late statements use different processes, but all participants use same for each

for(time in 1:2){
  xTicker = time*2
  yTicker = (time*2)+1
  x <- Allxy[,xTicker]
  y <- Allxy[,yTicker]
  n <- length(x)
  
  for(i in 1:n){
    if(y[i] == TRUE){
      y[i] <- 1
    }
    if(y[i] == FALSE){
      y[i] <- 0
    }
  }
  
  data <- list("x","y", "n") # to be passed on to JAGS
  myinits <- list(
    list("mu" = runif(1,0,1), "sigma" = .5))
  
  # parameters to be monitored:	
  parameters <- c("mu", "sigma")
  
  # The following command calls JAGS with specific options.
  #This is for running 1 chain (use code below for faster multiple chains)
  samples <- jags(data, inits=myinits, parameters,
                  model.file="NormCDF.txt", n.chains=1, n.iter=3000, 
                  n.burnin=1000, n.thin=1, DIC=T)
  
  #### For running multiple chains in parallel. Involves different jags call and initialization####
  #Ways this works is that when running chains in parallel, R basically initializes multiple jags calls with
  #a single chain, so you build an initialization function to create new initializations
  #each time R initiates the separate, parallel chain
  inits.jagsParallel=function(){
    return(list(list("mu" = runif(1,0,1),"sigma" = .5)))
  }
  
  #Now use jags.parallel to run multiple chains much quicker, adjust chains in n.chains
  
  samples =jags.parallel(data,inits=inits.jagsParallel(), parameters.to.save = parameters,
                         model.file="NormCDF.txt",n.chains=4,n.burnin=1000,n.iter=3000,n.thin = 2)
  
  # End parallel running section
  
  # Now the values for the monitored parameters are in the "samples" object, 
  # ready for inspection.
  
  if(time == 1){
    EarlyMu    <- samples$BUGSoutput$sims.list$mu
    EarlySigma <- samples$BUGSoutput$sims.list$sigma
  }
  if(time == 2){
    LateMu    <- samples$BUGSoutput$sims.list$mu
    LateSigma <- samples$BUGSoutput$sims.list$sigma
  }
}

chainlength = length(EarlyMu)/4

hist(EarlyMu)
hist(EarlySigma)
plot(1:chainlength,EarlyMu[1:chainlength],type = 'l', main = 'Early Mu Chains')
colors <- c('blue','red','green')
for(i in 1:3){
  low = (i*chainlength)+1
  hi = (i+1)*chainlength
  lines(1:chainlength,EarlyMu[low:hi],col = colors[i])
}

plot(1:chainlength,EarlySigma[1:chainlength],type = 'l', main = 'Early Sigma Chains')
colors <- c('blue','red','green')
for(i in 1:3){
  low = (i*chainlength)+1
  hi = (i+1)*chainlength
  lines(1:chainlength,EarlySigma[low:hi],col = colors[i])
}

EarlyStats <- matrix(NaN,nrow = 2, ncol = 3)
EarlyStats <- data.frame(EarlyStats)
names(EarlyStats) <- c('mode','Conf_5','Conf_95')
row.names(EarlyStats) <- c('mu','sigma')

EarlyStats$mode[1] <- mean(EarlyMu)
EarlyStats$mode[2] <- mean(EarlySigma)

EarlyMuSort <- sort(EarlyMu)
EarlySigmaSort <- sort(EarlySigma)
EarlyStats$Conf_5[1] <- EarlyMuSort[200]
EarlyStats$Conf_95[1] <- EarlyMuSort[3800]
EarlyStats$Conf_5[2] <- EarlySigmaSort[200]
EarlyStats$Conf_95[2] <- EarlySigmaSort[3800]

chainlength = length(LateMu)/4

hist(LateMu)
hist(LateSigma)
plot(1:chainlength,LateMu[1:chainlength],type = 'l', main = 'Late Mu Chains')
colors <- c('blue','red','green')
for(i in 1:3){
  low = (i*chainlength)+1
  hi = (i+1)*chainlength
  lines(1:chainlength,LateMu[low:hi],col = colors[i])
}

plot(1:chainlength,LateSigma[1:chainlength],type = 'l', main = 'Late Sigma Chains')
colors <- c('blue','red','green')
for(i in 1:3){
  low = (i*chainlength)+1
  hi = (i+1)*chainlength
  lines(1:chainlength,LateSigma[low:hi],col = colors[i])
}

LateStats <- matrix(NaN,nrow = 2, ncol = 3)
LateStats <- data.frame(LateStats)
names(LateStats) <- c('mode','Conf_5','Conf_95')
row.names(LateStats) <- c('mu','sigma')

LateStats$mode[1] <- mean(LateMu)
LateStats$mode[2] <- mean(LateSigma)

LateMuSort <- sort(LateMu)
LateSigmaSort <- sort(LateSigma)
LateStats$Conf_5[1] <- LateMuSort[200]
LateStats$Conf_95[1] <- LateMuSort[3800]
LateStats$Conf_5[2] <- LateSigmaSort[200]
LateStats$Conf_95[2] <- LateSigmaSort[3800]


q = seq(0,1,by = .01)

plot(q,pnorm(q,EarlyStats$mode[1],EarlyStats$mode[2]), type = 'l',main = 'Early Responses',
     ylab = 'Prob of True Response',xlab = 'Estimated Prevalence')
lines(q,pnorm(q,EarlyStats$mode[1],EarlyStats$Conf_5[2]), col = 'red', lty = 2)
lines(q,pnorm(q,EarlyStats$mode[1],EarlyStats$Conf_95[2]), col = 'red', lty = 2)
lines(q,pnorm(q,EarlyStats$Conf_5[1],EarlyStats$mode[2]), col = 'blue', lty = 3)
lines(q,pnorm(q,EarlyStats$Conf_95[1],EarlyStats$mode[2]), col = 'blue', lty = 3)
legend(.8, .3, col = c('black','red','blue'),
       legend = c('mean','5-95% Sigma','5-95% Mu'),lty = c(1,2,3))
abline(h=.5,v=.5,lty=4)
abline(v = mean(OnlineData$EarlyPrior),col = 'green')

plot(q,pnorm(q,LateStats$mode[1],LateStats$mode[2]), type = 'l', main = 'Late Responses',
     ylab = 'Prob of True Response',xlab = 'Estimated Prevalence')
lines(q,pnorm(q,LateStats$mode[1],LateStats$Conf_5[2]), col = 'red',lty = 2)
lines(q,pnorm(q,LateStats$mode[1],LateStats$Conf_95[2]), col = 'red',lty = 2)
lines(q,pnorm(q,LateStats$Conf_5[1],LateStats$mode[2]), col = 'blue',lty = 3)
lines(q,pnorm(q,LateStats$Conf_95[1],LateStats$mode[2]), col = 'blue', lty = 3)
legend(.8, .3, col = c('black','red','blue'),
       legend = c('mean','5-95% Sigma','5-95% Mu'),lty = c(1,2,3))
abline(h=.5,v=.5,lty=4)
abline(v = mean(OnlineData$LatePrior),col = 'green')

plot(q,pnorm(q,.5,.001),type = 'l')
lines(q,pnorm(q,.5,.05),lty=2, col = 'red')
abline(h = .5, v = .5,lty = 4)

DemoTrue <- pnorm(q,.55,.05)
DemoFalse <- 1-pnorm(q,.45,.05)
DemoNothing <- rep(NaN,length(DemoTrue))
for(i in 1:length(DemoTrue)){
  DemoNothing[i] <- 1-(DemoTrue[i]+DemoFalse[i])
}


plot(q,DemoTrue,type = 'l',xlab = 'Estimated Prevalence',ylab = 'Prob of Response')
lines(q,DemoFalse, col = 'red',lty = 2)
lines(q,DemoNothing, col = 'blue',lty = 3)
legend (.8,.5,col = c('black','red','blue'),lty = c(1,2,3),legend = c('True','False','Nothing'))