setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Illustrative graph --------------------------------------------------

par(mfrow=c(1,1))

par(mar = c(3.5, 3, 1, 2))

q = seq(0,1,by = .01)
q_odds <- q/(1-q)
q_logodds <- log(q_odds)

temp_hierarchical <- q



thresholds <- c(seq(.5,.99,.01), seq(.99,.01,-.01), seq(.01,.5,.01))

for(i in 1:length(thresholds)){
  temp_file <- paste("alpha_demo",i,".png",sep ="")
  png(temp_file,width = 902,height = 573) 
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
  
  lines(x = c(0,1),y = c(.5,.5),lty = 4)
  lines(x = c(.5,.5),y = c(0,1),lty = 4)

  example_threshold <- thresholds[i]
  example_alpha <- log(example_threshold/(1-example_threshold))
  example_beta <- 5
  temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
  lines(q,temp_hierarchical,type = 'l',lwd = 2)

  text(x = .2, y = .7,labels = expression(paste(alpha, " = ", "  ")), cex = 1.5)
  if(round(example_alpha,2) == 0){
    text(x = .22, y = .7,labels = round(example_alpha,2), cex = 1.5)
  }
  else{
  text(x = .24, y = .7,labels = round(example_alpha,2), cex = 1.5)
  }
  text(x = .2, y = .64,labels = expression(paste(beta, " = ", 5)), cex = 1.5)

  dev.off() 

}

wump = 40
q_beta <- wump:(100+wump)
betas_temp <- (q_beta/wump)^2
closest <- min(abs(betas_temp-5))
betas_temp_id <- which(abs(betas_temp-5)==closest)
betas_temp_id
betas_index <- c(seq(betas_temp_id,1,-1),seq(1,length(betas_temp),1),seq(length(betas_temp),betas_temp_id,-1))
betas <- betas_temp[betas_index]

plot(1:length(betas),betas)
abline(h=5)
for(i in 1:length(betas)){
  temp_file <- paste("beta_demo",i,".png",sep ="")
  png(temp_file,width = 902,height = 573) 
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
  
  lines(x = c(0,1),y = c(.5,.5),lty = 4)
  lines(x = c(.5,.5),y = c(0,1),lty = 4)
  
  example_threshold <- .5
  example_alpha <- log(example_threshold/(1-example_threshold))
  example_beta <- betas[i]
  temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
  lines(q,temp_hierarchical,type = 'l',lwd = 2)
  
  text(x = .2, y = .7,labels = expression(paste(alpha, " = ", 0)), cex = 1.5)


  
  text(x = .2, y = .64,labels = expression(paste(beta, " = ", "  ")), cex = 1.5)
  text(x = .24, y = .64,labels = round(example_beta,2), cex = 1.5)
  
  dev.off() 
  
}

png("illustrative0.png",width = 902,height = 573) 
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

dev.off() 


png("illustrative1.png",width = 902,height = 573) 
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

lines(x = c(0,1),y = c(.5,.5),lty = 4)
lines(x = c(.5,.5),y = c(0,1),lty = 4)

dev.off() 
# 
# example_threshold <- thresholds[i]
# example_alpha <- log(example_threshold/(1-example_threshold))
# example_beta <- 45
# temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
# lines(q,temp_hierarchical,type = 'l',lwd = 2,col="tomato",lty = 2)
# 
# example_threshold <- .5
# example_alpha <- log(example_threshold/(1-example_threshold))
# example_beta <- 1
# temp_hierarchical <- 1/(1 + exp(-example_beta*(q_logodds-example_alpha)))
# lines(q,temp_hierarchical,type = 'l',lwd = 2,col="tomato",lty = 2)
