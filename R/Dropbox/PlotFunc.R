MainPlot = function(Novice_inf,Novice_emp,Novice_empCon,Expert_inf,Expert_emp,
                    MainTitle,if_legend,leftmost,bottommost){
  plot(1:6,Expert_inf[1:6,1],ylim = c(0,1),main = MainTitle, xlab = ' ',ylab = '',
       col = rgb(red = 220/255, green = 20/255, blue = 60/255),bty = 'n',xlim = c(.75,6),frame.plot = FALSE,
       axes = FALSE,cex = 2,lwd =2,pch = 10)
  if(leftmost==1){
    title(ylab="Interpretation ( /100 games)", line=1.5, cex.lab=1.8)
  }
  
  if(bottommost==1){
    title(xlab="Composition Type", line=1.8, cex.lab=1.8)
  }
  
  xlabels = c('E+','E0','E-','L+','L0','L-')
  xat = c(1:6)
  axis(side=1,pos=0,cex.axis = 1.6,lwd = 2,labels = xlabels,at = xat)
  axis(side = 2, pos=.75,cex.axis = 1.6,lwd = 2)
  lines(x= c(.75,8),y = c(Novice_inf[1],Novice_inf[1]),col = rgb(red = 51/255, green = 153/255, blue = 0/255),lwd=2)
  lines(x= c(.75,8),y = c(Novice_inf[2],Novice_inf[2]),col = rgb(red = 51/255, green = 153/255, blue = 0/255),lty = 2,lwd=2)
  lines(x= c(.75,8),y = c(Novice_inf[3],Novice_inf[3]),col = rgb(red = 51/255, green = 153/255, blue = 0/255),lty = 2,lwd=2)
  for(i in 1:6){
    arrows(x0 = i, y0 = Expert_inf[i,2],
           x1 = i, y1 = Expert_inf[i,3],code = 3,
           col = rgb(red = 220/255, green = 20/255, blue = 60/255),
           angle = 90, length = .1,lwd = 2)
  }
  
  points(1.15:6.15,Expert_emp[1:6],col='black',pch = 16,cex = 2)
  if(if_legend==1){
  legend(3.1,.37,legend = c('Naive (Inferred)','Experienced (Inferred)','Naive by Condition (Empirical)',
                         'Naive Overall (Empirical)','Experienced (Empirical)'),
         pch = c(NA,10,17,NA,16),lty = c(1,NA,NA,3,NA),lwd = c(2,2,NA,3,NA),pt.cex = c(NA,2,2,NA,2),
         col = c(rgb(red = 51/255, green = 153/255, blue = 0/255),
                                                              rgb(red = 220/255, green = 20/255, blue = 60/255),
                                                              '#7B68EE',col=rgb(red=0,green=0,blue=0),
                                                              'black'),
         cex = 1.3,bty = 'n')
  }
  lines(x= c(.75,8),y=c(Novice_emp,Novice_emp),col=rgb(red=0,green=0,blue=0),lwd=3,lty=3)
  
  points(.85:5.85,Novice_empCon[1:6],col='#7B68EE',pch = 17,cex = 2)
  
  lines(x = c(.75,6.25),y=c(0,0),lwd =2)
  abline(v=0)
  
  
  return(recordPlot())
}