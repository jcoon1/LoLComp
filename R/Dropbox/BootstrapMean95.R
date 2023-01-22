BootstrapMean95 = function(dataVector, nSamples, nIterations){
  nData <- length(dataVector)
  collectMeans <- rep(NaN,nIterations)
  collectSamples <- rep(NaN,nSamples)
  rawRandom <- rep(NaN,nSamples)
  
  nSamples <- nData
  
  for(i in 1:nIterations){
    for(j in 1:nSamples){
    rawRandom[j] <- runif(1)
    sampleIndex <- round(rawRandom*nData)
    collectSamples <- dataVector[sampleIndex]
    collectMeans[i] <- mean(collectSamples)
  }
}
  lo = .025*nIterations
  hi = .975*nIterations
  sortMeans = sort(collectMeans)
  bootRange <- c(sortMeans[lo],sortMeans[hi])
  return(bootRange)
}