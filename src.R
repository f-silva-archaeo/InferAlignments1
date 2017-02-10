################################
### Skyscape Alignment Model ###
################################
### PARAMETERS:
## target    :: target orientation (in azimuth or declination)
## N         :: sample size
## stage1    :: stage 1 uncertainty, or total deviation if stage2 and stage 3 left empty
## stage2    :: (optional), stage 2 uncertainty, leave empty for simplified model
## stage3    :: (optional), stage 3 uncertainty, leave empty for simplified model
#################################

SimSkyscape  = function (target, N, stage1, stage2 = NaN, stage3 = NaN) {
  if (is.nan(stage2)) {
    # Simplified version
    noise <- stage1
  } else {
    # Full version
    noise <- sqrt(stage1^2+stage2^2+stage3^2)
  }
  
  out <- rnorm(N, mean=target, noise) 
  return(out)
}

###################################################################################################



#################################
##### Skyscape Curvigram Analysis #####
#################################
### PARAMETERS:
## data      :: data (in azimuth or declination)
## mes.unc   :: estimated measurement uncertainty
#################################

Curv_Skyscape = function (data, mes.unc) {
  require(numDeriv)
  require(rootSolve)
  require(MESS)
  
  xx <- seq(-90, 90, 0.01)
  N <- NROW(data)
  
  # Create SPD/curvigram
  spd <- array(0,NROW(xx))
  for (i in 1:N) {
    spd <- spd + dnorm(xx,data[i],mes.unc)
  }
  
  # Identifying maxima
  func0 <- splinefun(xx,spd)
  fd <- grad(func0,xx); func1 <- splinefun(xx,fd)
  sd <- grad(func1,xx); func2 <- splinefun(xx,sd)
  roots <- uniroot.all(func1,interval=c(min(xx),max(xx)))
  ind <- which(func2(roots)<0)
  maxima <- roots[ind]
  f.max <- func0(maxima)
  
  ind <- sort(f.max, index.return=T, decreasing = T)$ix
  
  # Uncertainty
  if (NROW(maxima)==1) {
    prec <- SPD_Prec(data, mes.unc)
  } else { prec <- c(NA,NA)}
  
  # Output
  output <- list(x = xx, y = spd/auc(xx,spd), peak = maxima[ind], unc = prec)
  return(output)
}

###################################################################################################



#################################
##### Skyscape ML Analysis #####
#################################
### PARAMETERS:
## data      :: data (in azimuth or declination)
## mes.unc   :: estimated measurement uncertainty
#################################

ML_Skyscape = function (data, mes.unc) {

  xx <- seq(mean(data)-10*sd(data),mean(data)+10*sd(data),0.01)
  N <- NROW(data)
  
  # MLE
  m <- weighted.mean(data,rep(mes.unc,N))
  sdom <- sd(data)/sqrt(N)
  
  # Output
  ml <- dnorm(xx, m, sdom)
  output <- list(x = xx, y = ml, peak = m, unc = 2*sdom)
  return(output)
}

###################################################################################################


################################
#### Estimate Curvigram Precision ####
################################
### PARAMETERS:
## data      :: data (in azimuth or declination)
## mes.unc   :: estimated measurement uncertainty
#################################

Curv_Prec = function (data, mes.unc) {
  
  N <- NROW(data)
  
  # Load fit
  load("sims.1-target.RData")
  
  prec <- predict(fit5, newdata=list(x=log(N), y=sd(data), z=mes.unc))
  
  return(prec)
}

###################################################################################################


################################
#### Estimate ML Precision ####
################################
### PARAMETERS:
## data      :: data (in azimuth or declination)
## mes.unc   :: estimated measurement uncertainty
#################################

ML_Prec = function (data, mes.unc) {
  
  N <- NROW(data)
  
  prec <- 2* sd(data)/sqrt(N)
  
  return(prec)
}

###################################################################################################