###########################################################
## Model evaluation metrics
## From Markus Reichstein
###########################################################
modelEval  <- function(obs, mod) {
  res <- obs - mod
  RMSE  <- sqrt(mean(res^2))
  r <- cor(obs, mod)
  r2  <- r^2
  
  meanObs  <- mean(obs)
  meanMod <- mean(mod)
  sdObs <- sd(obs)
  sdMod <- sd(mod)
  
  MEF  <- 1 - sum(res^2)/sum((obs-meanObs)^2)
  
  bias <- meanMod - meanObs
  bias_sq  <- (meanMod - meanObs)^2
  var_err <- (sdObs-sdMod)^2
  phase_err <- (1-r)*2*sdObs*sdMod
  
  c(RMSE=RMSE, r=r, r2=r2, MEF=MEF, bias=bias, bias_sq=bias_sq, var_err=var_err, phase_err=phase_err)
}
