modelTab6 <- function(data.subset){
  
  ## change command names for MASS and dplyr 
  dselect <- dplyr::select
  dfilter <- dplyr::filter
  
  
  y <- data.subset %>% dselect(lpatient_years)
  X <- data.subset %>% dselect(lhires,lstations, lstaff, inspection_rate, time_since_survey)
  X.shift <- data.subset %>% dselect(hirepos, for_profit, compLevel, davita, fresenius)
  
  ### Table 6 results ###
  
  ## Col 1: replace q.iv <- q (estimate coefficient assuming no measurement error, should see attenuation bias)
  
  first.stage.reg.col1 <- firstStage(y,X,X.shift,as.matrix(data.subset$q.res), as.matrix(data.subset$q.res))
  alpha.q.col1 <- first.stage.reg.col1[[1]]
  Phi.t <- first.stage.reg.col1[[2]]
  

  
  ### Add the recovered Phi.t to the original database
  data.subset.col1 <- cbind(data.subset, Phi.t)
  data.subset.col1 <- data.subset.col1 %>% filter(!(is.na(Phi.t)))
  
  ## Second Stage to estimate beta 
  second.stage.reg.col1 <- secondStageTab6(alpha.q.col1,data.subset.col1)
  
  
  #############################################################################################################
  ## Col 2: Do not control for patient characteristics. Use infect_rate as quality and death_ratio as instrument. 
  
  first.stage.reg.col2 <- firstStage(y,X,X.shift,as.matrix(data.subset$infect_rate), as.matrix(data.subset$death_ratio))
  alpha.q.col2 <- first.stage.reg.col2[[1]]
  Phi.t <- first.stage.reg.col2[[2]]
  

  ### Add the recovered Phi.t and Phi.t.nq to the original database
  data.subset.col2 <- cbind(data.subset, Phi.t)
  data.subset.col2 <- data.subset.col2 %>% filter(!(is.na(Phi.t)))
  
  
  ### Second Stage to estimate beta parameters [Create a function for this later] ###
  second.stage.reg.col2 <- secondStageTab6(alpha.q.col1,data.subset.col2)
  
  ################################################################################################################
  ## col 3: Get rid of competition indicator from X.shift just modify X.shift and perform the same computations 
  

  X.shift <- data.subset %>% dselect(hirepos, for_profit, davita, fresenius)
  
  first.stage.reg.col3 <- firstStage(y,X,X.shift,as.matrix(data.subset$q.res), as.matrix(data.subset$q.iv.res)) 
  alpha.q.col3 <- first.stage.reg.col3[[1]]
  Phi.t <- first.stage.reg.col3[[2]]
  
  ### Add the recovered Phi.t and Phi.t.nq to the original database 
  data.subset.col3 <- cbind(data.subset, Phi.t)
  data.subset.col3 <- data.subset.col3 %>% filter(!(is.na(Phi.t)))
  
  
  ### Second Stage to estimate beta and OLS/FE parameters [Create a function for this later] ### 
  second.stage.reg.col3 <- secondStageTab6(alpha.q.col3,data.subset.col3)
  
  tab6.output <- rbind(alpha.q.col1,as.matrix(second.stage.reg.col1[1:2,]), alpha.q.col2, as.matrix(second.stage.reg.col2[1:2,]),alpha.q.col3, as.matrix(second.stage.reg.col3[1:2,]))
  
  return(tab6.output)
}