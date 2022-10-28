model <- function(data.subset){
  
  ## change command names for MASS and dplyr 
  dselect <- dplyr::select
  dfilter <- dplyr::filter
  ### Call the function firstStage to recover alpha_q and Phi_t ### 
  
  y <- data.subset %>% dselect(lpatient_years)
  X <- data.subset %>% dselect(lhires,lstations, lstaff, inspection_rate, time_since_survey)
  X.shift <- data.subset %>% dselect(hirepos, for_profit, compLevel, davita, fresenius)
  
  first.stage.reg <- firstStage(y,X,X.shift,as.matrix(data.subset$q.res), as.matrix(data.subset$q.iv.res)) 
  alpha.q <- first.stage.reg[[1]]
  Phi.t <- first.stage.reg[[2]]

  ### no quality case Phi = y - y.hat ## 
  Phi.hat.nq <- locReg(y,X,X.shift)
  
  Phi.t.nq <- Phi.hat.nq[[1]]
  Phi.t.nq[Phi.hat.nq[[2]] > 0 ] <- NaN 
  
 
  ### Add the recovered Phi.t and Phi.t.nq to the original database 
  data.subset <- cbind(data.subset, Phi.t, Phi.t.nq)
  data.subset <- data.subset %>% filter(!(is.na(Phi.t)))
  
  
  ### Second Stage to estimate beta and OLS/FE parameters [Create a function for this later] ### 
  second.stage.reg <- secondStage(alpha.q,data.subset)
  
  tab5.output <- rbind(as.matrix(alpha.q),second.stage.reg) 
  
  return(tab5.output) 
  
  
}