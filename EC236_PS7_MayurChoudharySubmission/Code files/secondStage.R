secondStage <- function(alpha.q,data.subset){
  
  
  ### Subset data to get rid of observations for which lag is not available 
  ### Note that the subsetting works in such a way that I am left with very few observations (4.5k instead of around 6k) and hence the results vary somewhat from the original paper
  ### Uncomment the subestting of data.subset for removing NAs in Phi.t  resolves this problem and we are left with around 6k observations. 
  
  for.id <- data.subset %>% dselect(for_id)
  for.id <- unique(for.id)
  names(for.id) <- 'row_id'
  
  data.subset.t <- inner_join(data.subset, for.id, by = 'row_id')
  data.subset.t <- data.subset.t %>% filter(!(is.na(q.res)))
  
  ### Get the lagged Phi data ### 
  cur.id <- data.subset.t %>% dselect(row_id)
  names(cur.id) <- 'for_id'
  data.subset.lag <- inner_join(data.subset,cur.id, by = 'for_id')
  data.subset.lag <- unique(data.subset.lag) ## get rid of duplicates due to inner join 
  
  
  ### Will need to modify this for table 8 
  data.subset.t <- data.subset.t %>% mutate(tf.t = lpatient_years - alpha.q[1,1]*q.res)  ### The transformation function (LHS = y + alpha*q) 
  tf.t <- as.matrix(data.subset.t$tf.t)
  data.subset.t <- data.subset.t %>% mutate(tf.t.nq = lpatient_years) ### No quality case for tab 5
  tf.t.nq <- as.matrix(data.subset.t$tf.t.nq)
    
  ### Get lagged and time t data for input to the GMM function 
  factors.t  <- data.subset.t %>% dselect(lstations, lstaff)
  factors.lag <- data.subset.lag %>% dselect(lstations,lstaff)
  phi.lag <- data.subset.lag %>% dselect(Phi.t)
  phi.lag.nq <- data.subset.lag %>% dselect(Phi.t.nq)
  
  
  ### Create the types for non-parametric g function for table 6 column 3 
  type.data <- data.subset.lag %>% dselect(for_profit, davita, fresenius)
  type.unique <- unique(type.data)
  type.unique$type.id <- matrix(seq(1,nrow(type.unique)), nrow(type.unique), 1)
  
  type.data <- left_join(type.data, type.unique, by = c('for_profit','davita','fresenius'), by.x = TRUE)
  
  ### GMM optimization ### 
  # initialize weight matrix 
  W <- diag(1,ncol(factors.t), ncol(factors.t)) ## 2 moment conditions for 2 factors
  
  # initialize parameter beta 
  beta <- matrix(0.5, ncol(factors.t),1)
  second.stage.reg <- optim(beta, gmmFun , tf.t = tf.t , factors.t = factors.t , phi.lag = phi.lag , factors.lag = factors.lag , W = W, order = 4, type = NULL, method = "Nelder-Mead")
  beta.tab5 <- as.matrix(second.stage.reg$par) 
  
  # No quality estimate 
  second.stage.reg.nq <- optim(beta, gmmFun , tf.t = tf.t , factors.t = factors.t , phi.lag = phi.lag.nq , factors.lag = factors.lag , W = W, order = 4, type = NULL, method = "Nelder-Mead")
  beta.tab5.nq <- as.matrix(second.stage.reg.nq$par) 
  
  # Estimates for Table 7 column 3 (non-parametric)
  second.stage.np.reg <- optim(beta, gmmFun , tf.t = tf.t , factors.t = factors.t , phi.lag = phi.lag , factors.lag = factors.lag , W = W, order = 4, type = as.matrix(type.data$type.id), method = "Nelder-Mead")
  beta.np <- as.matrix(second.stage.np.reg$par) 
  
  # Estimates for Table 7 column 2 (parametric) 
  type.data <- type.data %>% mutate(type.1 = ifelse(type.id == 1,1,0), type.2 = ifelse(type.id == 2,1,0), type.3 = ifelse(type.id == 3,1,0), type.4 = ifelse(type.id == 4,0,1))
  type.dummy <- type.data %>% dselect(type.1, type.2, type.3, type.4) 
  
  second.stage.para.reg <- optim(beta, gmmFunPara , tf.t = tf.t , factors.t = factors.t , phi.lag = phi.lag , factors.lag = factors.lag , W = W, order = 4, type = as.matrix(type.dummy), method = "Nelder-Mead")
  beta.p <- as.matrix(second.stage.para.reg$par) 
  
  
  #### OLS and FE regressions ##### 
  data.ols <- data.subset.t %>% dselect(lpatient_years,lstations,lstaff,q.res,cms_code)
  
  ols.reg <- lm(lpatient_years ~  q.res + lstations + lstaff, data = data.ols)
  beta.ols <- as.matrix(coef(ols.reg))
  
  ols.reg.nq <- lm(lpatient_years ~ lstations + lstaff, data = data.ols)
  beta.ols.nq <- as.matrix(coef(ols.reg.nq)) 
  
  fe.reg <- plm(lpatient_years ~  q.res + lstations + lstaff, index = c('cms_code'), model = c('within'), data = data.ols)
  beta.fe <- as.matrix(coef(fe.reg))
  
  fe.reg.nq <- plm(lpatient_years ~ lstations + lstaff, index = c('cms_code'), model = c('within'), data = data.ols)
  beta.fe.nq <- as.matrix(coef(fe.reg.nq))

  second.stage.output <- rbind(beta.tab5, beta.tab5.nq, beta.ols, beta.ols.nq, beta.fe, beta.fe.nq, beta.np, beta.p)  
  return(second.stage.output)
  
  ### Note here I comment out the code without this subsetting and find that the coefficient on quality is not as low as in the paper. It is in fact in line with the one from the model. This is a curious fact.
  ### Please uncomment and check that the above assertion is indeed true, this implies that we need to think about why there is this discrepancy between OLS and FE regression estimates just on subsetting. Not filtering data.
  
  
  # data.ols <- data.phi.subset %>% dselect(lpatient_years,lstations,lstaff,q.res,cms_code)
  # ols.reg <- lm(lpatient_years ~  q.res + lstations + lstaff, data = data.ols)
  # beta.ols <- coef(ols.reg)
  # 
  # ols.reg.nq <- lm(lpatient_years ~ lstations + lstaff, data = data.ols)
  # beta.ols.nq <- coef(ols.reg.nq) 
  # 
  # fe.reg <- plm(lpatient_years ~  q.res + lstations + lstaff, index = c('cms_code'), model = c('within'), data = data.ols)
  # beta.fe <- coef(fe.reg)
  # 
  # fe.reg.nq <- plm(lpatient_years ~ lstations + lstaff, index = c('cms_code'), model = c('within'), data = data.ols)
  # beta.fe.nq <- coef(fe.reg.nq)
  # estimate.old <- list(beta.ols, beta.ols.nq, beta.fe, beta.fe.nq)
}