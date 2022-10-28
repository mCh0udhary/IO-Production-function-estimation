modelTab8 <- function(data.subset){
  
  ## change command names for MASS and dplyr 
  dselect <- dplyr::select
  dfilter <- dplyr::filter
  
  ### Set up for first stage ### 
  y <- data.subset %>% dselect(lpatient_years)
  X <- data.subset %>% dselect(lhires,lstations, lstaff, inspection_rate, time_since_survey)
  X.shift <- data.subset %>% dselect(hirepos, for_profit, compLevel, davita, fresenius)
  
  q.input <- data.subset %>% dselect(q.res,q.k,q.l)
  q.input.iv <- data.subset %>% dselect(q.iv.res,q.k.iv,q.l.iv)
  
  first.stage.reg <- firstStage(y,X,X.shift,as.matrix(q.input), as.matrix(q.input.iv))
  alpha.q <- first.stage.reg[[1]]
  Phi.t <- first.stage.reg[[2]]
  
  
  ### Add the recovered Phi.t and Phi.t.nq to the original database
  data.subset <- cbind(data.subset, Phi.t)
  data.subset <- data.subset %>% filter(!(is.na(Phi.t)))
  
  
  ### Second Stage to estimate beta and OLS/FE parameters  ###
 
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
  data.subset.lag <- inner_join(data.subset, cur.id, by = 'for_id')
  data.subset.lag <- unique(data.subset.lag) ## get rid of duplicates due to inner join 
  
  
  ### Will need to modify this for table 8 
  q.input.t <- data.subset.t %>% dselect(q.res,q.k,q.l)
  data.subset.t <- data.subset.t %>% mutate(tf.t = lpatient_years - as.matrix(q.input.t) %*% alpha.q)  ### The transformation function (LHS = y + alpha*q) 
  tf.t <- as.matrix(data.subset.t$tf.t)
  
  ### Get lagged and time t data for input to the GMM function 
  factors.t  <- data.subset.t %>% dselect(lstations, lstaff)
  factors.lag <- data.subset.lag %>% dselect(lstations,lstaff)
  phi.lag <- data.subset.lag %>% dselect(Phi.t)
  
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
  beta.tab8 <- as.matrix(second.stage.reg$par) 
  
  
  
  #### OLS and FE regressions ##### 
  data.ols <- data.subset.t %>% dselect(lpatient_years,lstations,lstaff,q.res,q.k,q.l,cms_code)
  
  ols.reg <- lm(lpatient_years ~  q.res + q.k + q.l + lstations + lstaff, data = data.ols)
  beta.ols <- as.matrix(coef(ols.reg))
  
  
  fe.reg <- plm(lpatient_years ~  q.res + q.k + q.l  + lstations + lstaff, index = c('cms_code'), model = c('within'), data = data.ols)
  beta.fe <- as.matrix(coef(fe.reg))
  
  tab8.output <- rbind(alpha.q, beta.tab8, as.matrix(beta.ols[2:6,]), beta.fe)
  
  return(tab8.output)
  
}