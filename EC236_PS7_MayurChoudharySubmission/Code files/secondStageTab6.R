secondStageTab6 <- function(alpha.q, data.subset){
  
  ### Second Stage to estimate beta and OLS/FE parameters [Create a function for this later] ###
  ### Subset data to get rid of observations for which lag is not available 
  
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
  beta.tab6.col1 <- as.matrix(second.stage.reg$par) 
  
  return(beta.tab6.col1)

}