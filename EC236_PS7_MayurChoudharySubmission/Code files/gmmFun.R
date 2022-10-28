gmmFun <- function(beta, tf.t, factors.t, phi.lag, factors.lag, W, order, type){
  
  omega.t <- as.matrix(tf.t) - as.matrix(factors.t) %*% beta  ## productivity = tf - beta*factors (this is the LHS = omega_t)
  omega.lag <- as.matrix(phi.lag) - as.matrix(factors.lag) %*% beta ## productivity = Phi - beta*factors (this is the argument of g(omega_{t-1}))
  
  if(is.null(type)){
    type <- matrix(1,nrow(tf.t),1) ## main specification g(omega) function same for all
  }
  
  type.num <- max(type) ## number of types 
   
  ## initialize eta (the residual) and the coefficients gamma of the g() polynomial function 
  eta <- matrix(0,nrow(factors.t),ncol = 1)
  gamma <- matrix(0,order + 1, 1)
  order.vector <- t(matrix(seq(0,order),(order+1),nrow(omega.lag)))
  ## compute the g function for each type (in the main specification only one type)
  
  for (t in 1:type.num){
    omega.t.mat <- as.matrix(omega.t[type == t])
    g <- (matrix(omega.lag[type == t], length(omega.lag[type == t]),(order+1)))^order.vector[type == t,]
    gamma <- ginv(t(g) %*% g) %*% (t(g) %*%  omega.t.mat)
    eta[type == t] <- omega.t.mat - g %*% gamma
    
  }
  
  gmm.mom <- t(factors.t) %*% eta 
  gmm.obj <- t(as.matrix(gmm.mom)) %*% W %*% as.matrix(gmm.mom)
  return(gmm.obj)
}