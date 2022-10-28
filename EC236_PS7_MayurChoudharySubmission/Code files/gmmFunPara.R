gmmFunPara <- function(beta, tf.t, factors.t, phi.lag, factors.lag, W, order, type){
  
  omega.t <- as.matrix(tf.t) - as.matrix(factors.t) %*% beta  ## productivity = tf - beta*factors (this is the LHS = omega_t)
  omega.lag <- as.matrix(phi.lag) - as.matrix(factors.lag) %*% beta ## productivity = Phi - beta*factors (this is the argument of g(omega_{t-1}))
  
  if(is.null(type)){
    type <- matrix(1,nrow(tf.t),1) ## main specification g(omega) function same for all
  }
  
  
  
  ## initialize eta (the residual) and the coefficients gamma of the g() polynomial function 
  eta <- matrix(0,nrow(factors.t),ncol = 1)
  gamma <- matrix(0,order + 1, 1)
  order.vector <- t(matrix(seq(0,order),(order+1),nrow(omega.lag)))
  
  ## compute the g function for each type (in the main specification only one type)
  
 
    g <- (matrix(omega.lag, nrow(omega.lag),(order+1)))^order.vector
    g.dummy <- as.matrix(cbind(type,g))
    gamma <- ginv(t(g.dummy) %*% g.dummy) %*% (t(g.dummy) %*% omega.t)
    eta <- omega.t - g.dummy %*% gamma
    
 
  gmm.mom <- t(factors.t) %*% eta 
  gmm.obj <- t(as.matrix(gmm.mom)) %*% W %*% as.matrix(gmm.mom)
  return(gmm.obj)
}