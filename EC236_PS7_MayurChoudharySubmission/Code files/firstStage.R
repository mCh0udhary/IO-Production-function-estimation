firstStage <- function(y,X,X.shift,q,q.iv){
  
  y.hat <- locReg(y,X,X.shift)
  
  q.num <- ncol(q)
  q.hat <- matrix(0, nrow(q), ncol(q))
  q.iv.hat <- matrix(0, nrow(q), ncol(q)) 
  
  for (j in 1:q.num){
   q.temp.reg <-  locReg(as.matrix(q[,j]),X,X.shift)
   q.hat[,j] <- q.temp.reg[[1]]
   
   q.iv.temp.reg <- locReg(as.matrix(q.iv[,j]), X, X.shift)
   q.iv.hat[,j] <- q.iv.temp.reg[[1]]
  }
  
  y.reg <- y - y.hat[[1]] 
  q.reg <- q - q.hat
  q.iv.hat.reg <- q.iv - q.iv.hat
  
  y.final <- y.reg[y.hat[[2]] == 0] ## only take observations where the sing.obs == 0 (no singular observations)
  q.final <- q.reg[y.hat[[2]] == 0,]
  q.iv.final <- q.iv.hat.reg[y.hat[[2]] == 0,]
  
  alpha.q <- ginv(t(q.iv.final) %*% q.final) %*% (t(q.iv.final) %*% y.final)
  Phi <- y - alpha.q[1,1] * q 
  ## Use the local regression sieve to predict phi hat 
  
  Phi.hat <- locReg(Phi,X,X.shift)
  Phi.t <- Phi.hat[[1]] 
  Phi.t[Phi.hat[[2]] > 0] <- NaN
  
  
  return(list(alpha.q, Phi.t))

  }