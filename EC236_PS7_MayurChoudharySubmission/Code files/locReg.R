locReg <- function(y,X,X.shift){
  
  type.unique <- unique(X.shift)
  n.type <- nrow(type.unique)
  X.k <- ncol(X) ## dimensionality of X 
  
  ### initialize output vectors 
  y.hat <- matrix(0,nrow(y),1)
  sing.obs <- matrix(0,nrow(y),1)
  
    for (t in 1:n.type){
        if (ncol(X.shift) == 5){
          type.match <- as.matrix(X.shift[,1] == type.unique[t,1] & X.shift[,2] == type.unique[t,2] & X.shift[,3] == type.unique[t,3] & X.shift[,4] == type.unique[t,4] & X.shift[,5] == type.unique[t,5])
          y.type <- as.matrix(y[type.match])
          X.type <- as.matrix(X[type.match,])
        } else if (ncol(X.shift) == 4){
          type.match <- as.matrix(X.shift[,1] == type.unique[t,1] & X.shift[,2] == type.unique[t,2] & X.shift[,3] == type.unique[t,3] & X.shift[,4] == type.unique[t,4]) 
          y.type <- as.matrix(y[type.match])
          X.type <- as.matrix(X[type.match,])
        }
        
        n.local <- nrow(y.type)
        X.mean <- t(apply(X.type,2,mean)*matrix(1,ncol(X.type),nrow(X.type)))
        X.sd <-  t(apply(X.type,2,sd)*matrix(1,ncol(X.type),nrow(X.type)))
        ### normalize X.type to create the weight matrix 
        X.norm <- (X.type - X.mean)/X.sd
        ### initialize bandwidth 
        h = 1.5 * nrow(y.type)^(-1/(X.k +4))
        
        ### initialize yhat.type and sing.obs.type 
        y.hat.type <- matrix(0,nrow(y.type),1)
        sing.obs.type <- matrix(0,nrow(y.type),1)
        
        for (i in 1:nrow(y.type)){
            X.diff <- X.norm -  do.call(rbind, replicate(nrow(X.norm), X.norm[i,], simplify = FALSE))  ## get the x - x_1 distance 
            w <- apply(dnorm(as.matrix(X.diff/h) ,1,1),1,prod)
            w.sort <- order(w)
            w.threshold <- w[w.sort[length(w.sort)]] ## ideally we could choose a threshold for local regressions I simply weight all observations. 
            w.norm <- w/w.threshold ## normalize w so that the highest weight is 1 
            
            W <- diag(w.norm)
            X.reg <- as.matrix(cbind(matrix(1,nrow(X.type)), X.type))
            y.reg <- as.matrix(y.type)
            
            if(abs(rcond(ginv(t(X.reg) %*% W %*% X.reg))) < 10^(-18)){
              sing.obs.type[i] <- 1
            }  
            
            y.hat.type[i,] = X.reg[i,] %*% (ginv(t(X.reg) %*% W %*% X.reg))%*%(t(X.reg)%*% W %*% y.reg)
               
            
        }
       
        y.hat[type.match] <- y.hat.type
        sing.obs[type.match] <- sing.obs.type
    }
  
    return(list(y.hat,sing.obs)) ## return predicted yhat and singular observations (where the localized X matrix is singular)
}