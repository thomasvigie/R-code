# this function implements the weighted minimum distance estimator proposed by Antoine and Lavergne (2014)
# Works with one instrument only for now


rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
rep.row<-function(x,n)    # functions to repeat a vector in the repmat way
{
  matrix(rep(x,each=n),nrow=n)
}

wmd <- function( x , y , z , intercept = c("TRUE", "FALSE"))
{
 if(missing(intercept)){
   intercept <- TRUE
 }
  n <- nrow(as.matrix(x))
  
  if(intercept == TRUE){
    e <- matrix(1, nrow = n, ncol = 1 )
    Y1star <- cbind( e, x )
  }
  else {
    Y1star <- x
  }
  
  Ystar <- cbind( y , Y1star )
  
  if (is.vector(z) == FALSE)
  {
  z<- z/rep.row((sqrt(diag(var(z)))), n )    # Doing some rescaling for the weights
  Ktilde <- matrix(1,nrow = nrow(z),ncol = n)
  for(i in 1:ncol(z))
  {
    K <- dnorm(rep.col(z[,i],nrow(z))-rep.row(t(z[,i]) , n))
    K <- K - diag(diag(K))
    Ktilde <- Ktilde*K
  }
  } else {  
    z <- z/sd(z)
    Ktilde <- dnorm(rep.col(z,n)-rep.row(z,n))
    }
  
  eig <- eigen(solve(t(Ystar)%*%Ystar)%*%(t(Ystar)%*%Ktilde%*%Ystar))
  lambdatilde <- min(eig$values)
  
  lambda_WMDF <- (lambdatilde - (1 - lambdatilde)/n)/(1 - (1 - lambdatilde)/n)
  
  WMD <- solve(t(Y1star)%*%(Ktilde-lambdatilde*diag(n))%*%Y1star)%*%(t(Y1star)%*%(Ktilde-lambdatilde*diag(n))%*%y)
  WMDF <- solve(t(Y1star)%*%(Ktilde-lambda_WMDF*diag(n))%*%Y1star)%*%(t(Y1star)%*%(Ktilde-lambda_WMDF*diag(n))%*%y)

  res <- list(WMD, WMDF)
  names(res) <- c("WMD", "WMDF")
  return(res)
}



# Example ------------------
# n <- 1000
# z1 <- rnorm( n )
# z2 <- rnorm( n , 2 )
# z <- cbind(z1, z2)
# 
# u <- rnorm( n )
# x <- rnorm( n ) + u^2 + .1*z1^2 + .1*exp(z2)
# y <-  1 + 2*x + u
# 
# wmd( x , y , z, intercept = TRUE )





