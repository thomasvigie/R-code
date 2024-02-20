# This function estimates the parameters from moment conditions using the 2-step method or the CUE ("continusously updated") method  

# Details ----------------------------
# This script contains the following functions:
#   - eq: functions that checks if two matrices x and y have the same elements everywhere
#   - gbar: for some value of the parameters "theta" and a moment function "momentfun", 
#           computes the average of each function in "momentfun" and returns a column vector of the size of the number of functions in "momentfun".
#   - Q: for some value theta, the moment function "momentfun", some data and the weighting matrix W, computes the GMM objective function.
#   - Qcue: for some value theta, the moment function "momentfun", and some data , computes the continuously updted GMM (CU-GMM) objective function. 
#           If it is a time series environment, the long run variance of the moments  can be used for the weighting matrix.
#   - jacobar: compute the average of the jacobian matrix of "momentfun". No numerical derivatives, "jacobian" must be provided.
#   - secondmom: computes the second moments of a matrix
#   - GMM: minimizes "Q" (if the method is "twostep", "onestep" or "iterative") or "Qcue" (if the method is "cue") using the "DEoptim" package from Ardia et al. (see https://cran.r-project.org/web/packages/DEoptim/DEoptim.pdf).
#          Returns the estimates of, the value of the objective function, the estimate of the variance-covariance matrix of the estimates, Hansen's J-test and whether a generalized inverse had to be used for the weighting matrix.  

  
list.of.packages <- c('DEoptim', 'sandwich', 'MASS','pracma')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library('DEoptim')   # algorithm used for the GMM estimation
library('sandwich')
library('MASS')
library('pracma')
source('lvar.R')



eq <- function(x,y,tol)   # this function works the best I believe, when the matrices contain numbers. Does not work for logical things or characters yet
{ 
  if(missing(tol)){tol <- 0.00001}
  x <- as.matrix(x)
  y <- as.matrix(y)
if( (length(x)!=length(y))==TRUE ) {stop("Dimensions don't match mate !")}
  iter <- 1
  while(any(abs(x[iter] - y[iter]) <= tol) && iter <= nrow(x)*ncol(x) )
  {
    iter <- iter +1          
  }
  if(iter < nrow(x)*ncol(x)) {out <- 'na'}
  else{out <- 'yes'}
  return(out)
}



gbar <- function (theta,momentfun,dat)    # returns the average value of ghat : gbar
{
  g <- apply( momentfun( theta = theta , x = dat ) , 2 , mean )     # mean of each column who represents a different moment
  
  return( data.matrix(g) )   # reports a column vector, with as many rows as moments
}


Q <- function (theta,momentfun,dat,W)    # Objective function
{
  Q <- crossprod( gbar( theta = theta , momentfun = momentfun , dat = dat ) , W )%*%( gbar( theta = theta , momentfun = momentfun , dat = dat) )
  return(Q)
}

Qcue <- function (theta, momentfun, dat, env = c("iid","ts"))    # Objective function
  
{
  if(missing(env)) {env <- "ts"}
  
  f1 <- momentfun(theta = theta , dat  )
  
  if(env=="iid"){    Wmat <-  cov( f1 )   + gbar( theta = theta , momentfun = momentfun , dat = dat )%*%t( gbar( theta = theta , momentfun = momentfun , dat = dat ) )       # we use a generalized inverse since there are singularity issues. The function cov comes the stats package
  
  }   else      {    Wmat <-  longrunvar(f1,center=TRUE,method='HACC_B')[[1]] }
 
  Wcue <- MASS::ginv( Wmat)

  Qcue <- t(gbar(theta,momentfun,dat)) %*% Wcue %*% (gbar(theta,momentfun,dat))
  
  return(Qcue)
}

jacobar <- function (jacobian , theta , dat)    # computes the average jacobian in a not too inefficient way
{
  J <- list()
  for (i in 1:nrow(data.matrix(dat)))
  {
    y <- as.matrix( as.matrix(dat)[i,] )
    J[[i]] <- jacobian( theta , x = t(y))
  }
  return((1/nrow(data.matrix(dat)))*Reduce('+',J))
}

secondmom <- function (momentfun , theta , dat)    # computes the second moments of a matrix

{
  J <- list()
  f <- momentfun(theta,dat)
  
  for (i in 1:nrow(data.matrix(dat)))
  {
    J[[i]] <- cbind(f[i,])%*%rbind(f[i,])
  }
  return((1/nrow(data.matrix(dat)))*Reduce('+',J))
}

GMM <- function( momentfun , x0  , jacobian , dat , method = c("twostep","cue","iterative", "onestep") , env = c("iid","ts") , optimethod = c('optim','DEoptim')  , lb , ub , maxiter , trace , cluster )
{
  
  library(pracma)
  
  options(warn=-1)
  if(missing(env)) {env <- "ts"}
  if(missing(cluster)) {cluster <- NULL}
  if(missing(maxiter)) {maxiter <- 500 }    
  if(missing(trace)) {trace <- FALSE }
  if(missing(lb) | missing(ub)  ){warnings('Upper or lower bound not supplied for optimethod DEoptim; optim is used instead')
    optimethod <- "optim"}
  if(missing(optimethod)){optimethod <- 'optim'}

  n <- nrow(data.matrix(dat) )
  
  nmoments <- ncol ( momentfun ( x0  , dat ) )     # number of moments used by evaluating the function at random points. We should find a better way to get that info
  
  
  if (method == "onestep")
  {
    if(optimethod=='optim') {
      step1 <- optim( par = x0 , fn = Q , gr = NULL , dat = dat , momentfun = momentfun ,  W = diag(nmoments) , method = 'Nelder-Mead', control = list(trace=trace) , hessian = FALSE)
      sol <- step1$par
      obj <- step1$value
    }
    
    
    else{step1 <- DEoptim( Q , lower = lb ,upper = ub , control = DEoptim.control(itermax = maxiter ,strategy = 2, trace = trace , cluster = cluster) , momentfun = momentfun , dat = dat , W = diag(nmoments)  )     # Differential evolution algorithm
    sol <- step1$optim$bestmem    # output when using DEoptim
    obj <- step1$optim$bestval    }# with DEoptim

  }
  # Step 1 : minimize gbar'gbar and get W1 with theta1
  if (method == "twostep")
  {
    if(optimethod=='optim') {
      step1 <- optim( par = x0 , fn = Q , gr = NULL , dat = dat , momentfun = momentfun  , W = diag(nmoments) , method = 'Nelder-Mead', control = list(trace=trace) , hessian = FALSE)
      theta1 <- step1$par
      obj1 <- step1$value
    }
    
  else {  
  step1 <- DEoptim(Q , lower = lb ,upper = ub , dat = dat , momentfun = momentfun , W = diag(nmoments) ,DEoptim.control(itermax = maxiter , trace = trace , strategy = 2, cluster = cluster) )     # Differential evolution algorithm
   
  theta1 <- step1$optim$bestmem    # output when using DEoptim
       }
    
  nparam <- nrow(data.matrix(theta1) )    # obtain the dimension of theta here
  f1 <- momentfun(theta1 , dat )
  
  if(env=="iid"){    Wmat <-  var( f1 )   + gbar(theta1 , momentfun ,dat)%*%t( gbar ( theta1 , momentfun , dat ) )       

  }   else      {    Wmat <-  longrunvar(f1,center=TRUE,method='HACC_B')[[1]] }
  
  Wmat <- matrix(as.numeric(Wmat),nmoments,nmoments)
  invindicator<-eq(solve(Wmat),MASS::ginv(Wmat))
  
  Wmat <- MASS::ginv(Wmat)

  # Step 2 : minimize gbar' %*% Wmat %*% gbar and get W2 with theta2
  if(optimethod=='optim') {
    step2 <- optim( par = x0 , fn = Q , gr = NULL , dat = dat , momentfun = momentfun , W = Wmat , method = 'Nelder-Mead', control = list(trace=trace) , hessian = FALSE)
    sol <- step2$par
    obj <- step2$value
  }
  else{
    step2 <- DEoptim(Q , lower = lb ,upper = ub , dat = dat , momentfun = momentfun , W = Wmat ,DEoptim.control(itermax = maxiter ,strategy = 2,trace = trace, cluster = cluster))     # Differential evolution algorithm
    sol <- step2$optim$bestmem    # when DEoptim is used
    obj <- step2$optim$bestval
  }
  
  f1 <- momentfun( sol , dat )
  if(env=="iid"){    Wmat <-  var( f1 )   + gbar(sol , momentfun ,dat)%*%t( gbar ( sol , momentfun , dat ) )       
  
  }   else      {    Wmat <-  longrunvar(f1,center=TRUE,method='HACC_B')[[1]] }
  
  
  Wmat <- matrix(as.numeric(Wmat),nmoments,nmoments)
  invindicator<-eq(solve(Wmat),MASS::ginv(Wmat))
  
  Wmat <- MASS::ginv(Wmat)
  
  }
  
  else if (method == "iterative")
  {
    if(optimethod=='optim') {
      step1 <- optim( par = x0 , fn = Q , gr = NULL , dat = dat , momentfun = momentfun ,  W = diag(nmoments) , method = 'Nelder-Mead', control = list(trace=trace) , hessian = FALSE)
      theta1 <- step1$par
      obj1 <- step1$value
                            }
    else{
    step1 <- DEoptim( Q , lower = lb , upper = ub , dat = dat , momentfun = momentfun , W = diag(nmoments) ,DEoptim.control(itermax = maxiter ,strategy = 2,trace = trace, cluster = cluster))     # Differential evolution algorithm
    
    theta1 <- step1$optim$bestmem
    obj1 <- step1$optim$bestval
        }
    nparam <- nrow(data.matrix(theta1) )    # obtain the dimension of theta here
    f1 <- momentfun (theta1 , dat )
    
    theta0 <- x0
    iter <- 0
    dev <- 1
    itermax <- 1000
    
    while (dev > 0.000001 && iter < itermax)
    {
      W1 <- MASS::ginv( cov(( f1 ) ) + gbar(theta1 , momentfun ,dat)%*%t( gbar ( theta1 , momentfun , dat ) )     , tol = sqrt(.Machine$double.eps) )    # we use a generalized inverse since there are singularity issues. The function cov comes the stats package
      
      # Step 2 : minimize gbar' %*% W1 %*% gbar and get W2 with theta2
      if(optimethod=='optim') {
        step2 <- optim( par = x0 , fn = Q , gr = NULL , dat = dat , momentfun = momentfun ,  W = W1 , method = 'Nelder-Mead', control = list(trace=trace) , hessian = FALSE)
        thetanew <- step2$par
        objnew <- step2$value
                              }
      else          {
        step2 <- DEoptim( Q , lower = lb , upper = ub , dat = dat , momentfun = momentfun , W = W1 , DEoptim.control(itermax = maxiter ,strategy = 2,trace = trace, cluster = cluster))     # Differential evolution algorithm
        thetanew <- step2$optim$bestmem
        objnew <- step2$optim$bestval
                    }
    dev <- abs(obj1 - objnew)
    fnew <- momentfun (thetanew , dat )      # update if the condition in the while statement is not violated 
    theta1 <- thetanew
    obj1 <- objnew
    f1 <- fnew
    iter <- iter + 1
    }
    
  sol <- thetanew
  obj <- objnew
  Wmat <- W1
  rm(f1)
  rm(fnew)
  rm(dev)
    
  }
  
  else if (method == "cue")
  {
    if(optimethod=='optim') {
      cuesol <- optim( par = x0 , fn = Qcue , gr = NULL , dat = dat , momentfun = momentfun , env = env , method = 'Nelder-Mead', control = list(trace=trace) , hessian = FALSE)
      sol <- cuesol$par
      obj <- cuesol$value
                            }
   cuesol <- DEoptim(Qcue , lower = lb ,upper = ub , dat = dat , momentfun = momentfun , env = env  , control = DEoptim.control(itermax = maxiter, strategy = 2 , trace = trace, cluster = cluster))     # Differential evolution algorithm
  
    sol <- cuesol$optim$bestmem
    obj <- cuesol$optim$bestval
    f1 <- momentfun(sol , dat )
    
    if(env=="iid"){    Wmat <-  cov( ( f1 ) ) + gbar(sol , momentfun ,dat)%*%t( gbar ( sol , momentfun , dat ) )       
    }   else{          Wmat <-  longrunvar(f1,center=TRUE,method='HACC_B')[[1]]  }
    
    
    invindicator<-eq(solve(Wmat),MASS::ginv(Wmat))
    
    Wmat <- ginv(Wmat)
  }
  
  f1 <- momentfun(sol , dat )
  if(env=="iid"){    Wmat <-  cov( ( f1 ) ) #+ gbar(sol , momentfun ,dat)%*%t( gbar ( sol , momentfun , dat ) )       
  }   else{          Wmat <-  longrunvar(f1,center=TRUE,method='HACC_B')[[1]]  }
  
  Wmat <- MASS::ginv(Wmat)
  
  Ejac <- jacobar(jacobian , sol , dat)    # Expectd Jacobian, nothing else...
  VARHAT <- solve(crossprod(Ejac,Wmat)%*%Ejac)/n
  rm(Ejac)
  
  
  if (method == "iterative")
  {
    Jobj <- Q(sol,momentfun,dat,W=Wmat)
    Jtest <- n*Jobj
    output <- list(sol,obj,Jtest,VARHAT,iter )
    names(output) <- c('Estimates','Objective function','Jtest','Estimated_variance','number of iterations')
    return (output)
  }
  else if (method == "onestep")
  {
    output <- list( Estimates = sol , Objective = obj )
    names(output) <- c("Estimates", "Objective function" )
    
    return(output)
  }
  
   else 
  {
    Jobj <- Q(sol,momentfun,dat,W=Wmat)
    Jtest <- n*Jobj
    output <- list(sol, obj, VARHAT, Jtest, invindicator)
    names(output) <- c("Estimates", "Objective_function", "Estimated_variance", "Jtest", "regular_inverse_used")
    
    return (output )
  }
}
