# This file contains functions that estimate parameters in moment restrictions through Generalized Empirical Likelihood methods. It include:
#   - The Empirical Likelihood (EL) estimator
#   - The Exponential Tilting (ET) estimaor
#   - The Hellinger Distance (HD) estimator from Kitamura, Otsu and Evdokimov (2013)
#   - The Exponential Tilting Empirical (ETEL) Likelihood estimator from Schennach (2007)
#   - The Exponential Tilting estimator Hellinger (ETHD) Distance estimator from Antoine and Dovonon (2018)
#   - The continuously updated estimator (CU)
# The function uses the global optimization algorithm "DEoptim" from Ardia et al. (2012)

# Details about the functions --------------------------

# The functions are:
#   - smoothG: in time series contexts, the moment conditions can be smoothed out to account for serial dependence. Same function as in the "gmm" package.
#   - inner: function that returns the value of the inner function, a function of the Lagrange multipliers lambdas and theta, the parameter vector of interest.
#   - innergradient: function that returns the gradient of the inner function for each type of likelihood. Used in "getlambdas".
#   - outer: returns the outer function for any value of theta, i.e. the "inner" function after it is minimized w.r.t. lambdas.
#   - getlambdas: for each type of likelihood function, computes the optimal Lagrange multipliers lambdas by solving the "inner" problem
#   - impliedprobs: computes the probabilities on each observation that are implied by the Lagrange multipliers found after solving the "inner" problem.
#   - likelihood: computes the likelihood function using the implied probabilities (via the function "impliedprobs") computed for each type of empirical likelihood estimator. 
#                 The implied probabilities are a function of the Lagrange multipliers, computed via "getlambdas".
#   - GEL: the final function that computes the estimates by minimizing the "outer" function. It uses the "DEoptim" package from David Ardia, Katharine Mullen, Brian Peterson, Joshua Ulrich, and Kris Boudt (https://cran.r-project.org/web/packages/DEoptim/index.html). 
#          It returns the estimates of theta, the value of the objective function, the corresponding Lagrange multipliers and the implied probabilities.
                 

list.of.packages <- c('DEoptim','sandwich','gmm')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library('DEoptim')   # algorithm used for the GMM and GEL estimation
library('sandwich')  # used to smooth the moment function in the GEl case (See Kitamura and Sutzer (1997), Smith (1997), anatolyev(2005) )
library('gmm')

smoothG <- function (dat, bw = bwNeweyWest, prewhite = 0, ar.method = "ols", 
                     weights = weightsAndrews, kernel = c("Bartlett", "Parzen", 
                                                          "Truncated", "Tukey-Hanning"), approx = c("AR(1)", "AR(1)"), 
                     tol = 1e-07) 
{
  kernel <- match.arg(kernel)
  approx <- match.arg(approx)
  n <- nrow(dat)
  if (is.function(weights)) {
    class(dat) <- "gmmFct"
    w <- weights(dat, bw = bw, kernel = kernel, prewhite = prewhite, 
                 ar.method = ar.method, tol = tol, verbose = FALSE, 
                 approx = approx)
  }
  else {
    w <- weights
  }
  if (is.numeric(w)) {
    rt <- length(w)
    if (rt >= 2) {
      w <- c(w[rt:2], w)
      w <- w/sum(w)
      w <- kernel(w[rt:length(w)])
    }
    else {
      w <- kernel(1)
    }
  }
  else {
    if (class(w) != "tskernel") 
      stop("Provided weights must be a numeric vector or an object of class 'tskernel'")
  }
  if (length(w$coef) > 1) 
    dat <- kernapply(dat, w)
  sx <- list(smoothx = dat, kern_weights = w, bw = bw)
  return(sx)
}

inner <- function(momentfun , lambda , theta , dat,type = c("ET","EL","HD","CU"),smooth = c(TRUE,FALSE) , delta)
  
{
  if(missing(smooth))
  { smooth = FALSE }
  
  m <- momentfun(theta,x = dat )
  if (smooth == TRUE )
  {
    m <- smoothG( m , bw = bwAndrews, prewhite = 0, ar.method = "ols", 
                  weights = weightsAndrews, kernel = c( "Truncated"),  
                  tol = 1e-07)$smoothx
  }
  
  options(warn=-1)
    n <- nrow(m)
    lambda <- matrix(lambda , ncol = 1)

    nu <- as.vector(m%*%lambda )
    if (type == "ET")
    {   
      h <- - exp(nu)
    } 
    if (type == "EL")
    {
      y <- 1 - nu
      h <- as.matrix(ifelse(y > delta , log(y) , log(delta)-1.5+2*y/delta-y^2/(2*delta^2) ))
    } 
    if (type == "HD")
    {
      h <- 2/(1-nu/2)
    } 
    if (type == "CU")
    {
      h <- -((1 + nu)^2)/2
    }
  return(-mean(h))
  
}


innergradient <- function(theta , lambda , momentfun ,dat, type = c("ET","EL","HD","CU"),smooth = c(TRUE,FALSE),delta )
{
  
  if(missing(smooth))
  { smooth = FALSE }
  
  m <- momentfun(theta,x = dat )
  if (smooth == TRUE )
  {
    m <- smoothG( m , bw = bwAndrews, prewhite = 0, ar.method = "ols",
                  weights = weightsAndrews, kernel = c( "Truncated"),
                  tol = 1e-07)$smoothx
  }
  n <- nrow(m)
  lambda <- matrix(lambda , ncol = 1)
  nmoments <- ncol(m)
  nu <-  as.vector(m%*%lambda)
    
    if (type == "EL")
    {
      y <- matrix( rep( 1 - nu , each = nmoments  ) , ncol = nmoments , byrow = TRUE ) # We repeat the column to match the dimensions of m 
      h <- as.matrix(ifelse(y>delta , -m/y , -(2/delta)*m + m*y/(delta^2) ) )
    }
    
    if (type == "ET")
    {
      h <- - m*exp(nu) 
    }
  
  if (type == "HD")
  {
    h <-  4*m/( (2 - nu)^2 )
  }
  if (type == "CU")
  {
    h <- -(1 + nu)*m
  }
  
  grad <- -apply(h,2,sum)/n
  return(grad)
}


outer <- function(theta,momentfun,dat,type = c("ET","EL","HD","CU","ETEL","ETHD") , smooth = c(TRUE,FALSE) , delta)
{ 
  if(missing(delta)) {delta<- 0.001}
  if(missing(smooth))
  { smooth = FALSE }
  
m <- momentfun( theta, dat )
  
  if (smooth == TRUE )
  {
    m <- smoothG( m , bw = bwAndrews, prewhite = 0, ar.method = "ols", 
                        weights = weightsAndrews, kernel = c( "Truncated"),  
                        tol = 1e-07)$smoothx
  }

  nmoments <- ncol(m) 
  n <- nrow(m)
  # The problem is easy to solve in lambda, so we use a basic algorithm
 
  if (type == "ETEL")
  {
    innersol <- optim(matrix(0,nmoments,1), inner , gr = innergradient , theta = theta ,method = "BFGS" ,delta = delta , momentfun = momentfun , type = "ET" , dat = dat , smooth = smooth )
    m <- momentfun(theta,x = dat )
    nu <- as.vector(m%*%innersol$par )

    w <- -exp(nu)
    w <- w/sum(w)
    w <- ifelse(is.nan(w)==TRUE,1/(n^2),w)
    # if ( any(w==0)) 
    #   {
    #   w <- 1/n
    # }
    obj <- - mean(log(n*w))
    
  }
  else if (type == "ETHD")
  {
    innersol <- optim(matrix(0,nmoments,1), inner , gr = innergradient , theta = theta ,method = "BFGS" ,delta = delta , momentfun = momentfun , type = "ET" , dat = dat , smooth = smooth )
    m <- momentfun(theta,x = dat )
    nu <- as.vector(m%*%innersol$par )
    w <- -exp(nu)

    w <- w/sum(w)
    w <- ifelse(is.nan(w)==TRUE,1/(n^2),w)
    obj <- -4*(sqrt(n*w) - 1 )   # HD likelihood function
    
    obj <- mean(obj)
  }
  else
   {
      innersol <- optim(matrix(0,nmoments,1), inner , gr = innergradient , theta = theta ,method = "BFGS" ,delta = delta , momentfun = momentfun , type = type , dat = dat , smooth = smooth )
      obj <- -innersol$val 
     }   # The inner problem returns the solutions to the first minimization 

  return(obj)
  }

getlambdas <- function(theta,momentfun,dat,type = c("ET","EL","HD","CU"),smooth = c(TRUE,FALSE), delta)
{    
  
  if(missing(smooth))
  { smooth = FALSE }
  
  m <- momentfun( theta, dat )
  
  if (smooth == TRUE )
  {
    m <- smoothG( m , bw = bwAndrews, prewhite = 1, ar.method = "ols", 
                  weights = weightsAndrews, kernel = c( "Truncated"), 
                  tol = 1e-07)$smoothx
  }
  
  nmoments <- ncol(m) 
  n <- nrow(m)
  
  # The problem is easy to solve in lambda, so we use a basic algorithm
  innersol <- optim(matrix(0,nmoments,1), inner , gr = innergradient, theta = theta ,method = "BFGS" ,delta = delta , momentfun = momentfun , type = type , dat = dat , smooth = smooth )
  return(innersol$par )   # The inner problem returns the solutions to the first minimization 
}

impliedprobs <- function(theta, lambda, momentfun, dat, type = c("ET","EL","HD","CU"), smooth = c(TRUE,FALSE))
{ 
  
  if(missing(smooth))
  { smooth = FALSE }
  
  m <- momentfun( theta, dat )
  
  if (smooth == TRUE )
  {
    m <- smoothG( m , bw = bwAndrews, prewhite = 1, ar.method = "ols", 
                  weights = weightsAndrews, kernel = c( "Truncated"),  
                  tol = 1e-07)$smoothx
  }
  
  nmoments <- ncol(m) 
  n <- nrow(m)
  
  nu <- m%*%lambda
  if (type == "EL")
  {
  rho1 <- 1/(1-nu)
  pi <- rho1/sum(rho1)
  }
  if (type == "ET")
  {
  pi <- exp(nu)/(sum(exp(as.vector(nu))))
  }
  if (type == "HD")
  {
    rho1 <- 1/(1-nu)^2
    pi <-  rho1/sum(rho1)
  }
  if (type == "CU")
  {
    rho1 <- -(1 + nu)
    pi <- rho1/sum(rho1) 
  }
  return(pi)
}

likelihood <- function(theta, momentfun, dat, type=c("EL","ET","ETEL","HD","CU","ETHD"), smooth = c(TRUE,FALSE), delta)
{

  
  if (type == "EL")
  {
    obj <- mean(log(impliedprobs(theta,getlambdas(theta,momentfun=momentfun,dat=dat,type = "EL",smooth = smooth , delta),momentfun=momentfun,dat=dat,type = type , smooth = smooth ))*nrow(dat))
  }
  if (type == "ET")
  {
    obj <- - mean( nrow(dat)*impliedprobs(theta = theta, lambda = getlambdas(theta, momentfun = momentfun, dat = dat, type = "ET",smooth = smooth,delta),momentfun=momentfun,dat=dat,type = type , smooth = smooth)*log(nrow(dat)*impliedprobs(theta,getlambdas(theta,momentfun=momentfun,dat=dat,type = "ET",smooth = smooth,delta),momentfun=momentfun,dat=dat,type = type , smooth = smooth)) )
  }
  if (type == "HD")   # it corresponds to the likelihood of the HD criterion function, but the ET weights are used 
  {
    obj <- 4*mean(sqrt( (nrow(dat)* (impliedprobs(theta,getlambdas(theta,momentfun=momentfun,dat=dat,type = "HD",smooth = smooth,delta),momentfun=momentfun,dat=dat,type = "HD",smooth = smooth)) )) - 1 )
  }
  if (type == "CU")
  {
    obj <- - mean((impliedprobs(theta,getlambdas(theta,momentfun=momentfun,dat=dat,type = "CU",smooth = smooth , delta),momentfun=momentfun,dat=dat,type = type , smooth = smooth )*nrow(dat))^2)
  }
  if (type == "ETEL")   # it corresponds to the likelihood of the EL criterion function, but the ET weights are used 
  {
    obj <- mean(log(impliedprobs(theta,getlambdas(theta,momentfun=momentfun,dat=dat,type = "ET",smooth = smooth,delta),momentfun=momentfun,dat=dat,type = "ET",smooth = smooth))*nrow(dat))
  }
  if (type == "ETHD")   # it corresponds to the likelihood of the HD criterion function, but the ET weights are used 
  {
    obj <- 4*mean(sqrt((nrow(dat)* (impliedprobs(theta,getlambdas(theta,momentfun=momentfun,dat=dat,type = "ET",smooth = smooth,delta),momentfun=momentfun,dat=dat,type = "ET",smooth = smooth)) )) - 1 )
  }
  
  return(-obj)
}


GEL <- function(x0 , momentfun, dat, type = c("EL" , "ET" , "ETEL" , "HD" , "CU" , "ETHD") , smooth = c( TRUE , FALSE ) , lb  , ub , delta , trace , cluster , maxiter )
{
  if(missing(trace)) {trace <- FALSE }
  if(missing(cluster)){cluster <- NULL}
  if(missing(maxiter)){maxiter <- 500}
  if(missing(delta)){delta <- 0.01}
  options(warn=-1)
  
  n <- nrow(dat)
  nmoments <- ncol(momentfun( x0, dat ) ) 

outersol <- DEoptim(outer , lower = lb , upper = ub , momentfun = momentfun , dat = dat , delta = delta , type = type ,smooth = smooth, control = DEoptim.control(itermax = maxiter , trace = trace , strategy = 2 ,CR = 1 , NP = 20 , cluster = cluster))
thetahat <- outersol$optim$bestmem
Qhat <- outersol$optim$bestval

if (type == "ETEL" | type == "ETHD")
{
  lambdahat <- getlambdas(theta=thetahat,momentfun=momentfun,dat=dat,type="ET", smooth = smooth , delta = delta)
  pihat <- impliedprobs( thetahat , lambdahat , momentfun , dat = dat , type = "ET" , smooth = smooth )
}

else {
lambdahat <- getlambdas(theta=thetahat,momentfun=momentfun,dat=dat,type=type,delta = delta)
pihat <- impliedprobs( thetahat , lambdahat , momentfun , dat = dat , type = type )
}

return(list( estimates = thetahat , objective = Qhat , multipliers = lambdahat , probabilities = pihat ) )
}


