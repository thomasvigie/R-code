# This function implements the P2SLS ("Polynomial two stage least squares") estimator. It contains the following functions:  
#    name_k: a function that returns a character in the form "char_k" where char and k are provided by the user
#    dummysator: a function that creates a vector of 1 and 0 corresponding to the a character in the "data" vector
#    dum_data: a function that creates all the dummies required in a dataset and adds them to the dataset
#    iv_k: a function that uses the ivreg function from AER package to compute the 2SLS estimates of X on Y using instruments Z. Polynomials of Z are created (to the order k), 
#          and if Z has dimension > 1 (more than one instrument), the polynomials of Z are orthogonalized using the Graham-Schmidt algorithm from the pracma package
#    pluginMSE: a function that computes the plug-in MSE of the two stage least squares estimates using the plug-in formula of Vigié (2021): "Improving 2SLS: polynomial-augmented 2SLS"
#    P2SLS: a function that finds the best polynomial combination of the instruments to include in the first stage for the estimation of slope parameters in a linear model using the two stage least squares estimator
#           Y is the dependent variable vector, X is the matrix of endogenous regressors, W is the matrix of exogenous regressors, K is the maximal order for polynomials of the instrument matrix Z
#           An example is provided at the end

name_k <- function(k, char)    # gives names according to some character and the numeric variable k
{
  paste(char, k, sep = "_")
}

# Dummysator creates a vector of 1 and 0 corresponding to the a character in the "data" vector
dummysator <- function(data, chr)
{
  res <- ifelse(data == chr, 1, 0)  
  return(res)
}

# A function that creates all the dummies required in a dataset and adds them to the dataset

dum_data <- function(data, var, chr, new_names)
{
  dum <- map(chr, dummysator, data = var)
  dum <- data.frame(reduce(dum, cbind))
  names(dum) <- new_names
  data <- cbind(data, dum)
  return(data)
}

iv_k <- function(Y, X, Z, W, k, weights = NULL)
{
  if(missing(W))
  {
    if (ncol(as.matrix(Z)) > 1) 
    { 
      psiz <- pracma::gramSchmidt(poly(Z, k, raw = F, simple = TRUE))$Q 
    }
    else 
    {
      psiz <-  poly(Z, k, raw = F, simple = TRUE)  
    }
    
    b_2sls <- AER::ivreg( Y ~ X | psiz, weights = weights)$coefficients["x"]  
  }
  else
  {
    if (ncol(as.matrix(Z)) > 1) 
    { 
      psiz <- pracma::gramSchmidt(poly(Z, k, raw = F, simple = TRUE))$Q 
    }
    else 
    {
      psiz <-  poly(Z, k, raw = F, simple = TRUE)  
    }
    b_2sls <- AER::ivreg( Y ~ X + W| psiz + W, weights = weights)$coefficients["X"]  # We use orthogonal polynomials
  }
  return(b_2sls)
}

pluginMSE <- function(Y, X, Z, k, W, weights = NULL, raw = F, int = c(TRUE, FALSE))
{
  if (missing(W))
  {
    if (ncol(as.matrix(Z)) > 1) 
    { 
      psiz <- pracma::gramSchmidt(poly(Z, k, raw = F, simple = TRUE))$Q 
    }
    else 
    { 
      psiz <- poly(Z, k, raw = raw, simple = TRUE)
    }
    nvar <-  ncol(as.matrix(X))
    psiz <- cbind(rep(1, nrow(as.matrix(Z))), psiz)
    PZ <- psiz %*% solve(t(psiz) %*% psiz) %*% t(psiz)
    xfit <- lm(X ~ psiz, weights = weights)$fitted
    vhat <- X - xfit
    
    muhat <- crossprod(xfit, xfit) / var(vhat)
    uhat <- AER::ivreg( Y ~ X | psiz, weights = weights)$residuals
    varbeta <- (((sum(uhat ^ 2))/(nrow(as.matrix(Y)) - nvar)) * solve(t(X) %*% PZ %*% X))[1, 1]
    
    if (int == TRUE)
    {
      X <- cbind(c(1), X) 
      nvar <- nvar + 1 # take the intercept into account
      varbeta <- (((sum(uhat ^ 2))/(nrow(as.matrix(Y)) - nvar)) * solve(t(X) %*% PZ %*% X))[2, 2] # the second diagonal element corresponds to x
    }
  }
  else 
  {
    if (int == TRUE)
    {
    Wi <- cbind(c(1), W)  
    }
    nvar <- ncol(as.matrix(Wi)) + ncol(as.matrix(X))
    MW <- diag(nrow(as.matrix(Wi))) - Wi %*% solve( t(Wi) %*% Wi ) %*% t(Wi)
   
    if (ncol(as.matrix(Z)) > 1) 
    { 
      psiz <- pracma::gramSchmidt(poly(Z, k, raw = F, simple = TRUE))$Q 
    }
    else 
    {
      psiz <-  poly(Z, k, raw = raw, simple = TRUE)  
    }
    
    psiW <- cbind(c(1), psiz, W)
    PZ <- psiW %*% solve(t(psiW) %*% psiW) %*% t(psiW)
    xfit <- lm(X ~ psiz + W, weights = weights)$fitted
    vhat <- X - xfit
  
    muhat <- ( t(xfit)%*%MW%*%xfit ) / ((sum(vhat^2))/( nrow(as.matrix(Y)) - ncol(psiW) ) )
    uhat <- AER::ivreg( Y ~ X + W| psiz + W, weigths = weights)$residuals
    
    varbeta <- (((sum(uhat ^ 2))/(nrow(as.matrix(Y)) - nvar)) * solve(t(X) %*% PZ %*% X))[1, 1]
    
    if (int == TRUE)
    {
      X <- cbind(c(1), X, W) 
      nvar <- nvar + 1 # take the intercept into account
      varbeta <- (((sum(uhat ^ 2))/(nrow(as.matrix(Y)) - nvar)) * solve(t(X) %*% PZ %*% X))[2, 2]
    }
  }

  k <- ncol(psiz)
  BWbias <- (cov(uhat, vhat) / var(vhat)) * (k / (muhat + k) - 2 * muhat ^ 2 / ((muhat + k) ^ 3) )
  MSE <- BWbias ^ 2 + varbeta 
  
  res <- data.frame(bias = abs(BWbias),
                    variance = varbeta,
                    MSE = MSE)
  return(res)
}

# Example
# n <- 500
# u <- rnorm(n)
# v <- rnorm(n)
# u <- v + rnorm(n) + rnorm(n)
# Z <- rnorm(n)
# X <- Z + Z^2+ v
# W <- rnorm(n)
# Y <- 1 + X + 2*W  + u
# pluginMSE(Y = y, X = x, Z = z, W = w, k = 2, weights = NULL, int = TRUE, raw = FALSE)

# The P2SLS estimator ------------------------------------ 

P2SLS <- function(Y, X, Z, K, W, option = c("MSE", "Bias", "Variance"), weights = NULL, raw = c(FALSE, TRUE), int = c(TRUE, FALSE))
{
  options(warn=-1)
  library(tidyverse)
  if(missing(W))
  {
    criterion <- 1:K %>% purrr::map_df(pluginMSE, Y = Y, X = X, Z = Z, weights = weights, raw = raw, int = int) 
    opti_k <- which.min(criterion[, option])
    
  if (ncol(as.matrix(Z)) > 1) 
  { 
    psi_z <- pracma::gramSchmidt(poly(Z, opti_k, raw = raw, simple = TRUE))$Q  # if one uses multiple instruments, I use teh Gram-Schmidt orthogonalization algorithm
  }
  else 
  {
    psi_z <- poly(Z, opti_k, raw = raw, simple = TRUE)
  }
    if (int == TRUE)
    {
      res <- AER::ivreg(Y ~  X  | psi_z )
    }
    else
    {
      res <- AER::ivreg(Y ~ -1 + X  | psi_z )
    }
    }
  else
  {
  criterion <- 1:K %>% purrr::map_df(pluginMSE, Y = Y, X = X, Z = Z, W = W, weights = weights, raw = raw, int = int) 
  opti_k <- which.min(criterion[, option])
  
   if (ncol(as.matrix(Z)) > 1) 
   { 
    psi_z <- pracma::gramSchmidt(poly(Z, opti_k, raw = raw, simple = TRUE))$Q  
   }
   else 
   {
    psi_z <- poly(Z, opti_k, raw = raw, simple = TRUE) 
   }
   if (int == TRUE)
   {
   res <- AER::ivreg(Y ~  X + W | psi_z + W)
   }
   else
   {
   res <- AER::ivreg(Y ~ -1 + X + W | psi_z + W)
   }
  }

  return(list(res = summary(res), 
              opti_k = opti_k, 
              criterion = criterion, 
              option = option) )
}

# Example
# n <- 500
# u <- rnorm(n)
# v <- rnorm(n)
# u <- v + rnorm(n) + rnorm(n)
# Z <- rnorm(n)
# X <- Z + Z^2+ v
# W <- rnorm(n)
# Y <- 1 + X + 2*W  + u
# P2SLS(Y = Y, X = X, Z = Z, W = W, K = 5, weights = NULL, option = "MSE", int = TRUE)
