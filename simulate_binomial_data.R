################################################################################
############################### Simulation study ###############################
################################################################################

# This is an adaptation of the script by Zhi Zhao (zhi.zhao@medisin.uio.no) 
# functions to simulate two heterogeneous data sources with bivariate response
# used in Zhao and Zucknick (2020).

library(Matrix)
library(MASS)
library(mvtnorm)

covariance_matrix <- function(rho, p, b, num_nonpen) {
  # Generate covariance matrix
  
  Ap1 <- matrix(rep(rho, (p[1] / b) ^ 2), nrow = p[1] / b)
  diag(Ap1) <- rep(1, p[1] / b)
  
  Ap2 <- matrix(rep(rho, (p[2] / b) ^ 2), nrow = p[2] / b)
  diag(Ap2) <- rep(1, p[2] / b)
  
  Bp12 <- matrix(rep(rho, p[1] / b * p[2] / b), nrow = p[1] / b)
  Bp21 <- matrix(rep(rho, p[1] / b * p[2] / b), nrow = p[2] / b)
  
  Xsigma1 <- Ap1
  Xsigma2 <- Ap2
  Xsigma12 <- Bp12
  Xsigma21 <- Bp21
  
  for (i in 2:b) {
    Xsigma1  <- bdiag(Xsigma1, Ap1)
    Xsigma12 <- bdiag(Xsigma12, Bp12)
    Xsigma2  <- bdiag(Xsigma2, Ap2)
    Xsigma21 <- bdiag(Xsigma21, Bp21)
  }
  Xsigma <-
    rbind(cbind(Xsigma1, Xsigma12), cbind(Xsigma21, Xsigma2))

  if(num_nonpen > 0){

    nonpenSigma <-
      matrix(rep(rho, num_nonpen^2),  nrow = num_nonpen)
    diag(nonpenSigma) <- rep(1, num_nonpen)

    Xsigma <- bdiag(nonpenSigma, Xsigma)
  }

  return(as.matrix(Xsigma))
}


simulation_data <-
  function(p = c(100, 100),
           p_relevant_set = c(10, 10),
           beta_set = c(0.2, 0.6),
           num_nonpen = 0,
           beta_nonpen = 0.7,
           n = 100,
           rho = .4,
           prop = .5,
           dicotomise_second_layer = FALSE,
           diagonal_covariance = FALSE) {

    b = 10
    
    if (is.na(p[2]) |
        is.na(p_relevant_set[2]) | is.na(beta_set[2])) {
      cat("Please specify bidimensional p, p_relevant_set, and beta_set
          vectors.\n")
    }
    
    Y <- rbinom(n, 1, prop)
    if (diagonal_covariance) {
      Xsigma <- diag(p[1] + p[2] + num_nonpen)
    } else{
      Xsigma <- covariance_matrix(rho, p, b, num_nonpen)
    }
    X <- rmvnorm(n,
                 mean = rep(0, num_nonpen + p[1] + p[2]),
                 sigma = as.matrix(Xsigma))
    
    X0 <- X[, 1:num_nonpen]
    X1 <- X[, (1 + num_nonpen):(p[1] + num_nonpen)]
    X2 <- X[, (p[1] + 1 + num_nonpen):(p[1] + p[2] + num_nonpen)]
    
    relevant_covariates <- c()
    
    if (p_relevant_set[1] > 0 & beta_set[1] != 0) {
      
      ind.rel <- sample(ncol(X1), p_relevant_set[1])
      
      relevant_covariates <- 1:p_relevant_set[1] + num_nonpen
      
      X1 <- cbind(X1[, ind.rel], X1[, -ind.rel])
      X1[Y == 1, 1:p_relevant_set[1]] <-
        X1[Y == 1, 1:p_relevant_set[1]] + beta_set[1]
    }
    
    if (p_relevant_set[2] > 0 & beta_set[2] != 0) {
      ind.rel <- sample(ncol(X2), p_relevant_set[2])
      
      relevant_covariates <- c(relevant_covariates,
                               1:p_relevant_set[2] + num_nonpen + p[1])
      
      X2 <- cbind(X2[, ind.rel], X2[, -ind.rel])
      X2[Y == 1, 1:p_relevant_set[2]] <-
        X2[Y == 1, 1:p_relevant_set[2]] + beta_set[2]
    }
    
    if (dicotomise_second_layer) {
      X2 <- data.matrix(X2 > 0) + 0 # Make this binary
    }

    X = cbind(X1, X2)
    
    if (num_nonpen > 0) {
      relevant_covariates <- c(relevant_covariates, 1:num_nonpen)
      
      X0[Y == 1,] <- X0[Y == 1,] + beta_nonpen
      X <- cbind(X0, X)
    }
    
    return(list(Y = Y, # Dependent variable
                X = X, # Independent variables
                p = p, # Number of variables in each layer
                relevant_covariates = sort(relevant_covariates)
                ))
}
