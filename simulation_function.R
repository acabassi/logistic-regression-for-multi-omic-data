################################################################################
############################### Simulation study ###############################
################################################################################

# This is an adaptation of the script by Zhi Zhao (zhi.zhao@medisin.uio.no) used
# for the analysis of simulated data in the Zhao & Zucknick (2020) 

#'@param num iteration number / seed for data generation (?).
#'@param p_set vector containing the number of covariates in each layer.
#'@param method can be \code{elastic-net}, \code{sIPF-elastic-net},
#'\code{separate-elastic-net}, or \code{fixed-alpha}.
#'@param lambda optional user-supplied lambda sequence.
#'@param bounds bounds for the interval-searching parameters.
#'@param N number of start points.
#'@param min.iter minimum number of iterations after first N.
#'@param parallel If \code{TRUE}, use parallel foreach to fit each fold except
#'parallelizing each lambda for the tree-lasso methods. If \code{c(TRUE,TRUE)},
#'use parallel foreach to fit each fold and each lambda.
#'@param verbose print the middle search information, default is \code{TRUE}.
#'@param dicotomise_second_layer Boolean. If TRUE, the second layer has binary
#'values, otherwise it has continuous values. Default is FALSE.
#'@param diagonal_covariance Boolean. If TRUE, the covariance is a diagonal
#'matrix, otherwise it is the function specified in Chapter 2/Appendix A
#'(depending on whether or not there are non-penalised covariates) of my thesis.
#'@param fixed_alpha value of the parameter alpha for the fixed-alpha method. 
#'Default is 0.1.

simulation_function <- function(num = 1,
                        p_set = c(100, 100),
                        p_relevant_set = c(10, 10),
                        beta_set = c(0.5, 0.5),
                        num_nonpen = 0,
                        beta_nonpen = 0.7,
                        method = "elastic-net",
                        lambda = NULL,
                        bounds = NULL,
                        N = 21,
                        min.iter = 10,
                        seed = 1234,
                        parallel = FALSE,
                        verbose = TRUE,
                        dicotomise_second_layer = FALSE,
                        diagonal_covariance = FALSE,
                        fixed_alpha = 0.1) {

    if (!length(p_set) == 2 |
      !length(p_relevant_set) == 2 |
      !length(beta_set) == 2) {
    cat("Please specify bidimensional p_set, p_relevant_set,and beta_set
        vectors.\n")
    }
     
    set.seed(num)
  
    # Generate dataset
    sim <- simulation_data(
      p_set,
      p_relevant_set,
      beta_set,
      num_nonpen,
      beta_nonpen,
      n = 200,
      dicotomise_second_layer = dicotomise_second_layer,
      diagonal_covariance = diagonal_covariance
    )
    foldid <- sample(rep(seq(5),length=dim(sim$X)[1]/2))
    relevant_covariates<-sim$relevant_covariates
    p<-sim$p
    
    x<-scale(sim$X[1:100,])[,]
    y <- sim$Y[1:100]

    x_test <- scale(sim$X[101:200,])[,]
    y_test <- sim$Y[101:200]
    
    fit <-
      logisticIPFStructPenaltyReg(
        x,
        y,
        x_test,
        y_test,
        p,
        foldid,
        num.nonpen = num_nonpen,
        method = method,
        fixed_alpha = fixed_alpha,
        lambda = lambda,
        bounds = bounds,
        N = N,
        min.iter = min.iter,
        seed = seed,
        parallel = parallel,
        verbose = verbose
      )
    
    selected_covariates <- fit$selected
    
    cat("relevant covariates", relevant_covariates, "\n")
    cat("selected covariates", selected_covariates, "\n")
    
    fit$intersection <-
      sum(selected_covariates %in% relevant_covariates) - num_nonpen
    fit$precision <- fit$intersection /
      (length(selected_covariates) - num_nonpen)
    fit$recall <- fit$intersection /
      (length(relevant_covariates) - num_nonpen)
  
  return(fit)
}