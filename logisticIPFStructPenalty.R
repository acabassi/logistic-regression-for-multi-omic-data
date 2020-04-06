#' logisticIPFStructPenalty
#' @title Structured penalized logistic regression 
#' @description
#' Function producing results of the structured penalized logistic regression
#'
#' @importFrom Matrix Diagonal bdiag
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats sd
#' @importFrom parallel detectCores makeCluster	stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom penalized penalized
#'
#' @param x,y \code{x} is the input design matrix; \code{y} is the input
#' response matrix.
#' @param x_test,y_test \code{x} is the input design matrix for validation;
#' \code{y} is the input response matrix for validation. If these are not
#' provided, \code{x} and \code{y} are used instead.
#' @param p Number of predictors in each data source (excluding the
#' non-penalised ones).
#' @param foldid Vector of values for the cross-validation, its length must be
#' equal to the number of observations in the training set.
#' @param num.nonpen Number of non-penalised covariates.
#' @param method Method used to optimize the penalty parameters. The penalty
#' parameters of \code{elastic-net}, \code{IPF-lasso}, \code{sIPF-elastic-net}
#' are optimzed by the EPSGO algorithm. The penalty parameter of \code{lasso} is
#' optimzed by cross-validation. The default method is \code{IPF-lasso}.
#' @param fixed_alpha Value of alpha for the fixed alpha method. Default is 0.1.
#' @param lambda Optional user-supplied \code{lambda} sequence; default is NULL,
#' and \code{epsgo} chooses its own sequence.
#' @param bounds Bounds for the interval-searching parameters.
#' @param search.path Boolean indicating whether to save the visited points;
#' default is \code{FALSE}.
#' @param EI.eps Convergence threshold for the expected improvement between fmin
#' and the updated point
#' @param fminlower Minimal value for the function Q.func, default is 0.
#' @param N Number of start points; depends on the dimensionality of the
#' parameter space.
#' @param min.iter Minimum number of iterations after the initial \code{N}
#' iterations.
#' @param seed Random seed.
#' @param parallel If \code{TRUE}, use parallel foreach to fit each fold except
#' parallelizing each lambda for the tree-lasso methods. If \code{c(TRUE,TRUE)},
#' use parallel foreach to fit each fold and each lambda.
#' @param verbose Print the middle search information, default is \code{TRUE}.
#' @return An object of list "\code{IPFStructPenaltyReg}" is returned:
#'  \item{cvm}{the mean cross-validated error}
#'  \item{cvm_cv}{the mean cross-validated error if providing external dataset
#'  "\code{x_test}" and "\code{y_test}".}
#'  \item{alpha}{optimized \code{alpha}}
#'  \item{lambda}{optimized \code{lambda}}
#'  \item{pred}{the prediction of the responses}
#'  \item{ipf}{optimzed penalty factors}
#'  \item{Beta}{estimate of the coefficients}
#'  \item{cv}{number of nonzero coefficients}
#' @references Zhao, Z. & Zucknick, M. (2019). \emph{Stuctured penalized
#' regression for drug sensitivity prediction.} arXiv: 1905.00095.
#' @export

library(Matrix)
library(stats)
library(parallel)
library(foreach)
library(doParallel)
library(penalized)
library(glmnet)

misclassification_rate <- function(y_true, y_pred){
  sum(1 - (y_true == y_pred))/dim(y_true)[1]
}

logistic_prediction = function(x, beta) {
  as.numeric(1 / (1 + exp(-(x %*% beta))) > 0.5)
}

logisticIPFStructPenaltyReg <-
  function(x,
           y,
           x_test = NULL,
           y_test = NULL,
           p,
           foldid,
           num.nonpen = 0,
           method = "IPF-lasso",
           fixed_alpha = NULL,
           lambda = NULL,
           bounds = NULL,
           search.path = FALSE,
           EI.eps = 0.01,
           fminlower = 0,
           N = NULL,
           min.iter = 20,
           family = "binomial",
           seed = 1234,
           parallel = FALSE,
           verbose = TRUE,
           ...) {
    
  if((method!="lasso") & (method!="fixed-alpha") & is.null(bounds)){
    
    if(method=="elastic-net" | method == "separate-elastic-net"){
        bounds <- t(data.frame(alpha=c(0,1)))
    }else if(method=="IPF-lasso"){
        # the default setting for bounds is only for two data sources
        bounds <- t(data.frame(ipf=c(0.1,10)))
    }else if(method=="sIPF-elastic-net"){
        bounds <- t(data.frame(alpha=c(0,1), ipf1 = c(0.1, 10)))
    }else {
        cat("Ooops! Please give searching bounds for EPSGO algorithm!")
    }
    
    colnames(bounds) <- c("lower", "upper")
  }
  
  # Training set 
  x <- data.matrix(x)
  y <- data.matrix(y)
  
  # Test set (if not provided, it is the same as the training set)
  if(is.null(x_test) & is.null(y_test)){
    x_test <- x
    y_test <- y
  }else{
    x_test <- data.matrix(x_test)
    y_test <- data.matrix(y_test)
  }
  
  # Build penalty factors
  Beta <- as.vector(rbind(matrix(0, nrow = num.nonpen, ncol = 1),
                          matrix(10, ncol = 1, nrow = sum(p, na.rm = T))))
  adpen <- rep(1, dim(x)[2])
  adpen[Beta == 0] <- 0
  # This is just a vector with 1 + num.nonpen zeros and sum(p) ones
  
  #=======
  # LASSO
  #=======
  if (method == "lasso") {
    
    # Select lambda values
    lambda <- glmnet(x, y, family = "binomial", nlambda = 20, intercept = F,
                     standardize.response = F, penalty.factor = adpen)$lambda
    
    # Use cross-validation to find optimal value of lambda
    la <- cv.glmnet(x, y, family = "binomial", lambda = lambda, intercept = F,
                    foldid = foldid, standardize.response = F,
                    penalty.factor = adpen, parallel = parallel)
    
    # Save alpha, lambda, and in-sample misclassification error
    mr_cv <- min(la$cvm)
    alpha <- 1
    lambda <- la$lambda.min
    
    # Get values of the estimated coefficients
    Beta <- glmnet(x, y, family = "binomial", intercept = F,
                  lambda = seq(lambda*0.8, lambda*1.2, length = 11),
                  standardize.response = F, penalty.factor = adpen)$beta[,6]
    
    # Get predicted values for the test set
    y_pred <- logistic_prediction(x_test, Beta)
    
    # Compute out-of-sample misclassification error
    mr_val <- misclassification_rate(y_test, y_pred)
    
    # Save number of selected variables
    vs <- sum(Beta != 0)
  }
  
  #=============
  # Elastic-net
  #=============
  if (method == "elastic-net") {
    
    # Use EPSGO algorithm to tune alpha and lambda
    fit <- epsgo(Q.func = "tune.glmnet.interval", bounds = bounds,
                 lambda = lambda, N = N, parms.coding = "none", seed = seed,
                 fminlower = fminlower, x = x, y = y, num.nonpen = num.nonpen,
                 family = "binomial", foldid = foldid, type.min = "lambda.min",
                 p = p, intercept = T, standardize.response = F,
                 type.measure = "class", min.iter = min.iter, verbose = verbose,
                 parallel = parallel, EI.eps = EI.eps)
    sumint <- summary(fit, verbose = F)

    # Save alpha, lambda, and in-sample misclassification error
    mr_cv <- sumint$opt.error
    alpha <- sumint$opt.alpha
    lambda <- sumint$opt.lambda
    
    # Standardise penalty factors 
    adpen <- adpen*(dim(x)[2])/sum(adpen)
    
    # Get values of the estimated coefficients
    Beta <- glmnet(x, y, family="binomial", intercept=F,
                   lambda=seq(lambda*0.8,lambda*1.2, length=11),
                   standardize.response=F, penalty.factor = adpen)$beta[,6]

    # Get predicted values for the test set

    y_pred <- logistic_prediction(x_test, Beta)

    # Compute out-of-sample misclassification error
    mr_val <- misclassification_rate(y_test, y_pred)

    # Save number of selected variables
    vs <- sum(Beta!=0)
  }

  #======================
  # Separate elastic-net
  #======================
  if(method == "separate-elastic-net") {
    
    colnames(bounds) <- c("lower", "upper")
    alpha <- vector(mode = "numeric", length = length(p))
    lambda_temp <- vector(mode = "numeric", length = length(p) + 1)
    vs <- 0
    selected <- c()
    
    # For each layer
    for(i in 1:length(p)){
      cat("Layer", i, "\n")
      
      # Find first element of layer
      if(i>1){
        start <- (sum(p[1:(i - 1)]) + num.nonpen + 1)
      }else{
        start<- num.nonpen + 1
      }
      
      # Find last element of layer
      end <- sum(p[1:i]) + num.nonpen
      
      # Use EPSGO algorithm to tune alpha_m and lambda_m
      fit <- epsgo(
        Q.func = "tune.glmnet.interval",
        bounds = bounds,
        lambda = lambda,
        N = N,
        parms.coding = "none",
        seed = seed,
        fminlower = fminlower,
        x = x[, start:end],
        y = y,
        family = "binomial",
        foldid = foldid,
        type.min = "lambda.min",
        p = p[i],
        intercept = T,
        standardize.response = F,
        type.measure = "class",
        min.iter = min.iter,
        verbose = verbose,
        parallel = parallel,
        EI.eps = EI.eps
      )
      sumint <- summary(fit, verbose=F)
      
      # Save alpha_m and lambda_m
      alpha[i] <- sumint$opt.alpha
      lambda_temp[i] <- sumint$opt.lambda
      
      # Get values of the estimated coefficients (just to check which ones are
      # zero)
      Beta_temp <-
        glmnet(
          x[, start:end],
          y,
          family = "binomial",
          intercept = F,
          alpha = alpha[i],
          lambda = seq(lambda_temp[i] * 0.8, lambda_temp[i] *
                         1.2, length = 11),
          standardize.response = F
        )$beta[, 6]
      
      # Save indices of selected elements of beta
      selected <- c(selected, which(Beta_temp != 0) + start - 1)
    }
    
    # Save total number of selected variables (excluding non-penalised ones)
    vs <- length(selected)
    
    # If necessary, add non-penalised features to the set of features used in
    # the final step
    if(num.nonpen>0)
      selected <- c(1:num.nonpen, selected)
    
    # Select lambda values
    lambda <-
      glmnet(
        x[, selected],
        y,
        family = "binomial",
        nlambda = 20,
        intercept = F,
        standardize.response = F,
        alpha = 0,
        penalty.factor = adpen
      )$lambda
    
    # Use cross-validation to select optimal value of lambda
    ridge <-
      cv.glmnet(
        x[, selected],
        y,
        family = "binomial",
        lambda = lambda,
        intercept = F,
        foldid = foldid,
        standardize.response = F,
        parallel = parallel,
        alpha = 0,
        type.measure = "class",
        penalty.factor = adpen
      )
    
    # Save in-sample misclassification error
    mr_cv <- min(ridge$cvm)
    
    # Save value of lambda used in the final step
    lambda_temp[i+1] <- ridge$lambda.min
    
    # Get estimated values of the coefficients in the final step
    Beta <-
      glmnet(
        x[, selected],
        y,
        family = "binomial",
        intercept = F,
        lambda = seq(lambda_temp[i+1]*0.8, lambda_temp[i+1]*1.2, length = 11),
        standardize.response = F,
        alpha = 0,
        penalty.factor = adpen
      )$beta[, 6]
    
    # Get predicted values for the test set
    y_pred <- logistic_prediction(x_test[,selected], Beta)
    
    # Compute out-of-sample misclassification error
    mr_val <- misclassification_rate(y_test, y_pred)
    
    # Save all values of lambda (one for each layer + the final one)
    lambda <- lambda_temp
  }
  
  #====================
  # sep-EN fixed alpha
  #====================
  
  if (method == "fixed-alpha") {
    
    alpha <- fixed_alpha
    
    lambda_temp <- vector(mode = "numeric", length = length(p) + 1)
    vs <- 0
    selected <- c()
    
    for (i in 1:length(p)) {

      # Find first element of layer
      if (i > 1) {
        start <- (sum(p[1:(i - 1)]) + num.nonpen + 1)
      } else{
        start <- num.nonpen + 1
      }
      
      # Find last element of layer
      end <- sum(p[1:i]) + num.nonpen
      
      # Select lambda values
      lambda <-
        glmnet(
          x[, start:end],
          y,
          family = "binomial",
          nlambda = 20,
          intercept = F,
          standardize.response = F,
          alpha = alpha,
        )$lambda
      
      # Use cross-validation to select optimal value of lambda
      fit_temp <-
        cv.glmnet(
          x[, start:end],
          y,
          family = "binomial",
          lambda = lambda,
          intercept = F,
          foldid = foldid,
          standardize.response = F,
          parallel = parallel,
          alpha = alpha,
          type.measure = "class"
        )
      
      lambda_temp[i] <- fit_temp$lambda.min
      
      # Get coefficients
      Beta_temp <-
        glmnet(
          x[, start:end],
          y,
          family = "binomial",
          intercept = F,
          lambda = seq(lambda_temp[i]*0.8, lambda_temp[i]*1.2, length = 11),
          standardize.response = F,
          alpha = alpha,
        )$beta[, 6]
      
      selected <- c(selected, which(Beta_temp != 0) + start - 1)
    }
    
    # Save total number of selected variables (excluding non-penalised ones)
    vs <- length(selected)
    
    # If necessary, add non-penalised features to the set of features used in
    # the final step
    if(num.nonpen>0)
      selected <- c(1:num.nonpen, selected)
    
    # Select lambda values
    lambda <-
      glmnet(
        x[, selected],
        y,
        family = "binomial",
        nlambda = 20,
        intercept = F,
        standardize.response = F,
        alpha = 0,
        penalty.factor = adpen
      )$lambda
    
    # Use cross-validation to select optimal value of lambda
    ridge <-
      cv.glmnet(
        x[, selected],
        y,
        family = "binomial",
        lambda = lambda,
        intercept = F,
        foldid = foldid,
        standardize.response = F,
        parallel = parallel,
        alpha = 0,
        type.measure = "class",
        penalty.factor = adpen
      )
    
    # Save in-sample misclassification error
    mr_cv <- min(ridge$cvm)
    
    # Save value of lambda used in the final step
    lambda_temp[i+1] <- ridge$lambda.min
    
    # Get estimated values of the coefficients in the final step
    Beta <-
      glmnet(
        x[, selected],
        y,
        family = "binomial",
        intercept = F,
        lambda = seq(lambda_temp[i+1]*0.8, lambda_temp[i+1]*1.2, length = 11),
        standardize.response = F,
        alpha = 0,
        penalty.factor = adpen
      )$beta[, 6]
    
    # Get predicted values for the test set
    y_pred <- logistic_prediction(x_test[,selected], Beta)
    
    # Compute out-of-sample misclassification error
    mr_val <- misclassification_rate(y_test, y_pred)
    
    # Save all values of lambda (one for each layer + the final one)
    lambda <- lambda_temp
  }
  
  #====================
  # IPF-lasso, sIPF-EN
  #====================
  if((method == "IPF-lasso") | (method == "sIPF-elastic-net")){
    
    # Use EPSGO algorithm to tune the regression parameters
    fit <-
      epsgo(
        Q.func = "tune.glmnet.interval",
        bounds = bounds,
        lambda = lambda,
        N = N,
        parms.coding = "none",
        seed = seed,
        fminlower = fminlower,
        x = x,
        y = y,
        num.nonpen = num.nonpen,
        family = "binomial",
        foldid = foldid,
        type.min = "lambda.min",
        p = p,
        intercept = T,
        standardize.response = F,
        type.measure = "class",
        min.iter = min.iter,
        verbose = verbose,
        parallel = parallel,
        EI.eps = EI.eps
      )
    sumint <- summary(fit, verbose=F)
    
    # Save alpha, lambda, and in-sample misclassification error
    mr_cv <- sumint$opt.error
    alpha <- sumint$opt.alpha
    lambda <- sumint$opt.lambda
    
    # Get penalty factors
    ipf <- sumint$opt.ipf
    
    # Rescale penalty factors
    beta <- matrix(0, nrow = num.nonpen + sum(p), ncol = 1)
    for (s in 1:length(p))
      beta[num.nonpen + ifelse(s == 1, 0, sum(p[1:(s - 1)])) + 1:p[s], ] <-
      10 + 10 * s
    adpen <- rep(1, dim(x)[2])
    adpen[as.vector(beta) == 0] <- 0
    adpen[as.vector(beta) == 20] <- 1
    for (i in 1:length(ipf))
      adpen[as.vector(beta) == (10 + 10 * (i + 1))] <- ipf[i]
    adpen <- adpen * (dim(x)[2]) / sum(adpen)
    
    # Get values of the estimated coefficients
    Beta <-
      glmnet(
        x = x,
        y = y,
        family = "binomial",
        alpha = alpha,
        offset = NULL,
        lambda = seq(lambda * 0.8, lambda * 1.2, length = 11),
        penalty.factor = adpen,
        intercept = F,
        standardize.response = F
      )$beta[, 6]
    
    # Get predicted values for the test set
    y_pred <- logistic_prediction(x_test, Beta)
    
    # Compute out-of-sample misclassification error
    mr_val <- misclassification_rate(y_test, y_pred)
    
    # Save number of predicted variables
    vs <- sum(Beta != 0)
  }

  output = list(cvm = mr_val, # Out-of-sample misclassification rate
                cvm_cv = mr_cv, # Within-sample misclassification rate
                alpha = alpha, # Value(s) of alpha
                lambda = lambda, # Value(s) of lambda
                pred = y_pred, # Estimated values of Y for the test set
                Beta = Beta, # Regression coefficients
                cv = vs) # Number of selected variables (excluding the
                         # non-penalised ones)
  
  if ((method != "lasso") & (method != "fixed-alpha")) {
    output$search.fit = fit
  }
  
  if ((method == "IPF-lasso") | (method == "sIPF-elastic-net")) {
    output$ipf = ipf
  }
  
  if ((method == "separate-elastic-net") | (method == "fixed-alpha")) {
    output$selected <- selected
  }else{
    output$selected <- which(Beta!=0)
  }
  
  return(output)
}
