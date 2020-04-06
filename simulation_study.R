################################################################################
############################### Simulation study ###############################
################################################################################

# This is an adaptation of the script by Zhi Zhao (zhi.zhao@medisin.uio.no) used
# for the analysis of simulated data in Zhao & Zucknick (2020) JRSSC.

# ==============================================================================
# The following settings are taken as input when running this on the HPC
# ==============================================================================

# Read input arguments from command line
# args <- commandArgs(trailingOnly = TRUE)
setting_number <- as.integer(args[1])
method <- args[2]
diagonal_covariance <- as.logical(args[3])
nonpenalised_covariates <- args[4]

# Covariance type
if(diagonal_covariance == TRUE){
    diag_string <- "diagonalcov"
} else {
    diag_string <- "nondiagonal"
}
cat("Covariance:", diag_string, "\n")

# Method
cat("Method:", method, "\n")
# Must be one of 
# "elastic-net", "sIPF-elastic-net", "separate-elastic-net", "fixed-alpha"

# Value of alpha, if it is fixed by the user
fixed_alpha <- 0.1
if(method == "fixed-alpha")
    cat("alpha =", fixed_alpha, "\n")

# Setting number
settings <- c("A", "B", "C", "D", "E", "F")
setting <- settings[setting_number]
cat("Setting:", setting, "\n")

# Non-penalised covariates
cat("Non-penalised covariates:", nonpenalised_covariates, "\n")

# Experiment number
# num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cat("Experiment number:", num, "\n")
# ==============================================================================

library(IPFStructPenalty)
source("logisticIPFStructPenalty.R") # Use alternative function modified by me

# Load the simulation functions
source("simulate_binomial_data.R")
source("simulation_function.R")

# Same settings as in Boulesteix et al. (2017)
p_sets <- matrix(c(1000, 1000, # Setting A
                   100, 1000,  # Setting B
                   100, 1000,  # Setting C
                   100, 1000,  # Setting D
                   20, 1000,   # Setting E
                   20, 1000),  # Setting F
                 ncol = 2,
                 byrow = T)

p_relevant_sets <- matrix(c(10, 10, # Setting A
                            3, 30,  # Setting B
                            10, 10, # Setting C
                            20, 0,  # Setting D
                            3, 10,  # Setting E
                            15, 3), # Setting F
                          ncol = 2,
                          byrow = T)

beta_values <- matrix(
  c(0.5, 0.5,  # Setting A
    0.5, 0.5,  # Setting B
    0.5, 0.5,  # Setting C
    0.3, 0,    # Setting D
    1, 0.3,    # Setting E
    0.5, 0.5), # Setting F
  ncol = 2,
  byrow = T
)

max_num <- 100
max_attempts <- 200
n_settings <- 6
n_methods <- length(methods)
n_layers <- dim(p_sets)[2]

lambda <- alpha <- array(NA, c(n_layers + 1))

names(alpha) <- c("Layer1", "Layer2", "Global")
names(lambda) <- c("Layer1", "Layer2", "Final")

p_set <- p_sets[setting_number, ]
cat("Number of covariates:", p_set, "\n")

p_relevant_set <- p_relevant_sets[setting_number,]
cat("Number of relevant covariates:", p_relevant_set, "\n")

beta_set <- beta_values[setting_number,]
cat("Values of the beta coefficients:", beta_set, "\n")

if(nonpenalised_covariates == "few"){
  num_nonpen = p_set[1]/10
  beta_nonpen = beta_set[1]
}else if(nonpenalised_covariates == "many"){
  num_nonpen = p_set[1]
  beta_nonpen = beta_set[1]
}else if(nonpenalised_covariates == "none"){
  num_nonpen = 0
  beta_nonpen = NA
}else{
  stop("Value of 'nonpenalised_covariates' can only be 'none', 'few', or
       'many.")
}
cat("Number of non-penalised covariates:", num_nonpen, "\n")
cat("Value of the beta coefficients for the non-pen cov:", beta_nonpen, "\n")

cat("Method:", method, "\n")

count <- 1
fit_successful <- FALSE

while (fit_successful == FALSE & count < max_attempts) {
    try(fit <- simulation_function(
        num,
        p_set,
        p_relevant_set,
        beta_set,
        num_nonpen = num_nonpen,
        beta_nonpen = beta_nonpen,
        method = method,
        parallel = FALSE,
        verbose = FALSE,
        seed = count,
        diagonal_covariance = diagonal_covariance,
        fixed_alpha = fixed_alpha
        )
    )
    
    # if fit exists AND is not null
    if (exists("fit"))
        if (!is.null(fit))
          fit_successful <- TRUE
    
    # if fails to fit change seed
    if (!fit_successful) {
        cat("Failed attempt number", count, "\n")
        count <- count + 1
    }
}

if (count <= max_attempts) {
    cvm <- fit$cvm
    cvm_cv <- fit$cvm_cv
    cv <- fit$cv
    
    if (method == "separate-elastic-net") {
        alpha <- c(fit$alpha, 0)
        lambda <- fit$lambda
    }else if (method == "fixed-alpha") {
        alpha <- c(fit$alpha, fit$alpha, 0)
        lambda <- fit$lambda
    }else{
        alpha <- c(NA, NA, fit$alpha)
        lambda <- c(NA, NA, fit$lambda)
    } # end if (count <= max_attempts)
    
    precision <- fit$precision
    recall <- fit$recall
} # end while (fit_successful == FALSE & count < max_attempts)

cat("Out-of-sample misclassification rate:", cvm, "\n")
cat("Within-sample misclassification rate:", cvm_cv, "\n")
cat("Value(s) of alpha:", alpha, "\n")
cat("Value(s) of lambda:", lambda, "\n")
cat("Number of selected variables:", cv, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("Intersection between relevant and selected:", fit$intersection, "\n")

# save(cvm, cvm_cv, alpha, cv, lambda,  precision, recall, count,
#      file =  paste("results/setting_", setting, "/nonpen_",
# 		   nonpenalised_covariates, "/method_",
# 		   method, "/", diag_string, "/", num,
#                    ".RData", sep = ""))