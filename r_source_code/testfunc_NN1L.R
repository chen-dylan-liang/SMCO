# Load the compiler package
library(compiler)

# ReLU activation function (Ensuring it returns a matrix)
relu <- function(x) {
  matrix(pmax(0, x), nrow = nrow(x), ncol = ncol(x))
}

# Initialize model parameters
initialize_nn <- function(input_dim, hidden_dim) {
  list(
    W1 = matrix(runif(input_dim * hidden_dim, -4, 4), nrow = input_dim, ncol = hidden_dim),
    b1 = matrix(runif(hidden_dim, 0, 8), nrow = 1, ncol = hidden_dim),
    W2 = matrix(runif(hidden_dim, -4, 4), nrow = hidden_dim, ncol = 1),
    b2 = matrix(runif(1,-4,4), nrow = 1, ncol = 1)
  )
}

# Forward pass
forward_pass <- function(model, x) {
  # Compute first layer activation
  Z1 <- x %*% model$W1 + matrix(rep(model$b1, nrow(x)), nrow = nrow(x), byrow = TRUE)
  A1 <- relu(Z1)  # Ensure relu() returns a matrix
  
  # Compute output layer
  Z2 <- A1 %*% model$W2 + matrix(rep(model$b2, nrow(A1)), nrow = nrow(A1), byrow = TRUE)
  
  list(Z1 = Z1, A1 = A1, Z2 = Z2)
}

# # Generate synthetic data
generate_data <- function(n, true_model, noise_sigma) {
  input_dim <- nrow(true_model$W1)
  hidden_dim <- ncol(true_model$W1) 
  x <- matrix(runif(n * input_dim, -4, 4), nrow = n, ncol = input_dim)
  y <- forward_pass(true_model, x)$Z2 + noise_sigma * rnorm(n) #sin(rowSums(x^2)) + log(abs(apply(x, 1, prod)) + 1)
  y <- matrix(y, nrow = n, ncol = 1)  # Reshape for matrix operations
  list(x = x, y = y)
}

model_to_param <- function(model) {
  input_dim = nrow(model$W1)
  hidden_dim = ncol(model$W1)
  W1vec = matrix(model$W1, nrow = input_dim*hidden_dim, ncol = 1)
  b1vec = matrix(model$b1, nrow = hidden_dim, ncol = 1)
  W2vec = matrix(model$W2, nrow = hidden_dim, ncol = 1)
  params = rbind(W1vec, b1vec, W2vec, model$b2)
  return(params)
}

param_to_model <- function(params, input_dim, hidden_dim) {
  d = length(params)
  d_W1 = input_dim * hidden_dim
  d_b1 = hidden_dim
  d_W2 = hidden_dim
  d_b2 = 1
  if (d != d_W1 + d_b1 + d_W2 + d_b2) {
    print("Error: dimension of parameter does not match model dimension")
    break
  } else {
    W1vec = params[1:d_W1]
    W1 = matrix(W1vec, nrow = input_dim, ncol = hidden_dim)
    b1vec = params[(d_W1+1):(d_W1+d_b1)]
    b1 = matrix(b1vec, nrow = 1, ncol = hidden_dim)
    W2 = matrix(params[(d_W1+hidden_dim+1):(d-1)], nrow = hidden_dim, ncol = 1)
    b2 = matrix(params[d], nrow = 1, ncol = 1)
    
    model = list(
      W1 = W1,
      b1 = b1,
      W2 = W2,
      b2 = b2
    )
    
    return(model)
  }
}

MSE_ANN <- function(params, input_dim, hidden_dim, x, y) {
  model <- param_to_model(params=params, input_dim=input_dim, hidden_dim=hidden_dim)
  model_predict <- forward_pass(model, x)$Z2 
  model_MSE <- mean((model_predict-y)^2)
  return(model_MSE)
}


setup_ANN <- function(input_dim, hidden_dim, n, use_truth, noise_sigma, seed) {
  # Initialize model
  true_model <- initialize_nn(input_dim, hidden_dim)
  
  # Generate synthetic data
  data <- generate_data(n, true_model, noise_sigma)
  
  # Compute true model prediction
  true_predict <- forward_pass(true_model, data$x)$Z2
  
  return(list(
    data = data,
    true_model = true_model,
    true_predict = true_predict
  ))
}

config_ANN <- function(input_dim, hidden_dim, n, use_truth, noise_sigma, seed) {
  
  ANN_setup <- setup_ANN(input_dim, hidden_dim, n, use_truth, noise_sigma, seed)
  data <- ANN_setup$data
  y_true <- if (use_truth) ANN_setup$true_predict else data$y
  dim = input_dim * hidden_dim + hidden_dim + hidden_dim + 1
  
  func <- function(params) MSE_ANN(params, input_dim=input_dim, hidden_dim=hidden_dim, x=data$x, y=y_true)
  func <- cmpfun(func)
  f_min = if (use_truth) 0 else noise_sigma^2
  
  config <- list(
    input_dim = input_dim,
    hidden_dim = hidden_dim,
    n = n,
    use_truth = use_truth,
    noise_sigma = noise_sigma,
    dim = dim,
    data = ANN_setup$data,
    true_predict = ANN_setup$true_predict,
    true_model = ANN_setup$true_model,
    f = func,
    f_min = f_min,
    f_max = NA,
    bounds_lower = rep(-10, dim),
    bounds_upper = rep(10, dim)
  )
  return(config)
}


# Compile the core computation functions
relu <- cmpfun(relu)
forward_pass <- cmpfun(forward_pass)
param_to_model <- cmpfun(param_to_model)
MSE_ANN <- cmpfun(MSE_ANN)

# Run test case using SBMS

if (F) {
  
  set.seed(123)
  
  # ANN Configuration
    input_dim = 3            # Number of input features
    hidden_dim = 5           # Number of hidden neurons
    n = 1000                 # Number of data points simulated
    use_truth = F            # Use simulated outcome y
    noise_sigma = 1          # Standard deviation of noise in y
    
    MSE_true <- if (use_truth) 0 else MSE_true = noise_sigma^2
    
    ANN_config <- config_ANN(input_dim = input_dim, hidden_dim = hidden_dim, n = n,
                             use_truth = use_truth, noise_sigma = noise_sigma, seed = NULL)
    
    func_to_max <- function(params) {- ANN_config$f(params)}
  
  # SBMS Configuration
    source("SBMS_250228.R")
    SBMS_options <- list(
      seed = NULL,
      use_parallel = F,
      bounds_buffer = 0.05,
      buffer_rand = F,
      tol_convergence = 1e-5,
      refine_search = T,
      refine_ratio = 0.25,
      boost_iterations = 1000,
      partial_option = "center",
      max_iterations = 300,
      n_starts = 150
    )
    
    # Run SBMS Optimization
    time_start = Sys.time()
    result_multi <- SBMS_multi(f = func_to_max,
                               bounds_lower = bounds_lower,
                               bounds_upper = bounds_upper,
                               opt_control = SBMS_options)
    best_result <- result_multi$best_result
    time_end = Sys.time()
    
    time_run = difftime(time_end, time_start, units = "secs")
    MSE_min_found <- -best_result$f_optimal
    params_found <- best_result$x_optimal
    
    # Saving result to file
    filename = paste0("test_NN1L_in", input_dim, "_hi", hidden_dim, ".txt")
    sink(filename)
    
    cat("ANN Configuration:\n")
    cat("input_dim = ", input_dim, "\n")
    cat("hidden_dim = ", hidden_dim, "\n")
    cat("params_dim = ", d_params, "\n")
    cat("n = ", n, "\n")
    cat("use_truth = ", use_truth, "\n")
    cat("noise_sigma = ", noise_sigma, "\n")
    cat("\n")
    
    cat("SBMS Configuration:\n")
    print(SBMS_options)
    
    cat("\n")
    cat("Performance:\n")
    cat("time_run = ", time_run, "\n")
    cat("\n")
    cat("MSE_min_true = ", MSE_true, "\n")
    cat("MSE_min_found = ", MSE_min_found, "\n")
    cat("MSE_diff = ", abs(MSE_min_found - MSE_true), "\n")
    
    sink()
    
    cat("Optimization Finished:\n")
    cat("time_run = ", time_run, "\n")
    cat("\n")
    cat("MSE_min_true = ", MSE_true, "\n")
    cat("MSE_min_found = ", MSE_min_found, "\n")
    cat("MSE_diff = ", abs(MSE_min_found - MSE_true), "\n")
}
