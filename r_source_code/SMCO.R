# Strategic Monte Carlo Optimization (SMCO) Algorithm

# "Optimization via Strategic Law of Large Numbers"
# By: Xiaohong Chen, Zengjing Chen, Wayne Yuan Gao, Xiaodong Yan, Guodong Zhang, and Yu Zhang
# Date: July 15, 2025
# GitHub Repository Maintained by: Wayne Yuan Gao

# This Rscript is a self-contained R Implementation of the SMCO Algorithm
# This Rscript contains the following functions:

# Core optimizer functions: 
# (1) "SMCO_single": executes the basic iterations of the SMCO algorithm
# (2) "SMCO_single_refine": calls SMCO_single to run refined search (with options)
# (3) "SMCO_single_boost": calls SMCO_refine to run additional boosted search (with ptions)
# (4) "SMCO_multi": runs SMCO_single_boost from multiple starting points and return the best result
#     Note: The options/arguments of the SMCO functions can be set to implement the raw "SMCO_single" algorithm
#           without additional refinement/boosting.

# Helper functions:
# (5) "generate_sobol_points": generate multiple Sobol starting points
# (6) "check_bounds": check and enforce bound constraints
# (7) "compute_partial_signs": compute partial finite differences used in the SMCO algorithm

# Example usage at the end of Rscript

#####################################
# Generate diverse Sobol starting points
generate_sobol_points <- function(n_starts, bounds_lower, bounds_upper, seed = NULL) {
  library(qrng)
  if (!is.null(seed)) set.seed(seed)
  
  d = length(bounds_lower)
  sobol_points <- sobol(n_starts, d = d)
  points <- matrix(NA, nrow = n_starts, ncol = d)
  for (j in 1:d) {
    points[,j] <- bounds_lower[j] + sobol_points[,j] * (bounds_upper[j] - bounds_lower[j])
  }
  return(points)
}

# Helper function to check bounds
check_bounds <- function(x, bounds_lower, bounds_upper) {
  is_out <- (sum((x > bounds_upper) + (x < bounds_lower)) > 0)
  x_in <- if (is_out) pmax(pmin(x, bounds_upper), bounds_lower) else x
  return(list(is_out = is_out, x_in = x_in))
}

# Helper function to:
# (1) Compute discrete partial derivative signs; and
# (2) Return max function evaluation in addition
compute_partial_signs <- function(f, x, fx, h_step, bounds_lower, bounds_upper, 
                                  partial_option = "center", use_runmax  = T) {
  d <- length(x)
  
  if (is.null(h_step)) {
    h_step <- rep(1e-8, d)
  }
  
  partial_signs <- numeric(d)
  f_partial_best <- fx
  x_partial_best <- x
  
  if (partial_option == "center") { # Do two-sided perturbation
      
      x_plus <- x
      x_minus <- x
    
      for(j in 1:d) {
        x_plus[j] <- min(x[j] + h_step[j], bounds_upper[j])
        x_minus[j] <- max(x[j] - h_step[j], bounds_lower[j])
        
        f_plus <- f(x_plus)
        f_minus <- f(x_minus)
        sign_f_diff <- (f_plus > f_minus)
        partial_signs[j] <- sign_f_diff
        
        if (use_runmax && sign_f_diff && (f_plus > f_partial_best)) {
          f_partial_best <- f_plus
          x_partial_best <- x_plus
        } else if (use_runmax && (!sign_f_diff) && (f_minus > f_partial_best)) {
          f_partial_best <- f_minus
          x_partial_best <- x_minus
        }
        
        # Reset vectors for next j
        x_plus[j] <- x[j]
        x_minus[j] <- x[j]
      }
  } else { # if partial_option = "forward", do one-sided perturbation

    x_perturb <- x
    for(j in 1:d) {
      if (x[j] >= bounds_upper[j]) {
        x_perturb[j] <- max(x[j] - h_step[j], bounds_lower[j])
        f_perturb <- f(x_perturb)
        sign_f_diff <- (fx > f_perturb)
      } else {
        x_perturb[j] <- min(x[j] + h_step[j], bounds_upper[j])
        f_perturb <- f(x_perturb)
        sign_f_diff <- (f_perturb > fx)
      }
      
      partial_signs[j] <- sign_f_diff
      if (use_runmax && (f_perturb > f_partial_best)) {
        f_partial_best <- f_perturb
        x_partial_best <- x_perturb
      }
      
      # Reset vector for next j
      x_perturb[j] <- x[j]
    }
  }
  
  return(list(signs = partial_signs, f_partial_best = f_partial_best, x_partial_best = x_partial_best))
}

# Optimized single-start SMCO function
SMCO_single <- function(f, bounds_lower, bounds_upper, start_point,
                        bounds_buffer, buffer_rand, iter_max, iter_nstart, iter_boost, 
                        tol_conv, partial_option, use_runmax) 
{
  d <- length(bounds_lower)
  
  # Initialize
  x_current <- start_point
  f_current <- f(x_current)
  f_runmax <- f_current
  x_runmax <- x_current
  
  # Initialize loop counters and running sum S
  bounds_diff <- bounds_upper - bounds_lower
  n_boost_1 <- iter_boost + iter_nstart
  n_boost_max <- iter_boost + iter_nstart + iter_max
  
  S <- start_point * n_boost_1
  
  # Main optimization loop
  for (n in n_boost_1:n_boost_max) {
    
    # Compute discrete partial derivative signs
    h_step <- bounds_diff / (n + 1)
    partial_sign_results <- compute_partial_signs(f, x_current, f_current, h_step, 
                                                  bounds_lower, bounds_upper, 
                                                  partial_option, use_runmax)
    
    # Implement deterministic pushout or random deviation from bounds
    bound_out <- bounds_buffer * bounds_diff
    if (buffer_rand) { # for random deviation from bounds
      bound_out <- bound_out * runif(d, -1, 1)
    }
    bounds_upper_out <- bounds_upper + bound_out
    bounds_lower_out <- bounds_lower - bound_out
    
    # Update sum S and point x
    Z <- partial_sign_results$signs * bounds_upper_out + (1-partial_sign_results$signs) * bounds_lower_out
    S <- S + Z
    x_next <- S / (n + 1)
    f_next <- f(x_next)
    
    if (use_runmax) {
      # Update the running max record
      f_partial_best <- partial_sign_results$f_partial_best
      f_next_best <- max(f_partial_best, f_next)
      if (f_next_best > f_runmax) {
        f_runmax <- f_next_best
        x_runmax <- if (f_partial_best > f_next) partial_sign_results$x_partial_best else x_next
      }
    }
    # Update x_current and f_current
    f_prev <- f_current
    f_current <- f_next
    x_current <- x_next
    
    # Check convergence after at least 50% of allocated iterations
    if ((n >= n_boost_1 + iter_max/2) && (abs(f_current - f_prev) < tol_conv)) {
      break
    }
  }
  
  if (use_runmax) {
    return(list(
      x_optimal = x_current,
      f_optimal = f_current,
      iterations = n - iter_boost,
      f_runmax = f_runmax,
      x_runmax = x_runmax
    ))} else {
    return(list(
      x_optimal = x_current,
      f_optimal = f_current,
      iterations = n - iter_boost
    ))
    }
}

# Single-start with refinement search
SMCO_single_refine <- function(f, bounds_lower, bounds_upper, start_point,
                               bounds_buffer, buffer_rand, iter_max, iter_nstart, iter_boost, tol_conv, 
                               refine_search, refine_ratio, partial_option, use_runmax) 
{ 
  # Run initial search
  if (refine_search == F) { refine_ratio = 0 }
  iter_max_initial <- round(iter_max * (1 - refine_ratio), digits = 0)
  
  result <- SMCO_single(f, bounds_lower, bounds_upper, start_point = start_point,
                        bounds_buffer = bounds_buffer, buffer_rand = buffer_rand, 
                        iter_max = iter_max_initial, iter_nstart = iter_nstart, iter_boost = iter_boost, 
                        tol_conv = tol_conv, partial_option = partial_option,
                        use_runmax = use_runmax)
  
  # Handle boundary constraints
  check_optimal <- check_bounds(result$x_optimal, bounds_lower, bounds_upper)
  if (check_optimal$is_out) {
    x_optimal_in <- check_optimal$x_in
    f_optimal_in <- f(x_optimal_in)
    result$x_optimal <- x_optimal_in
    result$f_optimal <- f_optimal_in
  } else {
    x_optimal_in <- result$x_optimal
    f_optimal_in <- result$f_optimal
  }
  
  if (use_runmax) {
    check_runmax <- check_bounds(result$x_runmax, bounds_lower, bounds_upper)
    if (check_runmax$is_out) {
      x_runmax_in <- check_runmax$x_in
      f_runmax_in <- f(x_runmax_in)
      result$x_runmax <- x_runmax_in
      result$f_runmax <- f_runmax_in
    } else {
      x_runmax_in <- result$x_runmax
      f_runmax_in <- result$f_runmax
    }
  }
  
  # Run refined search
  if (refine_search) {
    iter_max_refine <- round(iter_max * refine_ratio, digits = 0)
    
    start_point_refine <- if (use_runmax && (f_runmax_in > f_optimal_in)) x_runmax_in else x_optimal_in
    iter_boost_refine <- iter_boost + 1000
    
    refine_result <- SMCO_single(f, bounds_lower, bounds_upper, 
                                 start_point = start_point_refine,
                                 bounds_buffer = 0, buffer_rand = buffer_rand, 
                                 iter_max = iter_max_refine, iter_nstart = iter_nstart, iter_boost = iter_boost_refine, 
                                 tol_conv = tol_conv, partial_option = partial_option, 
                                 use_runmax = use_runmax)
    
    if (use_runmax) {
      # Handle boundary constraints for refined result
      f_refine_runmax <- refine_result$f_runmax
      x_refine_runmax <- refine_result$x_runmax
      check_refine_runmax <- check_bounds(x_refine_runmax, bounds_lower, bounds_upper)
      
      if (check_refine_runmax$is_out) {
        x_refine_runmax <- check_refine_runmax$x_in
        f_refine_runmax <- f(x_refine_runmax)
      }
      
      if (f_refine_runmax > refine_result$f_optimal) {
        refine_result$f_optimal <- f_refine_runmax
        refine_result$x_optimal <- x_refine_runmax
      }
    }
    
  } else {
    
    # If no refine_search, use the original results with bounds checking
    refine_result <- result
    
    
    
    if (use_runmax) {
      refine_result$f_runmax <- f_runmax_in
      refine_result$x_runmax <- x_runmax_in
      
      if (f_runmax_in > f_optimal_in) {
        refine_result$f_optimal <- f_runmax_in
        refine_result$x_optimal <- x_runmax_in
      } else {
        refine_result$f_optimal <- f_optimal_in
        refine_result$x_optimal <- x_optimal_in
      }
    }
  }
  
  return(refine_result)
}

# Single-start with/without additional boosted search
SMCO_single_boost <- function(f, bounds_lower, bounds_upper, start_point,
                              bounds_buffer, buffer_rand, iter_max, iter_nstart, iter_boost, tol_conv, 
                              refine_search, refine_ratio, partial_option, use_runmax)
{
  
  # Run regular search with non-boosted starting index n = 1
  refine_result <- SMCO_single_refine(f, bounds_lower, bounds_upper, start_point = start_point,
                                      bounds_buffer = bounds_buffer, buffer_rand = buffer_rand, 
                                      iter_max = iter_max, iter_nstart = iter_nstart, iter_boost = 0, 
                                      tol_conv = tol_conv, refine_search = refine_search, refine_ratio = refine_ratio, 
                                      partial_option = partial_option, use_runmax = use_runmax)
  
  # Run additional boosted search with "boosted" starting index n = iter_boost
  if (iter_boost > 0) {
    refine_result_boost <- SMCO_single_refine(f, bounds_lower, bounds_upper, start_point = start_point,
                                              bounds_buffer = bounds_buffer, buffer_rand = buffer_rand, 
                                              iter_max = iter_max, iter_nstart = iter_nstart, 
                                              iter_boost = iter_boost, 
                                              tol_conv = tol_conv, refine_search = refine_search, refine_ratio = refine_ratio, 
                                              partial_option = partial_option, use_runmax = use_runmax)
    
    if (refine_result_boost$f_optimal > refine_result$f_optimal) {
      refine_result <- refine_result_boost
    }
  }
  
  return(refine_result)
}

# Optimized multi-start function
SMCO_multi <- function(f, bounds_lower, bounds_upper, start_points = NULL,
                       opt_control = list(
                         n_starts = 100,
                         iter_max = 200,
                         iter_nstart = 1,
                         iter_boost = 0,
                         bounds_buffer = 0.05,
                         buffer_rand = F,
                         tol_conv = 1e-8,
                         refine_search = T,
                         refine_ratio = 0.5,
                         partial_option = "center",
                         use_runmax = T,
                         use_parallel = F,         
                         seed = 123)
)
{
  # Generate diverse starting points
  # if (!is.null(opt_control$seed)) {
  #   set.seed(opt_control$seed)
  # }
  
  # if (is.null(start_points)) {
  #   n_starts = opt_control$n_starts
  #   start_points <- generate_sobol_points(n_starts, bounds_lower, bounds_upper, opt_control$seed)
  # } else {
    n_starts = nrow(start_points)
    # start_points <- as.matrix(start_points)
  # }
  
  if (is.null(opt_control$iter_nstart)) {
    opt_control$iter_nstart <- opt_control$n_starts
  }
  
  # Prepare function to process each starting point
  process_start_point <- function(start_point) {
    SMCO_single_boost(f = f, bounds_lower, bounds_upper, start_point = start_point,
                      bounds_buffer = opt_control$bounds_buffer, buffer_rand = opt_control$buffer_rand, 
                      iter_max = opt_control$iter_max, iter_nstart = opt_control$iter_nstart,
                      iter_boost = opt_control$iter_boost, 
                      tol_conv = opt_control$tol_conv, refine_search = opt_control$refine_search, 
                      refine_ratio = opt_control$refine_ratio, partial_option = opt_control$partial_option,
                      use_runmax = opt_control$use_runmax)
  }
  
  if (opt_control$use_parallel) {
    # Parallel execution with optimized cluster setup
    n_cores <- min(parallel::detectCores() - 1, opt_control$n_starts)
    cat(n_cores, "\n")
    cl <- parallel::makeCluster(n_cores)
    
    # Export necessary functions and objects
    parallel::clusterExport(cl, c("SMCO_single_boost", "SMCO_single_refine", "SMCO_single", 
                                  "check_bounds", "compute_partial_signs", "Rast"))
    
    # Run optimization in parallel
    results <- parallel::parApply(cl, start_points, 1, process_start_point)
    
    parallel::stopCluster(cl)
    
  } else {
    # Sequential execution
    results <- apply(start_points, 1, process_start_point)
  }
  
  # Find best result
  f_values <- sapply(results, function(x) x$f_optimal)
  best_idx <- which.max(f_values)
  best_result <- results[[best_idx]]
  
  # Compile summary statistics
  all_endpoints <- t(sapply(results, function(x) x$x_optimal))
  all_iterations <- sapply(results, function(x) x$iterations)
  
  return(list(
    best_result = best_result,
    all_results = results,
    opt_control = opt_control,
    summary = list(
      n_starts = n_starts,
      #success_rate = sum(abs(f_values - max(f_values)) < opt_control$tol_conv) / opt_control$n_starts,
      mean_iterations = mean(all_iterations),
      std_values = sd(f_values),
      endpoints = all_endpoints,
      values = f_values
    )
    ))
}

# Load the compiler package
library(compiler)

# Compile the core computation functions
generate_sobol_points <- cmpfun(generate_sobol_points)
compute_partial_signs <- cmpfun(compute_partial_signs)
check_bounds <- cmpfun(check_bounds)
SMCO_single <- cmpfun(SMCO_single)
SMCO_single_refine <- cmpfun(SMCO_single_refine)
SMCO_single_boost <- cmpfun(SMCO_single_boost)
SMCO_multi <- cmpfun(SMCO_multi)


# Example usage
# run_example <- function() 
if (F) {
  # Define the objective function
  f <- function(x) {
    return(-sum(x^2))
  }
  
  # Define the bounds
  bounds_lower <- c(-5, -5)
  bounds_upper <- c(5, 5)
  
  # Run the optimization
  time_start = Sys.time()
  result <- SMCO_multi(f, bounds_lower, bounds_upper,  
                       opt_control = list(
                         n_starts = 100,
                         iter_max = 200,
                         iter_nstart = 1,
                         iter_boost = 0,
                         bounds_buffer = 0.05,
                         buffer_rand = F,
                         tol_conv = 1e-8,
                         refine_search = T,
                         refine_ratio = 0.55,
                         partial_option = "center",
                         use_runmax = T,
                         use_parallel = F,         
                         seed = 123))
  time_end = Sys.time()
  time_run = difftime(time_end, time_start, units = "secs")
  
  # Print the best result
  print(time_run)
  print(result$best_result$f_optimal)
  
}
