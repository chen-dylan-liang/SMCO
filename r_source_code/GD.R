# Numerical Gradient Descent Optimizer

#' Calculate numerical gradient
#' @param f Function to optimize
#' @param x Current point
#' @param eps Small value for numerical derivative
#' @return Gradient vector at point x
numerical_gradient <- function(f, x, eps = 1e-8) {
  n <- length(x)
  grad <- numeric(n)
  
  for (i in 1:n) {
    # Create perturbation vector
    h <- rep(0, n)
    h[i] <- eps
    
    # Central difference formula
    grad[i] <- (f(x + h) - f(x - h)) / (2 * eps)
  }
  
  return(grad)
}

#' Gradient descent optimizer with adaptive learning rate
#' @param f Function to optimize
#' @param x0 Initial point
#' @param max_iter Maximum number of iterations
#' @param tol Tolerance for convergence
#' @param learning_rate Initial learning rate
#' @param momentum Momentum coefficient
#' @param verbose Whether to print progress
#' @return List containing optimization results
gradient_descent <- function(f, x0,
                           bounds_lower,
                           bounds_upper,
                           max_iter = 1000, 
                           tol = 1e-6, 
                           learning_rate = 0.1,
                           momentum = 0.9,
                           verbose = F) {
  
  # Initialize
  x <- x0
  n <- length(x0)
  velocity <- rep(0, n)
  history <- matrix(NA, nrow = max_iter, ncol = n + 1)  # Store x and f(x)
  
  # Initialize best solution
  best_x <- x
  best_value <- f(x)
  
  # For adaptive learning rate
  previous_value <- f(x)
  
  for (iter in 1:max_iter) {
    # Store current state
    history[iter, 1:n] <- x
    history[iter, n + 1] <- f(x)
    
    # Calculate gradient
    grad <- numerical_gradient(f, x)
    
    # Update velocity with momentum
    velocity <- momentum * velocity - learning_rate * grad
    
    # Update position
    x_new <- x + velocity
    x_new <- pmax(pmin(x_new, bounds_upper), bounds_lower)
    current_value <- f(x_new)
    
    # Adaptive learning rate
    if (current_value > previous_value) {
      # Reduce learning rate if we're not improving
      learning_rate <- learning_rate * 0.5
      velocity <- velocity * 0.5
    } else if (current_value < best_value) {
      # Store best solution
      best_x <- x_new
      best_value <- current_value
    }
    
    # Check convergence
    if (sqrt(sum((x_new - x)^2)) < tol) {
      if (verbose) {
        cat("Converged after", iter, "iterations\n")
      }
      break
    }
    
    # Update for next iteration
    x <- x_new
    previous_value <- current_value
    
    # Print progress
    if (verbose && iter %% 100 == 0) {
      cat("Iteration:", iter, "Value:", current_value, 
          "Learning rate:", learning_rate, "\n")
    }
  }
  
  # Trim history to actual number of iterations
  history <- history[1:iter, , drop = FALSE]
  
  return(list(
    solution = best_x,
    value = best_value,
    iterations = iter,
    history = history,
    learning_rate = learning_rate
  ))
}
