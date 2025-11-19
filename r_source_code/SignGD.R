# Sign Gradient Descent Optimizer Implementation

#' Calculate numerical gradient
#' @param f Function to optimize
#' @param x Current point
#' @param eps Small value for numerical derivative
#' @return Gradient vector at point x
numerical_gradient <- function(f, x, eps = 1e-8) {
  n <- length(x)
  grad <- numeric(n)
  
  for (i in 1:n) {
    h <- rep(0, n)
    h[i] <- eps
    # Central difference formula
    grad[i] <- (f(x + h) - f(x - h)) / (2 * eps)
  }
  
  return(grad)
}

#' Sign Gradient Descent optimizer
#' @param f Function to optimize
#' @param x0 Initial point
#' @param max_iter Maximum number of iterations
#' @param tol Tolerance for convergence
#' @param learning_rate Learning rate (step size)
#' @param decay_rate Learning rate decay
#' @param verbose Whether to print progress
#' @return List containing optimization results
signGD <- function(f, x0, bounds_lower, bounds_upper,
                                max_iter = 1000, 
                                tol = 1e-6, 
                                learning_rate = 0.1,
                                decay_rate = 0.995,
                                verbose = FALSE) {
  
  # Initialize
  x <- x0
  n <- length(x0)
  history <- matrix(NA, nrow = max_iter, ncol = n + 1)  # Store x and f(x)
  
  # Initialize best solution tracking
  best_x <- x
  best_value <- f(x)
  no_improve_count <- 0
  
  for (iter in 1:max_iter) {
    # Store current state
    history[iter, 1:n] <- x
    current_value <- f(x)
    history[iter, n + 1] <- current_value
    
    # Calculate gradient
    grad <- numerical_gradient(f, x)
    
    # Get sign of gradient
    grad_sign <- sign(grad)
    
    # Update position using only the sign
    x_new <- x - learning_rate * grad_sign
    
    x_new <- pmax(pmin(x_new, bounds_upper), bounds_lower)
    
    # Evaluate new position
    new_value <- f(x_new)
    
    # Update best solution if improved
    if (new_value < best_value) {
      best_x <- x_new
      best_value <- new_value
      no_improve_count <- 0
    } else {
      no_improve_count <- no_improve_count + 1
    }
    
    # Decay learning rate
    learning_rate <- learning_rate * decay_rate
    
    # Check for convergence
    if (sqrt(sum((x_new - x)^2)) < tol) {
      if (verbose) {
        cat("Converged after", iter, "iterations\n")
      }
      break
    }
    
    # Additional stopping criterion: if no improvement for many steps
    if (no_improve_count > 50) {
      if (verbose) {
        cat("Stopped after", iter, "iterations due to no improvement\n")
      }
      break
    }
    
    # Update for next iteration
    x <- x_new
    
    # Print progress
    if (verbose && iter %% 100 == 0) {
      cat("Iteration:", iter, 
          "Value:", current_value, 
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

# Example usage with Rosenbrock function
rosenbrock <- function(x) {
  return((1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2)
}

# Test the optimizer
test_sign_gd <- function() {
  # Starting point
  x0 <- c(-1, 1)
  bounds_lower = rep(-5, 2)
  bounds_upper = rep(10,2)
  # Run optimization
  result <- signGD(rosenbrock, x0, bounds_lower, bounds_upper,
                                max_iter = 5000,
                                learning_rate = 0.01,
                                verbose = TRUE)
  
  # Print results
  cat("\nOptimization Results:\n")
  cat("Solution:", result$solution, "\n")
  cat("Value:", result$value, "\n")
  cat("Iterations:", result$iterations, "\n")
  
  # Plot optimization path if running interactively
  if (interactive()) {
    library(ggplot2)
    
    # Create data frame from history
    history_df <- data.frame(
      x1 = result$history[, 1],
      x2 = result$history[, 2],
      value = result$history[, 3]
    )
    
    # Create plot
    p <- ggplot(history_df, aes(x = x1, y = x2)) +
      geom_path(aes(color = value), size = 1) +
      geom_point(data = data.frame(x1 = result$solution[1], 
                                  x2 = result$solution[2]),
                 color = "red", size = 3) +
      scale_color_viridis_c(name = "Function\nValue") +
      labs(title = "Sign Gradient Descent Optimization Path",
           x = "x1",
           y = "x2") +
      theme_minimal()
    
    print(p)
  }
  
  return(result)
}

# Run test if desired
# test_result <- test_sign_gd()