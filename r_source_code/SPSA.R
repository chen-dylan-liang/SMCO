# Simultaneous Perturbation Stochastic Approximation (SPSA)
# Implemented according to Spall (1992), Guo and Fu (2022)
# Parameters:
# f: objective function (takes vector input, returns scalar value)
# startpoint: initial parameter vector
# A: step decay constant
# a: step scaling factor (ak = a / (k+1)^alpha)
# c: disturbance scaling factor (ck = c / (k+1)^gamma)
# alpha: step decay rate (default 0.602)
# gamma: perturbation decay rate (default 0.101)
# max_iter: maximum number of iterations

spsa <- function(f, startpoint, bounds_lower, bounds_upper, A, a, c, alpha = 0.602, gamma = 0.101, max_iter = 1e+3, tol = 1e-7) {
  # dimension
  d <- length(startpoint)
  lower_bounds <- bounds_lower
  upper_bounds <- bounds_upper
  loss_history <- numeric(max_iter)
  # initial starting point
  x_current <- startpoint
  for (k in 1:max_iter) {
    # compute ak and ck
    ak <- a/(A + k)^alpha
    ck <- c/k^gamma
    # Generate simultaneous perturbation vectors (±1 Bernoulli distribution)
    delta <- 2*rbinom(d, 1, 0.5) - 1
    # compute x_plus and x_minus
    x_plus <- x_current + ck*delta
    x_minus <- x_current - ck*delta
    # compute f value
    f_plus <- f(x_plus)
    f_minus <- f(x_minus)
    # compute grad
    g_hat <- (f_plus - f_minus)/(2*ck*delta)
    # update x
    x_current <- x_current - ak * g_hat
    # control
    x_current[x_current > upper_bounds] <- upper_bounds[x_current > upper_bounds]
    x_current[x_current < lower_bounds] <- lower_bounds[x_current < lower_bounds]
    # record f value
    loss_history[k] <- f(x_current)
    if(k > 1 && abs(loss_history[k]-loss_history[k-1]) < tol) {break}
  }
  return(list(f_optimal = loss_history[k], x_optimal = x_current, loss_history = loss_history))
}

# Example

if (FALSE) {
  startpoint <- c(2, -2)
  a <- 0.10
  c <- 1e-3
  bounds <- list(c(-3, 3))
  f <- function(x) {
    return(-x[1]^2 - 2*x[2]^2)
  }
  
  result <- spsa(f, 
                 startpoint = startpoint, 
                 A = 50,
                 bounds, 
                 a = a, 
                 c = c,
                 alpha = 0.602,
                 gamma = 0.101,
                 max_iter = 1e+3, 
                 tol = 1e-7)
  
  print(paste("result:", paste(round(result$x_optimal, 4), collapse = ", ")))
}
