# Supplemental Test Functions

# "Optimization via Strategic Law of Large Numbers"
# By: Xiaohong Chen, Zengjing Chen, Wayne Yuan Gao, Xiaodong Yan, and Guodong Zhang
# Date: July 15, 2025
# GitHub Repository Maintained by: Wayne Yuan Gao

# Coding of the following test functions are retrieved from:
# Sonja Surjanovic and Derek Bingham, (2013), Virtual Library of Simulation Experiments
# https://www.sfu.ca/~ssurjano/optimization.html
# Ackley, Bukin No.6, Dixon-Price, Dropwave, EggHolder, Griewank, McCormick, Michalewicz, 
# Rastrigin, Rosenbrock, Shubert, Six-Hump, Zakharov

# Codoing of the following test functions provided by the authors:
# Damavanndi, (Minus) Squared Norm, Mishra No.6, Qing


########################
# Squared Norm Function

norm2x <- function(z) {sum(z^2)}

config_norm2x <- function(d) {
  config <- list(
    name = "Norm Squared",
    dim = d,
    bounds_lower = rep(-1,d),
    bounds_upper = rep(1, d),
    f = norm2x,
    f_min = 0,
    f_max = norm2x(rep(1, d))
  )
  return(config)
} 

norm2x <- function(z) {sum(z^2)}

config_norm2x_minus <- function(d) {
  config <- list(
    name = "Minus Norm Squared",
    dim = d,
    bounds_lower = rep(-1,d),
    bounds_upper = rep(1, d),
    f = function(x) {-norm2x(x)},
    f_max = 0,
    f_min = -norm2x(rep(1, d))
  )
  return(config)
} 

dropwave <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  frac1 <- 1 + cos(12*sqrt(x1^2+x2^2))
  frac2 <- 0.5*(x1^2+x2^2) + 2
  
  y <- -frac1/frac2
  return(y)
}

config_drop <- list(
  name = "Dropwave",
  dim = 2,
  bounds_lower = rep(-5.12,2),
  bounds_upper = rep(5.12, 2),
  f = dropwave,
  f_max = NA,
  f_min = -1
)

# Shubert 

shubert <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  ii <- c(1:5)
  
  sum1 <- sum(ii * cos((ii+1)*x1+ii))
  sum2 <- sum(ii * cos((ii+1)*x2+ii))
  
  y <- sum1 * sum2
  return(y)
}

config_shubert <- list(
  name = "Shubert",
  dim = 2,
  bounds_lower = rep(-5.12,2),
  bounds_upper = rep(5.12, 2),
  f = shubert,
  f_max = NA,
  f_min = -186.7309
)

# McCormick

mccorm <- function(xx)
{

  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- sin(x1 + x2)
  term2 <-(x1 - x2)^2
  term3 <- -1.5*x1
  term4 <- 2.5*x2
  
  y <- term1 + term2 + term3 + term4 + 1
  return(y)
}

config_mccorm <- list(
  name = "McCormick",
  dim = 2,
  bounds_lower = c(-1.5,-3),
  bounds_upper = c(4, 4),
  f = mccorm,
  f_max = NA,
  f_min = -1.9133
)

# Zakharov

zakharov <- function(xx)
{
  ii <- c(1:length(xx))
  sum1 <- sum(xx^2)
  sum2 <- sum(0.5*ii*xx)
  
  y <- sum1 + sum2^2 + sum2^4
  return(y)
}

config_zak <- function(d) {
  config <- list(
    name = "Zakharov",
    dim = d,
    bounds_lower = rep(-5,d),
    bounds_upper = rep(10,d),
    f = zakharov,
    f_max = NA,
    f_min = 0
  )
  return(config)
}

# Dixon-Price

dixon <- function(xx)
{
  x1 <- xx[1]
  d <- length(xx)
  term1 <- (x1-1)^2
  
  ii <- c(2:d)
  xi <- xx[2:d]
  xold <- xx[1:(d-1)]
  sum <- sum(ii * (2*xi^2 - xold)^2)
  
  y <- term1 + sum
  return(y)
}

config_dixon <- function(d) {
  config <- list(
    name = "DixonPrice",
    dim = d,
    bounds_lower = rep(-10,d),
    bounds_upper = rep(10,d),
    f = dixon,
    f_max = NA,
    f_min = 0
  )
  return(config)
}


# Easom  
easom <- function(xx)
{
 
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- -cos(x1)*cos(x2)
  fact2 <- exp(-(x1-pi)^2-(x2-pi)^2)
  
  y <- fact1*fact2
  return(y)
}

config_easom <- list(
    name = "Easom",
    dim = 2,
    bounds_lower = rep(-100,2),
    bounds_upper = rep(100,2),
    f = easom,
    f_max = NA,
    f_min = -1
  )

# Rastrigin Function

Rast <- function(z) {
  d <- length(z)
  10*d + sum(z^2-10*cos(2*pi*z))
}

config_rast <- function(d) {
  x_max = rep(4.52299366,d)
  nameconfig <- paste0("Rastrigin",  d, "d")
  config <- list(
    name = nameconfig,
    dim = d,
    bounds_lower = rep(-5.12,d),
    bounds_upper = rep(5.12, d),
    f = Rast,
    x_min = rep(0, d),
    f_min = 0,
    x_max = x_max,
    f_max = Rast(x_max)
  )
  return(config)
}

# Bukin No. 6
bukin6 <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- 100 * sqrt(abs(x2 - 0.01*x1^2))
  term2 <- 0.01 * abs(x1+10)
  
  y <- term1 + term2
  return(y)
}

config_bukin6 <- list(
  name = "Bukin6",
  dim = 2,
  bounds_lower = c(-15,-5),
  bounds_upper = c(-3, 3),
  f = bukin6,
  f_max = NA,
  f_min = 0
)


# Damavandi

damava <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  term1 <- 1 - abs(sin(pi*(x1-2)) *  sin(pi*(x2-2)))^5
  term2 <- 2 + (x1-7)^2 + 2 * (x2-7)^2
  return(term1*term2)
}


config_dama <- list(
  name = "Damavandi",
  dim = 2,
  bounds_lower = rep(0,2),
  bounds_upper = rep(14,2),
  f = damava,
  f_max = NA,
  f_min = 0
)

# Mishra No 6

mish06 <- function(x) {
  x1 = x[1]
  x2 = x[2]
  term1 = (sin((cos(x1) + cos(x2))^2))^2
  term2 = (cos((sin(x1) + sin(x2))^2))^2
  term3 = 0.1 * ((x1 - 1)^2 + (x2-1)^2)
  final = - log((term1 - term2 + x1)^2) + term3
  return(final)
}

config_mish06 <- list(
  name = "Mishra06",
  dim = 2,
  bounds_lower = rep(-10,2),
  bounds_upper = rep(10,2),
  f = mish06,
  f_max = NA,
  f_min = -2.2839498
)

# Qing

qing <- function(xx) {
  d <- length(xx)
  ind_xx <- c(1:d)
  term_xx <- (xx^2-ind_xx)^2
  y <- sum(term_xx)
}

config_qing <- function(d) {
  config <- list(
    name = "Qing",
    dim = d,
    bounds_lower = rep(-500,d),
    bounds_upper = rep(500,d),
    f = qing,
    f_max = NA,
    f_min = 0
  )
}

# Ackley Function
  ackley <- function(xx, a=20, b=0.2, c=2*pi)
    {
      d <- length(xx)
  
      sum1 <- sum(xx^2)
      sum2 <- sum(cos(c*xx))
  
      term1 <- -a * exp(-b*sqrt(sum1/d))
      term2 <- -exp(sum2/d)
      
      y <- term1 + term2 + a + exp(1)
      return(y)
    }

  config_ackley <- function(d){
    config <- list(
    name = paste0("Ackley", d, "d"),
    dim = d,
    bounds_lower = rep(-32.768,d),
    bounds_upper = rep(32.768, d),
    # bounds = rep(list(c(-32.768, 32.768)), d),
    # bounds = rep(list(c(-2, 2)), d),
    # f_max = 21.71828, # From ChatGPT on 20-dim
    # f_max = 22.3504, # From ChatGPT on 30-dim
    f = ackley,
    f_min = 0,
    x_min = rep(0, d),
    f_max = NA)
    return(config)
}

# Griewank Function
  
  griewank <- function(xx)
  {
    ii <- c(1:length(xx))
    sum <- sum(xx^2)/4000
    prod <- prod(cos(xx/sqrt(ii)))
    
    y <- sum - prod + 1
    return(y)
  }
  
  config_griewank <- function(d){
    config <- list(
      name = paste0("Griewank", d, "d"),
      dim = d,
      f = griewank,
      #bounds = rep(list(c(-600, 600)), d),
      bounds_lower = rep(-600,d),
      bounds_upper = rep(600, d),
      # f_max = 902, # From Claude
      # bounds = rep(list(c(-60, 60)), d),
      # f_max = 10, # From Claude
      # bounds = rep(list(c(-6, 6)), d),
      # f_max = 1.774, # From Claude      
      f_min = 0,
      f_max = NA
    )
    return(config)
  }

schwef <- function(xx)
  {
    
    d <- length(xx)
    
    sum <- sum(xx*sin(sqrt(abs(xx))))
    
    y <- 418.9829*d - sum
    return(y)
}

config_schwef <- function(d) {
  config <- list(
    name = paste0("Schwefel", d, "d"),
    dim = d,
    bounds_lower = rep(-500,d),
    bounds_upper = rep(500, d),
    f = schwef,
    f_min = 0,
    x_min = rep(420.9687, d),
    f_max = NA)
  return(config)
}
  

rosen <- function(xx)
{
  d <- length(xx)
  xi <- xx[1:(d-1)]
  xnext <- xx[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  
  y <- sum
  return(y)
}

config_rosen <- function(d) {
  config <- list(
    name = paste0("Rosenbrock", d, "d"),
    dim = d,
    bounds_lower = rep(-5,d),
    bounds_upper = rep(10, d),
    f = rosen,
    f_min = 0,
    x_min = rep(1, d),
    f_max = NA)
  return(config)
}

    
# Michalewicz Function

  michal <- function(xx, m=10)
    {
     ii <- c(1:length(xx))
     sum <- sum(sin(xx) * (sin(ii*xx^2/pi))^(2*m))
  
     y <- -sum
      return(y)
    }

  michal_dim = c(2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
  michal_min = c(
            -1.8013034,
            -2.7603947,
            -3.6988571,
            -4.6876582,
            -5.6876582,
            -6.6808853,
            -7.6637574,
            -8.6601517,
            -9.6601517,
            -14.6464002,
            -19.6370136,
            -24.6331947,
            -29.6308839,
            -34.6288550,
            -39.6267489,
            -44.6256251,
            -49.6248323,
            -54.6240533,
            -59.6231462,
            -64.6226167,
            -69.6222202,
            -74.6218112)
   # Reference: Tables 3 & 6 in https://arxiv.org/abs/2003.09867


  config_michal <- function(d){
    
    if (d %in% michal_dim) {
      match_dim <- match(d, michal_dim)
      f_min <- michal_min[match_dim]
    } else{
      f_min = NA
    }
    # if ( match_dim == 0){
    #   print("Invalid dimension for Michalewicz function")
    #   return(NULL)
    # }
    config <- list(
      name = paste0("Michalewicz", d, "d"),
      dim = d,
      bounds_lower = rep(0,d),
      bounds_upper = rep(pi, d),
      #bounds = rep(list(c(0, pi)), d),
      f = michal,
      f_min = f_min,
      f_max = 0
    )
    return(config)
  }


# Six-Hump Camel Function

  camel6 <- function(xx)
  {
    x1 <- xx[1]
    x2 <- xx[2]
    
    term1 <- (4-2.1*x1^2+(x1^4)/3) * x1^2
    term2 <- x1*x2
    term3 <- (-4+4*x2^2) * x2^2
    
    y <- term1 + term2 + term3
    return(y)
  }

  config_camel6 <- list(
      name = "SixHumpCamel",
      bounds_lower = c(-3,-2),
      bounds_upper = c(3,2),
      #bounds = list(c(-3, 3), c(-2,2)),
      dim = 2,
      f = camel6,
      f_min = -1.0316,
      f_max = NA
    )


  # Egg Holder Function
  
  egg <- function(xx)
  {
    x1 <- xx[1]
    x2 <- xx[2]
    
    term1 <- -(x2+47) * sin(sqrt(abs(x2+x1/2+47)))
    term2 <- -x1 * sin(sqrt(abs(x1-(x2+47))))
    
    y <- term1 + term2
    return(y)
  }
  
  config_egg <- list(
    name = "EggHolder",
    dim = 2,
    bounds_lower = rep(-512,2),
    bounds_upper = rep(512,2),
    f = egg,
    f_min = -959.6407,
    x_min = c(512, 404.2319),
    f_max = NA
  )
  
  
# Allocate configuration based on name and dimension

  assign_config <- function(name, d, opt_config){
    if (name == "SqNorm") {
      test_config <- config_norm2x(d)
    } else if (name == "-SqNorm") {
      test_config <- config_norm2x_minus(d)
    } else if (name == "Ackley"){
      test_config <- config_ackley(d)
    } else if (name == "Six-Hump Camel"){
      test_config <- config_camel6
    } else if (name == "Rastrigin"){
      test_config <- config_rast(d)
    } else if (name == "Michalewicz"){
      test_config <- config_michal(d)
    } else if (name == "Griewank"){
      test_config <- config_griewank(d)
    } else if (name == "EggHolder") {
      test_config <- config_egg
    } else if (name == "Rosenbrock") {
      test_config <- config_rosen(d)
    } else if (name == "Dropwave") {
      test_config <- config_drop
    } else if (name == "Shubert") {
      test_config <- config_shubert
    } else if (name == "McCormick") {
      test_config <- config_mccorm
    } else if (name == "DixonPrice") {
      test_config <- config_dixon(d)
    } else if (name == "Zakharov") {
      test_config <- config_zak(d)
    } else if (name == "Bukin6") {
      test_config <- config_bukin6
    } else if (name == "Qing") {
      test_config <- config_qing(d)
    } else if (name == "Damavandi") {
      test_config <- config_dama
    } else if (name == "DVGlasser02") {
      test_config <- config_DVG02
    } else if (name == "Mishra06") {
      test_config <- config_mish06
    } else if (name == "ANN") {
      source("testfunc_NN1L.R")
      test_config <- config_ANN(opt_config$input_dim, opt_config$hidden_dim, opt_config$n,
                                opt_config$use_truth, opt_config$noise_sigma, opt_config$seed)
    } else {
      test_config <- print("Invalid test function name")
    }
    return(test_config)
  }
