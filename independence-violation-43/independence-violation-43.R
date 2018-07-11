# Clear the current data; set the correct working directory
rm(list = ls())
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("../Regression Calibration/rCalibration.R")


# Seed and set parameters
set.seed(314159265)
ns <- c(100, 500, 1000, 10000)
runs <- 10000
alpha_pairs <- list(c(0, 0.5),c(0, -0.5),c(0, 1),c(0, -1),c(0, 2.5),c(0, -2.5),c(0, 5),
  c(0, -5),c(0, 10),c(0, -10),c(0, 15),c(0, -15), # P(A=1) = 0.5 ends
  c(-2.3, 0.5),c(-2.3, -0.5),c(-2.55, 1),c(-2.55, -1),c(-3.9, 2.5),c(-3.9, -2.5),c(-6.75, 5),
  c(-6.75, -5),c(-13, 10),c(-13, -10),c(-19.25, 15),c(-19.25, -15), # P(A=1) = 0.1 ends
  c(2.3, 0.5),c(2.3, -0.5),c(2.55, 1),c(2.55, -1),c(3.9, 2.5),c(3.9, -2.5),c(6.75, 5),
  c(6.75, -5),c(13, 10),c(13, -10),c(19.25, 15),c(19.25, -15)) # P(A=1) = 0.9 ends
num_p_n <- 4# How many trials are run per n

# Initialize counter and final matrix
nc <- 0
results <- matrix(nrow=runs, ncol=num_p_n*length(alpha_pairs)*length(ns))

for (n in ns) {
  cat("Working on ", n, ".\n")
  for (alpha_pair in alpha_pairs){
    for (ii in 1:runs) {
      if (ii %% 100 == 0) { cat("\t", ii, "; ")}
      
      X <- rnorm(n)
      W1 <- X + rnorm(n)
      W2 <- X + rnorm(n)
      
      alpha_0 <- alpha_pair[1]
      alpha_1 <- alpha_pair[2]
      
      # 
      A <- rbinom(n,1,1/(1+exp(-(alpha_0 + alpha_1*W1))))
      A_x <- rbinom(n,1,1/(1+exp(-(alpha_0 + alpha_1*X))))
      
      Y <- 4 + 2*exp(X) + A*(3 - 2*X) + rnorm(n)
      Y_x <- 4 + 2*exp(X) + A_x*(3 - 2*X) + rnorm(n)
      
      # Rcal
      Wr <- RC(cbind(W1, W2))
      
      # Treatment Models
      tr.X <- fitted(glm(A_x~X, family = 'binomial'))
      tr.W <- fitted(glm(A~Wr, family = 'binomial'))
      
      # Weighted
      wt.X <- abs(A_x - tr.X)
      wt.W <- abs(A - tr.W)
        
      # Model Coefficients
      coeffs.X <- lm(Y_x ~ A_x*X, weights = wt.X)$coeff
      coeffs.W <- lm(Y ~ A*Wr, weights = wt.W)$coeff
      
      results[ii, (nc*num_p_n + 1):((nc+1)*num_p_n)] <- c(coeffs.X[2], coeffs.X[4],
                                                          coeffs.W[2], coeffs.W[4])
      
    }
    save.image()
    nc <- nc+1
    cat("\n")
  }
}
save.image()