# Clear the current data; set the correct working directory
rm(list = ls())
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("../Regression Calibration/rCalibration.R")

# Seed and set parameters
set.seed(314159265)
ns <- c(100, 500, 1000, 10000)
runs <- 10000
num_p_n <- 8 # How many trials are run per n

# Initialize counter and final matrix
nc <- 0
results <- matrix(nrow=runs, ncol=num_p_n*length(ns))

for (n in ns) {
  cat("Working on ", n, ".\n")
  for (ii in 1:runs) {
    if (ii %% 100 == 0) { cat("\t", ii, "; ")}
    
    X <- rnorm(n)
    W1 <- X + rnorm(n)
    W2 <- X + rnorm(n)
    A <- rbinom(n,1,1/(1+exp(-W1)))
    Y <- 1 + X + A*(1 + X) + rnorm(n)
    
    # NA the relevant W2s
    W2_10 <- W2
    W2_10[sample(c(1:n), n*0.9)] <- NA
    W2_50 <- W2
    W2_50[sample(c(1:n), n*0.5)] <- NA
    
    # Rcal
    Wr_full <- RC(cbind(W1, W2))
    Wr_10 <- RC(cbind(W1, W2_10))
    Wr_50 <- RC(cbind(W1, W2_50))
    
    # Treatment Models
    naive_fit <- fitted(glm(A~W1, family = 'binomial'))
    full_fit <- fitted(glm(A~Wr_full, family = 'binomial'))
    part_10_fit <- fitted(glm(A~Wr_10, family = 'binomial'))
    part_50_fit <- fitted(glm(A~Wr_50, family = 'binomial'))
    
    # Weighted
    wt.n <- abs(A - naive_fit)
    wt.f <- abs(A - full_fit)
    wt.10 <- abs(A - part_10_fit)
    wt.50 <- abs(A - part_50_fit)
    
    # Model Coefficients
    coeff.n <- lm(Y ~ A*W1, weights = wt.n)$coeff
    coeff.f <- lm(Y ~ A*Wr_full, weights = wt.f)$coeff
    coeff.10 <- lm(Y ~ A*Wr_10, weights = wt.10)$coeff
    coeff.50 <- lm(Y ~ A*Wr_50, weights = wt.50)$coeff
    
    results[ii, (nc*num_p_n + 1):((nc+1)*num_p_n)] <- c(coeff.n[2], coeff.n[4],
                                                        coeff.f[2], coeff.f[4],
                                                        coeff.10[2], coeff.10[4],
                                                        coeff.50[2], coeff.50[4])
    
  }
  nc <- nc+1
  
  # Periodically Save
  save.image()
  cat("\n")
}

save.image()