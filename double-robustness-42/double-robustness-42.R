# Clear the current data; set the correct working directory
rm(list = ls())
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("../Regression Calibration/rCalibration.R")

# Seed and set parameters
set.seed(314159265)
ns <- c(100, 500, 1000, 10000)
runs <- 10000
num_p_n <- 8# How many trials are run per n

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
    A <- rbinom(n,1,1/(1+exp(-exp(W1))))
    Y <- 1 - X + exp(X) + A*(2 + 2*X) + rnorm(n)
    
    # NA the relevant W2s
    W2[sample(c(1:n), n*0.8)] <- NA
    
    # Rcal
    Wr <- RC(cbind(W1, W2))
    
    # Treatment Models
    tr.corr <- fitted(glm(A~exp(Wr), family = 'binomial'))
    tr.incorr <- fitted(glm(A~Wr, family = 'binomial'))
    
    # Weighted
    wt.corr <- abs(A - tr.corr)
    wt.incorr <- abs(A - tr.incorr)
    
    # Model Coefficients
    coeff.an1 <- lm(Y ~ A*Wr, weights = wt.incorr)$coeff
    coeff.an2 <- lm(Y ~ A*Wr, weights = wt.corr)$coeff
    coeff.an3 <- lm(Y ~ A*Wr + exp(Wr), weights = wt.incorr)$coeff
    coeff.an4 <- lm(Y ~ A*Wr + exp(Wr), weights = wt.corr)$coeff
    
    results[ii, (nc*num_p_n + 1):((nc+1)*num_p_n)] <- c(coeff.an1[2], coeff.an1[4],
                                                        coeff.an2[2], coeff.an2[4],
                                                        coeff.an3[2], coeff.an3[5],
                                                        coeff.an4[2], coeff.an4[5])
    
  }
  save.image()
  nc <- nc+1
  cat("\n")
}
save.image()