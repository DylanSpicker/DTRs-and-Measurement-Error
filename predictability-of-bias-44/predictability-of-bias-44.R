# Clear the current data; set the correct working directory
rm(list = ls())
rm(list = ls())
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("../Regression Calibration/rCalibration.R")

# Seed and set paraeters
set.seed(314159265)
ns <- c(25000) # c(100, 500, 1000, 10000)
runs <- 10000
num_p_n <- 8# How many trials are run per n
rho <- 0.5

# Initialize counter and final matrix
nc <- 0
results <- matrix(nrow=runs, ncol=num_p_n*length(ns))

for (n in ns) {
  cat("Working on ", n, ".\n")
  for (ii in 1:runs) {
    if (ii %% 1000 == 0) { cat("\t", ii, "; ")}

    # Situation 1
    X <- rnorm(n, sd=sqrt(3))
    W1 <- X + rnorm(n, sd=sqrt(9))
    W2 <- X + rnorm(n, sd=sqrt(9))
    W2[sample(1:n, 0.5*n, replace = F)] <- NA
    A <- rbinom(n,1,1/(1+exp(W1)))
    Y <- 3 + 0.5*X + A*(0.5*X-2) + rnorm(n)
    Wr <- RC(cbind(W1, W2))
    tr.W <- fitted(glm(A~Wr, family = 'binomial'))
    tr.W1 <- fitted(glm(A~W1, family = 'binomial'))
    wt.W <- abs(A - tr.W)
    wt.W1 <- abs(A - tr.W1)
    coeffs.W <- lm(Y ~ A + A:Wr, weights = wt.W)$coeff
    coeffs.W1 <- lm(Y ~ A + A:W1, weights = wt.W1)$coeff
      
    results[ii, (nc*num_p_n+1):(nc*num_p_n+4)] <- c(coeffs.W[2], coeffs.W[3], coeffs.W1[2], coeffs.W1[3])
    
    # Situation 2
    # X <- rnorm(n, sd=sqrt(3))
    # W1 <- X + rnorm(n, sd=sqrt(9))
    # W2 <- X + rnorm(n, sd=sqrt(9))
    # A <- rbinom(n,1,1/(1+exp(W1)))
    # Y <- 3 - X + A*(X - 2) + rnorm(n)
    # Wr <- RC(cbind(W1, W2))
    # tr.W <- fitted(glm(A~Wr, family = 'binomial'))
    # tr.W1 <- fitted(glm(A~W1, family = 'binomial'))
    # wt.W <- abs(A - tr.W)
    # wt.W1 <- abs(A - tr.W1)
    # coeffs.W <- lm(Y ~ A + A:Wr, weights = wt.W)$coeff
    # coeffs.W1 <- lm(Y ~ A + A:W1, weights = wt.W1)$coeff
    # 
    # results[ii, (nc*num_p_n+7):(nc*num_p_n+12)] <- c(coeffs.W[2], coeffs.W[3], coeffs.W1[2], coeffs.W1[3])
    
    # Situation 3
    # X <- rnorm(n, sd=sqrt(9))
    # W1 <- X + rnorm(n, sd=sqrt(5))
    # W2 <- X + rnorm(n, sd=sqrt(5))
    # A <- rbinom(n,1,1/(1+exp(W1)))
    # Y <- 1 + 15*X + exp(X) + A*(1 - 15*X) + rnorm(n)
    # Wr <- RC(cbind(W1, W2))
    # tr.W <- fitted(glm(A~Wr, family = 'binomial'))
    # tr.W1 <- fitted(glm(A~W1, family = 'binomial'))
    # wt.W <- abs(A - tr.W)
    # wt.W1 <- abs(A - tr.W1)
    # coeffs.W <- lm(Y ~ A*Wr, weights = wt.W)$coeff
    # coeffs.W1 <- lm(Y ~ A*W1, weights = wt.W1)$coeff
    # 
    # 
    # results[ii, (nc*num_p_n+13):(nc*num_p_n+18)] <- c(coeffs.W[2], coeffs.W[4], coeffs.W1[2], coeffs.W1[4])
    
    # Situation 4
    X <- rnorm(n, mean=3, sd=sqrt(3))
    W1 <- X + rnorm(n, sd=sqrt(9))
    W2 <- X + rnorm(n, sd=sqrt(9))
    W2[sample(1:n, 0.5*n, replace = F)] <- NA
    A <- rbinom(n,1,1/(1+exp(-10*W1)))
    Y <- 1 + X + exp(X) + A*(1 + X) + rnorm(n)
    Wr <- RC(cbind(W1, W2))
    tr.W <- fitted(glm(A~Wr, family = 'binomial'))
    tr.W1 <- fitted(glm(A~W1, family = 'binomial'))
    wt.W <- abs(A - tr.W)
    wt.W1 <- abs(A - tr.W1)
    coeffs.W <- lm(Y ~ A*Wr, weights = wt.W)$coeff
    coeffs.W1 <- lm(Y ~ A*W1, weights = wt.W1)$coeff
    
    
    results[ii, (nc*num_p_n+5):(nc*num_p_n+8)] <- c(coeffs.W[2], coeffs.W[4], coeffs.W1[2], coeffs.W1[4])
    
  }
  save.image()
  nc <- nc+1
  cat("\n")

}
save.image()

# Generate Plots (WIP - will be moved eventually)

lower<- round((runs/2) - 1.96*sqrt(runs)/2)
upper<- round(1 + (runs/2) + 1.96*sqrt(runs)/2)

sorted_results <- results
sorted_results[,2] <- sort(results[,2])
sorted_results[,4] <- sort(results[,4])
sorted_results[,6] <- sort(results[,6])
sorted_results[,8] <- sort(results[,8])


medians <- c(median(sorted_results[,2]), 
             median(sorted_results[,4]),
             median(sorted_results[,6]),
             median(sorted_results[,8]))

CI.up <- c(sorted_results[upper,2],
           sorted_results[upper,4],
           sorted_results[upper,6],
           sorted_results[upper,8])
CI.dn <- c(sorted_results[lower,2],
           sorted_results[lower,4],
           sorted_results[lower,6],
           sorted_results[lower,8])

par(mfrow=c(1,2), mar=c(5.1, 4.1, 1.1, 2.1), mgp=c(3, 1, 0), las=0)
x_pos <- c(1,2,1,2)
plot(medians[1:2]~x_pos[1:2], ylim=c(min(CI.dn[1:2]),max(CI.up[1:2])),
     xaxt='n',ylab=expression(paste(psi[1], ' Estimates')),  xlim=c(0.5,2.5), xlab='Estimator')
axis(1, at=x_pos[1:2], labels=c("RC", "Naive"))
arrows(x_pos[1:2],CI.dn[1:2],x_pos[1:2],CI.up[1:2],code=3,length=0.2,angle=90)
abline(h=0.5, lty=1, lwd=2)

plot(medians[3:4]~x_pos[3:4], ylim=c(min(min(CI.dn[3:4]),1),max(CI.up[3:4])),
     xaxt='n', ylab=expression(paste(psi[1], ' Estimates')), xlim=c(0.5,2.5), xlab='Estimator')
axis(1, at=x_pos[3:4], labels=c("RC", "Naive"))
arrows(x_pos[3:4],CI.dn[3:4],x_pos[3:4],CI.up[3:4],code=3,length=0.2,angle=90)
abline(h=1, lty=1, lwd=2)
