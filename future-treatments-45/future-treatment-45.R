# Clear the current data; set the correct working directory
rm(list = ls())
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("../Regression Calibration/RCibration.R")

# Seed and set parameters
set.seed(314159265)
ns <- c(100, 500, 1000, 10000)
runs <- 10000
sig.u <- sqrt(1) # SD of Error
num_p_n <- 10
psi_1s <- c(0, 1,-2,15)
psi_2s <- c(1,-2,10)

# Initialize counter and final matrix
nc <- 0
results <- matrix(nrow=runs, ncol=num_p_n*length(psi_1s)*length(psi_2s)*length(ns))

for (n in ns) {
  cat("Working on ", n, ".\n")
  for (psi_1 in psi_1s){
    for(psi_2 in psi_2s){
      cat("\t Psi1: ", psi_1, " and Psi2: ", psi_2, "\n\t\t")
      for (ii in 1:runs) {
        if (ii %% 1000 == 0) { cat(ii, "; ")}
        X <- rnorm(n, 0, 1)
        W1 <- X + rnorm(n, 0, sig.u)
        W2 <- X + rnorm(n, 0, sig.u)
        # True Treatments
        A <- sign(psi_1 + psi_2*X)
        
        # Scenario 0 (Naive)
        A_naive <- sign(psi_1 + psi_2*W1)
        results[ii, 10*nc + 1] <- sum(A_naive==A)/n
        
        # Scenario 1 (100% Replication)
        Wr <- RC(cbind(W1, W2))
        A_1 <- sign(psi_1 + psi_2*Wr)
        results[ii, 10*nc + 2] <- sum(A_1==A)/n
        
        # Scenario 2 (50% Replication Rand)
        W2_sp <- W2
        W2_sp[sample(c(1:n), n*0.5)] <- NA
        Wr <- RC(cbind(W1, W2_sp))
        A_1 <- sign(psi_1 + psi_2*Wr)
        results[ii, 10*nc + 3] <- sum(A_1==A)/n
        
        # Scenario 3 (10% Replication Rand)
        W2_sp <- W2
        W2_sp[sample(c(1:n), n*0.9)] <- NA
        Wr <- RC(cbind(W1, W2_sp))
        A_1 <- sign(psi_1 + psi_2*Wr)
        results[ii, 10*nc + 4] <- sum(A_1==A)/n
        
        # Scenario 4 (50% Replication Small)
        W2_sp <- W2
        W2_sp[which(abs(psi_1 + psi_2*W1) > quantile(abs(psi_1 + psi_2*W1),0.5))] <- NA
        Wr <- RC(cbind(W1, W2_sp))
        A_1 <- sign(psi_1 + psi_2*Wr)
        results[ii, 10*nc + 5] <- sum(A_1==A)/n
        
        # Scenario 5 (10% Replication Small)
        W2_sp <- W2
        W2_sp[which(abs(psi_1 + psi_2*W1) > quantile(abs(psi_1 + psi_2*W1),0.1))] <- NA
        Wr <- RC(cbind(W1, W2_sp))
        A_1 <- sign(psi_1 + psi_2*Wr)
        results[ii, 10*nc + 6] <- sum(A_1==A)/n
        
        # Scenario 6 (50% Replication Large)
        W2_sp <- W2
        W2_sp[which(abs(psi_1 + psi_2*W1) < quantile(abs(psi_1 + psi_2*W1),0.5))] <- NA
        Wr <- RC(cbind(W1, W2_sp))
        A_1 <- sign(psi_1 + psi_2*Wr)
        results[ii, 10*nc + 7] <- sum(A_1==A)/n
        
        # Scenario 7 (10% Replication Large)
        W2_sp <- W2
        W2_sp[which(abs(psi_1 + psi_2*W1) < quantile(abs(psi_1 + psi_2*W1),0.9))] <- NA
        Wr <- RC(cbind(W1, W2_sp))
        A_1 <- sign(psi_1 + psi_2*Wr)
        results[ii, 10*nc + 8] <- sum(A_1==A)/n
        
        # Scenario 8 (50% Replication PH)
        W2_sp <- W2
        W2_sp[sample(c(1:n), n*0.5)] <- NA
        r_cal <- RC(cbind(W1, W2_sp))
        
        # Wr <- W1*(r_cal$sig.uu + r_cal$sig.xx)/r_cal$sig.xx
        
        # A_1 <- sign(psi_1 + psi_2*Wr)
        # results[ii, 10*nc + 9] <- sum(A_1==A)/n
        
        # Scenario 9 (10% Replication PH)
        W2_sp <- W2
        W2_sp[sample(c(1:n), n*0.9)] <- NA
        r_cal <- RC(cbind(W1, W2_sp))
        
        # Wr <- W1*(r_cal$sig.uu + r_cal$sig.xx)/r_cal$sig.xx
        
        # A_1 <- sign(psi_1 + psi_2*Wr)
        # results[ii, 10*nc + 10] <- sum(A_1==A)/n
      }
      cat("\n")
      save.image()
      nc <- nc+1
    }
  }
}
save.image()


