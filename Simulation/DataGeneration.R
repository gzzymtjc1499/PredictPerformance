current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

library(MASS)

dir.create("SIMU_DATA")

n_simulations <- 100 # Number of Simulations
obs_values <- seq(100, 2000, by = 100)

for (i in 1:length(obs_values)) {
  nobs <- obs_values[i]
  dir.create(paste0("SIMU_DATA/", nobs))
  
  for (j in 1:n_simulations) {
    cat('i=',nobs,'j=',j,sep='\n')
    
    loading1 <- c(0, rep(sqrt(.6), 10), rep(0, 9))
    loading2 <- c(sqrt(.6), rep(0, 10), rep(sqrt(.6), 9))
    L <- cbind(loading1, loading2)
    
    PSY <- diag(1 - rowSums(L ^ 2))
    SIGMA <- L %*% t(L) + PSY
    
    X <- mvrnorm(
      n = nobs,
      mu = rep(0, 20),
      Sigma = SIGMA,
      empirical = TRUE
    )
    
    # Generate Y
    
    # define regression coefficients
    beta0 <- 0
    beta1 <- 1
    beta2 <- -1
    
    Score <- X %*% L %*% solve(t(L) %*% L) # obtain factor scores of X
    Y <- rep(beta0, nobs) + beta1 * Score[, 1] + beta2 * Score[, 2] + rnorm(n =
                                                                              nobs, mean = 0, sd = 0.01)
    
    out <- cbind(X, Y)
    
    save(out, file = paste0("SIMU_DATA/",nobs,"/data_",nobs,"_",j, ".RData"))
  }
}

