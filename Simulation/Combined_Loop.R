#############################################################################
# Preparation
#############################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(lavaan)
library(MASS)
library(psych)
library(dplyr)
library(RGCCA)

source("00_Function_LSAM_PartitionedEstimation.R")
source("00_Function_LSAM_ComputeSummaryStats.R")
source("00_Function_LSAM_RegressionCoefs.R")
source("SIMU_FUNCTION.R")
source("LSLVLASSO.R")

#############################################################################
# Simulation
#############################################################################

# Sample size
sample_sizes <- seq(20000, 20000, 0)
# Number of simulation
max_j <- 100

# Create matrix to store the results
results_out1 <- array(0, dim = c(5, length(sample_sizes), max_j))
results_out2 <- array(0, dim = c(5, length(sample_sizes), max_j))

# Outer loop: different sample sizes
for (n_idx in seq_along(sample_sizes)) {
  N <- sample_sizes[n_idx]
  
  # Inner loop: number of simulations
  for (j in 1:max_j) {
    cat(sprintf("Sample size: %s\n", N))
    cat(sprintf("Iteration: %s\n", j))
    
    # Load the data
    file_path <- sprintf("SIMU_DATA/%d/data_%d_%d.RData", N, N, j)
    if (file.exists(file_path)) {
      load(file_path)  
      
      # Split them into X and Y
      X <- out[, 1:20]
      Y <- out[, 21]
      
      # Simulation, N is the sample size, 2 is the number of factors
      result <- SIMU(X, Y, N, 2)
      out1 <- result$out1
      out2 <- result$out2
      
      # Store the result
      results_out1[, n_idx, j] <- out1
      results_out2[, n_idx, j] <- out2
    } else {
      cat(sprintf("File not found: %s\n", file_path))
    }
  }
}

# Save the result
save(results_out1, results_out2, file = "simulation_results_20000.RData")

# Print
print("Simulation completed.")
print("Results stored in simulation_results.RData.")


