# Load the data
# simulation_results_0 is the simulation with n. of observations from 50 - 200
# simulation_results_1: n. of observations from 100 - 2000
# simulation_restuls_20000: n. of observations = 20000

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Only load the data that you are interested in
load("simulation_results_0.RData")
load("simulation_results_1.RData")
load("simulation_results_20000.RData")

# Part 0. ANOVA in total
library(reshape2)
library(dplyr)
library(lme4)

methods <- c("SumScore", "SEM", "SAM", "RGCCA", "LS-LV")
observations <- seq(50, 200, by = 10)
repeats <- 1:100

# Covert the data type
results_long <- melt(results_out1, varnames = c("Method", "Observations", "Repeat"), value.name = "Error")

# Name the method and n. of observations
results_long$Method <- factor(results_long$Method, labels = methods)
results_long$Observations <- factor(results_long$Observations, levels = 1:16, labels = observations)

# str(results_long)

model <- lmer(Error ~ Method + Observations + (1 | Repeat), data = results_long)
summary(model)
anova(model)

aov_model <- aov(Error ~ Method * Observations + Error(Repeat), data = results_long)
summary(aov_model)



# Part 1. Relationship between n. of observations and error
# This is for simulation_results_0
par(mfrow=c(1,1))
boxplot(t(results_out2[5,,]), names=seq(50, 200, by = 10), main = "Error2 of LS-LV Method", xlab = "N. of observations", ylab = "Error", outline = FALSE)

# Part 2. Draw the boxplot
# Graph Function
plot_boxplot <- function(data, title) {
  df <- as.data.frame(t(data))
  colnames(df) <- c("SumScore", "SEM", "SAM", "RGCCA", "LS-LV")
  boxplot(df, main = title, ylab = "lg(error)", xlab = "Method", outline = FALSE)
}

# Select a target sample size
# This is for simulation_results_0
target_sample_size <- 100
sample_size_col <- (target_sample_size - 50) / 10 + 1 
data1 <- log10(results_out1[, sample_size_col, ])
data2 <- log10(results_out2[, sample_size_col, ])

par(mfrow=c(1,2))
plot_boxplot(data1, "Error1")
plot_boxplot(data2, "Error2")


# 3. Perform ANOVA
# Here we excluded RGCCA, you can also include it
# Data preparation
dat <- results_out2[, sample_size_col, ]
error_matrix <- rbind(dat[1:3,],dat[5,])
rownames(error_matrix) <- c("sumscore", "sem", "sam", "lslv")
error_data <- data.frame(
  Method = rep(rownames(error_matrix), each = max_j),
  Error = as.vector(t(error_matrix))
)
error_data$Method <- factor(error_data$Method, levels = c("sumscore", "sem", "sam", "lslv"))

# Perform anova
anova_result <- aov(Error ~ Method, data = error_data)
summary(anova_result)

# ad-hoc comparison
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# plot(tukey_result)

