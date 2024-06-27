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
source("LSLVLASSO.R")


#############################################################################
# Data Generation
#############################################################################

dat <- scale(PoliticalDemocracy)
N <- nobs <- nrow(dat)  # Sample size
R <- 3  # Nr. of factors

random_order <- sample(1:nobs)
split_point <- round(nobs/2)
data_train <- dat[random_order[1:split_point], ]
data_test <- dat[random_order[(split_point + 1):nobs], ]


#############################################################################
# Method 1: Sum Score
#############################################################################

# 1. Transfrom data
data_train_sumscore <- as.data.frame(data_train)
data_test_sumscore <- as.data.frame(data_test)

# 2. Perform linear regression
fit <- lm(cbind(y5,y6,y7,y8) ~ cbind(y1,y2,y3,y4,x1,x2,x3), data = data_train_sumscore)
Y_hat <- predict(fit, newdata = data_test_sumscore)
Y_test <- data_test[,5:8]
error_sumscore <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )

#############################################################################
# Method 2: SEM
#############################################################################

# 1. Estimate a lavaan model
myModel <- '
  # measurement model
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
    ind60 =~ x1 + x2 + x3
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit <- sem(model=myModel, data=data_train)

# 2. Predict the Y value
Y_hat <- lavPredictY(fit, newdata = data_test, ynames = c('y5','y6','y7','y8'),xnames = c('x1','x2','x3','y1','y2','y3','y4'))
Y_test <- data_test[,5:8]
error_sem <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )


#############################################################################
# Method 3: SAM
#############################################################################

# 1. Rename the column names to fit the need of SAM code
data_sam <- data_train
colnames(data_sam)[1:4] <- c("x1", "x2", "x3", "x4")
colnames(data_sam)[5:8] <- c("y1", "y2", "y3", "y4")
colnames(data_sam)[9:11] <- c("z1", "z2", "z3")
S <- cov(data_sam) * (nrow(data_sam) - 1) / nrow(data_sam)

# 2. Estimate a SAM model
model <- '
  # measurement model
    fx =~ x1 + x2 + x3 + x4
    fy =~ y1 + y2 + y3 + y4
    fz =~ z1 + z2 + z3
  # regressions
    fx ~ fz
    fy ~ fx + fz
  # residual correlations
    x1 ~~ y1
    x2 ~~ x4 + y2
    x3 ~~ y3
    x4 ~~ y4
    y2 ~~ y4
'

sub.est <- partitioned.estimation(data = data_sam, method = "SAM_MLB")

# 3. Obtain Statistics
summ.stats <- LSAM_SumStats(
  S = S,
  sample.nobs = nrow(data_sam),
  LAMBDA = sub.est$LAMBDA,
  THETA = sub.est$THETA,
  mapping = "ML"
)

beta.coefs <- LSAM_regcoef(model = model, sumstats = summ.stats)

# 3. Calculate loadings and factor scores (TMatrix and PMatrix)
loadings <- sub.est$LAMBDA
loadings_x <- sub.est$LAMBDA[c(1:4,9:11), c(1,3)]
#loadings_y <- sub.est$LAMBDA[5:8, 2]

x_score <- data_test[,c(1:4, 9:11)] %*% loadings_x %*% solve(t(loadings_x) %*% loadings_x)
y_score <- x_score[,1] * beta.coefs[2] + x_score[,2] * beta.coefs[3]
total_score <- cbind(x_score[,1], y_score, x_score[,2])

# 4. Estimate Y
data_hat <- total_score %*% t(loadings)
Y_hat <- data_hat[,5:8]
Y_test <- data_test[,5:8]

error_sam <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )

#############################################################################
# Method 4: RGCCA
#############################################################################

# 1. Transform the data
data_train_dem60 <-data_train[,1:4]
data_train_dem65 <- data_train[,5:8]
data_train_ind60 <- data_train[,9:11]

# 2. Estimate the model
blocks <- list(bl1 = data_train_dem60, bl2 = data_train_dem65, bl3 = data_train_ind60)
sparsity <- as.data.frame(matrix(c(1, 1, 1), ncol = 3, byrow = T)) # Again, how should I set it?
fit_sgcca <- rgcca(blocks = blocks, connection =  1-diag(3), method = "sgcca", sparsity = sparsity, superblock = FALSE, ncomp = c(1,1,1),
                   scheme = "factorial",
                   comp_orth = TRUE,
                   verbose = TRUE)

# 3. Obtain loadings and factor scores
w <- rbind(cbind(fit_sgcca$astar$bl1,rep(0,4),rep(0,4)), cbind(rep(0,4),fit_sgcca$astar$bl2,rep(0,4)), cbind(rep(0,3),rep(0,3), fit_sgcca$astar$bl3))  # assigning non-loadings to 0
train_score <- as.matrix(data_train)%*%w  # factor scores of training set
res_loading <- t(solve(t(train_score)%*%train_score)%*%t(train_score)%*%as.matrix(data_train)) # loadings

# 4. Run the regression
res_reg <- lm(train_score[,2] ~ train_score[,1] + train_score[,3]) # do a regression and get coefficients
reg_coef <- unname(res_reg$coefficients)

# 5. Estimate factor scores of Y and the scores of the indicators
w_x <- rbind(cbind(fit_sgcca$astar$bl1,rep(0,4)), cbind(rep(0,3),fit_sgcca$astar$bl3))
test_score_x <- as.matrix(cbind(data_test[,1:4],data_test[,9:11])) %*% w_x
test_score_y <- rep(reg_coef[1],nrow(test_score_x)) + test_score_x[,1]*reg_coef[2] + test_score_x[,2] * reg_coef[3] # get the prediction results
total_score <- cbind(test_score_x[,1], test_score_y, test_score_x[,2])

data_hat <- total_score %*% t(res_loading)
Y_hat <- data_hat[,5:8]
Y_test <- data_test[,5:8]

error_rgcca <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )


#############################################################################
# Method 5: LSLV
#############################################################################

# 1. Run the function
data_lslv <- as.matrix(data_train)
res_method <- LSLVLASSO(as.matrix(data_train),3,2,1e3,1e-4)
train_score <- res_method$scores # the factor score in training set
res_loading <- res_method$loadings

L <- matrix(0, nrow = 11, ncol = 3) # Create the ideal loading matrix
L[1:4, 1] <- 1
L[5:8, 2] <- 1
L[9:11, 3] <- 1

perm <- gtools::permutations(R, R)
absdiff <- c()
for (p in 1:nrow(perm)) {
  corsign <- sign(diag(cor(L, res_loading[, perm[p,]])))
  L_res <- (res_loading[, perm[p,]]) %*% diag(corsign)
  absdiff[p] <- sum(rowSums(abs(L - L_res)))
}
bestperm <- which.min(absdiff)
loadings <- res_loading[, perm[bestperm,]]
corsign <- sign(diag(cor(L, loadings)))
res_loading <- loadings %*% diag(corsign)

# do a parallel transformation in train_score
scores <- train_score[, perm[bestperm,]]
train_score <- scores %*% diag(corsign)

# 2. Run the regression
res_reg <- lm(train_score[,2] ~ train_score[,1] + train_score[,3]) # do a regression and get coefficients
reg_coef <- unname(res_reg$coefficients)


loadings_x <- res_loading[c(1:4, 9:11),cbind(1,3)]

# 3. Test
score_x <- as.matrix(data_test[,c(1:4, 9:11)]) %*% loadings_x %*% solve(t(loadings_x) %*% loadings_x)
score_y <- rep(reg_coef[1],nrow(score_x)) + score_x[,1] * reg_coef[2] + score_x[,2] * reg_coef[3]
total_score <- cbind(score_x[,1], score_y, score_x[,2])

data_hat <- total_score %*% t(res_loading)
Y_hat <- data_hat[,5:8]
Y_test <- data_test[,5:8]

error_lslv <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )
