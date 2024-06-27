SIMU <- function(X, Y, N, R){
  nobs <- N
  
  # Set up the loadings
  loading1 <- c(0,rep(sqrt(.6),10),rep(0,9))
  loading2 <- c(sqrt(.6),rep(0,10),rep(sqrt(.6),9))
  L <- cbind(loading1,loading2)
  
  # Split the data
  random_order <- sample(1:nobs)
  split_point <- round(nobs/2)
  X_train <- X[random_order[1:split_point], ]
  X_test <- X[random_order[(split_point + 1):nobs], ]
  Y_train <- Y[random_order[1:split_point]]
  Y_test <- Y[random_order[(split_point + 1):nobs]]
  
  #############################################################################
  # Method 1: Sum Score
  #############################################################################
  
  # 1. Transform the data type
  df_train = data.frame(cbind(X_train))
  colnames(df_train) = paste("x", c(1:20), sep = "")
  
  # 2. Factor analysis
  result <- fa(df_train, nfactors=2, rotate='varimax') 
  loadings <- as.matrix(result$loadings)
  
  # 3. Transform loadings matrix to 0/1/-1 format for the calculation of sum scores
  # For each row, find the index of the element with the largest absolute value
  max_abs_indices <- apply(loadings, 1, function(x) which.max(abs(x)))
  
  # Create a matrix of the same size filled with zeros
  transformed_loadings <- matrix(0, nrow = nrow(loadings), ncol = ncol(loadings))
  
  # Assign +1 or -1 to the element with the largest absolute value in each row
  for (i in 1:nrow(loadings)) {
    max_index <- max_abs_indices[i]
    max_value <- loadings[i, max_index]
    transformed_loadings[i, max_index] <- ifelse(max_value < 0, -1, 1)
  }
  
  res_loading <- transformed_loadings
  
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
  
  # 4. Calculate sum scores
  factor_score <- X_train %*% res_loading
  
  # 5. Regression
  res_reg <- lm(Y_train ~ factor_score[,1] + factor_score[,2])
  reg_coef <- unname(res_reg$coefficients)
  
  # 6. Obtain test scores and errors
  test_score <- X_test %*% res_loading
  Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score[,1]*reg_coef[2] + test_score[,2] * reg_coef[3]
  error_sumscore <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )
  
  
  #############################################################################
  # Method 2: SEM
  #############################################################################
  
  
  # 1 Transform the data type
  df_train = data.frame(cbind(X_train, Y_train))
  df_test = data.frame(cbind(X_test, Y_test))
  colnames(df_train) = c(paste("x", c(1:20), sep = ""), "y")
  colnames(df_test) = c(paste("x", c(1:20), sep = ""), "y")
  
  # 2. Estimate a lavaan model
  myModel <- '
y ~ factor1 + factor2
factor1 =~ NA*x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11
factor2 =~ NA*x1 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20
factor1 ~~ 1 * factor1
factor2 ~~ 1 * factor2
'
  
  fit <- sem(model=myModel, data=df_train)
  # constrain variance of factors to 1
  # orthoganol
  # sign of loadings
  # check intercept
  
  # 1.6 Predict the Y value
  Y_hat <- lavPredictY(fit, newdata = df_test, xnames = (paste("x", c(1:20), sep = "")), ynames = c("y"))
  error_sem <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )
  
  
  #############################################################################
  # Method 3: SAM
  #############################################################################
  
  
  # 1 Transform the data type
  df_train = data.frame(cbind(X_train, Y_train))
  df_test = data.frame(cbind(X_test, Y_test))
  colnames(df_train) = c('y1',paste("x", c(1:10), sep = ""), paste("y", c(2:10), sep = ""), "z")
  colnames(df_test) = c('y1',paste("x", c(1:10), sep = ""), paste("y", c(2:10), sep = ""), "z")
  
  data <- df_train
  S <- cov(data) * (nrow(df_train) - 1) / nrow(df_train)
  
  # 2. Estimate a lavaan model
  model <- '
    fz ~ fy + fx
    fy =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
    fx =~ NA*y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10
    fx ~~ 1 * fx
    fy ~~ 1 * fy
  '
  
  sub.est <- partitioned.estimation(data = as.data.frame(df_train), method = "SAM_MLB")
  
  summ.stats <- LSAM_SumStats(
    S = S,
    sample.nobs = nobs,
    LAMBDA = sub.est$LAMBDA,
    THETA = sub.est$THETA,
    mapping = "ML"
  )
  
  beta.coefs <- LSAM_regcoef(model = model, sumstats = summ.stats)
  
  res_loading <- sub.est$LAMBDA[1:20, 1:2]
  
  perm <- gtools::permutations(R, R)
  absdiff <- c()
  for (p in 1:nrow(perm)) {
    corsign <- sign(diag(cor(L, res_loading[, perm[p, ]])))
    L_res <- (res_loading[, perm[p, ]]) %*% diag(corsign)
    absdiff[p] <- sum(rowSums(abs(L - L_res)))
  }
  bestperm <- which.min(absdiff)
  loadings <- res_loading[, perm[bestperm, ]]
  corsign <- sign(diag(cor(L, loadings)))
  res_loading <- loadings %*% diag(corsign)
  beta.coefs <- beta.coefs %*% diag(corsign)  
  # doing a parallel transformation for beta, correct?
  # if the sign reverses for loadings, they should also be reversed for each column of beta
  
  test_score <- X_test %*% res_loading %*% solve(t(res_loading) %*% res_loading)
  Y_hat <- test_score[, 1] * beta.coefs[1] + test_score[, 2] * beta.coefs[2]
  error_sam <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )
  
  
  #############################################################################
  # Method 4: RGCCA
  #############################################################################
  
  # 1. Transform the data
  X_train_f1 <- X_train[,2:11]
  X_train_f2 <- cbind(X_train[,1],X_train[,12:20])
  X_test_f1 <- X_test[,2:11]
  X_test_f2 <- cbind(X_test[,1],X_test[,12:20])
  
  # 2. Estimate the model
  blocks <- list(bl1 = X_train, bl2 = Y_train)
  connection <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  sparsity <- c(1/sqrt(4),1)
  fit_sgcca <- rgcca(blocks = blocks, connection =  connection, method = "sgcca", sparsity = sparsity, superblock = FALSE, ncomp = c(2,1),
                     scheme = "factorial",
                     comp_orth = TRUE,
                     verbose = TRUE)

  # 3. Obtain factor scores
  w <- fit_sgcca$astar$bl1
  train_score <- X_train%*%w  # factor scores of training set
  res_loading <- t(solve(t(train_score)%*%train_score)%*%t(train_score)%*%X_train) # loadings
  
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
  
  # 4. Run the regression
  res_reg <- lm(Y_train ~ train_score[,1] + train_score[,2]) # do a regression and get coefficients
  reg_coef <- unname(res_reg$coefficients)
  
  # 5. Test
  test_score <- X_test %*% res_loading %*% solve(t(res_loading) %*% res_loading) # get the factor scores in test set
  Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score[,1]*reg_coef[2] + test_score[,2] * reg_coef[3] # get the prediction results
  
  error_rgcca <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )
  
  
  #############################################################################
  # Method 5: LSLV
  #############################################################################
  
  # 1. Run the function
  res_method <- LSLVLASSO(X_train,2,2,1e3,1e-4) # before running this, run LSLVLASSO.R first
  train_score <- res_method$scores # the factor score in training set
  res_loading <- res_method$loadings
  
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
  res_reg <- lm(Y_train ~ train_score[,1] + train_score[,2]) # do a regression and get coefficients
  reg_coef <- unname(res_reg$coefficients)
  
  # 3. Test
  test_score <- X_test %*% res_loading %*% solve(t(res_loading) %*% res_loading) # get the factor scores in test set
  Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score[,1]*reg_coef[2] + test_score[,2] * reg_coef[3] # get the prediction results
  error_lslv <- c( sum(Y_test-Y_hat)^2 / floor(nobs/2), sum(Y_test-Y_hat)^2 / sum(Y_test^2) )
  
  #############################################################################
  # Output
  #############################################################################
  out1 <- c(error_sumscore[1], error_sem[1], error_sam[1], error_rgcca[1], error_lslv[1])
  out2 <- c(error_sumscore[2], error_sem[2], error_sam[2], error_rgcca[2], error_lslv[2])
  return(list(out1 = out1, out2 = out2))
  
}