#### Title: Regularized LS-LV method
#### Author: Tra Le 
#### Supervisor: Dr. Katrijn Van Deun
#### Created: February 1, 2023
#### Last modified: May 29, 2023

#########################################################################################
#####                  RLSLV with LASSO for correlated components                   #####
#########################################################################################

#1. LOSS FUNCTION WITH THE LASSO PENALTY
LOSS <- function(DATA, SCORES, LOADINGS, LAMBDA){
  XHAT <- SCORES%*%t(LOADINGS)
  res <- sum(rowSums((XHAT-DATA)^2))
  penalty <- sum(abs(LOADINGS))
  loss <- res+LAMBDA*penalty
  return(loss)
}

#2. REGULARIZED LSLV WITH THE LASSO
# R: number of factors
# eps: convergence criteria usually 1e-4
# lambda: tuning parameter, non-negative, >=0
LSLVLASSO <- function(DATA, R, lambda, MaxIter, eps){
  I <- dim(DATA)[1]
  J <- dim(DATA)[2]
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc   #Used to compute convergence criterium
  svd1 <- svd(DATA, R, R)
  P2 <- svd1$v %*% diag(svd1$d[1:R])/sqrt(I)
  P1 <- matrix(rnorm(J*R), nrow = J, ncol = R) #initialize P 
  loadings <- .3*P1 + .7*P2
  scores <- svd1$u 
  diffT <- 0
  diffP <- 0
  while (convAO == 0) {
    iter0 <- 1
    Losst <- 1
    Lossvec0 <- Losst
    convT0 <- 0
    Lossvec1 <- 1
    #1. Update component scores 
    while(convT0 == 0){
      Lossu1old <- LOSS(DATA,scores,loadings,lambda)
      E <- DATA - scores%*%t(loadings) 
      for (r in 1:R){
        Er <- E + scores[,r]%*%t(loadings[,r])
        num <- Er%*%loadings[,r]
        scores[,r] <- sqrt(I)*num/sqrt(sum(num^2))
        Lossu1 <- LOSS(DATA,scores,loadings,lambda)
        diffT <- c(diffT,Lossu1old-Lossu1)
        Lossu1old <- Lossu1
        Lossvec1 <- c(Lossvec1, Lossu1)
      }
      #t(scores)%*%scores
      #Calculate loss
      Lossu0 <- LOSS(DATA,scores,loadings,lambda)
      Lossvec0 <- c(Lossvec0,Lossu0)
      # check convergence
      if (iter0 > MaxIter) {
        convT0 <- 1
      }
      if (abs(Losst-Lossu0) < eps){
        convT0 <- 1
      }
      iter0 <- iter0 + 1
      Losst <- Lossu0
    }
    Loss <- LOSS(DATA, scores, loadings, lambda)
    
    #2. Update loadings
    Lossu1old <- LOSS(DATA,scores,loadings,lambda)
    E <- DATA - scores%*%t(loadings) 
    for (r in 1:R){
      Er <- E+scores[,r]%*%t(loadings[,r])
      crosstEr <- t(Er)%*%scores[,r]
      loadings[,r]<-sign(crosstEr)*apply(cbind(abs(crosstEr)-lambda/2,0),1,max)/I
    }
    
    #Calculate loss
    Lossu <- LOSS(DATA,scores,loadings,lambda)
    Lossvec <- c(Lossvec,Lossu)
    if (iter > MaxIter) {
      convAO <- 1
    }
    if (abs(Lossc-Lossu) < eps) {
      convAO <- 1
    }
    iter <- iter + 1
    Lossc <- Lossu
  }
  return_rlslv <- list()
  return_rlslv$loadings <- loadings
  return_rlslv$scores <- scores
  return_rlslv$Loss <- Loss
  return_rlslv$Lossvec <- Lossvec
  
  return(return_rlslv)
}

#3. MULTISTART PROCEDURE

MULTISTART_LSLVLASSO <- function(DATA, R, MaxIter, eps, nstarts, lambda){
  if(missing(nstarts)){
    nstarts <- 20
  } 
  
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  
  for (n in 1:nstarts){
    result <- LSLVLASSO(DATA, R, lambda, MaxIter, eps)
    
    Pout3d[[n]] <- result$loadings
    Tout3d[[n]] <- result$scores
    LOSS[n] <- result$Loss
    LOSSvec[[n]] <- result$Lossvec
  }
  # choose solution with lowest loss value
  k <- which(LOSS == min(LOSS))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  
  return_varselect <- list()
  return_varselect$loadings <- Pout3d[[k]]
  return_varselect$scores <- Tout3d[[k]]
  return_varselect$Lossvec <- LOSSvec
  return_varselect$Loss <- LOSS[k]
  
  return(return_varselect)
}

###function for recovery rate
num_correct <- function (TargetP, EstimatedP){
  total_vnumber <- dim(TargetP)[1] * dim(TargetP)[2]
  TargetP[which(TargetP != 0)] <- 1
  sum_select <- sum(TargetP)
  sum_zero <- total_vnumber - sum_select
  EstimatedP[which(EstimatedP != 0)] <- 1
  total_correct <- sum(TargetP == EstimatedP) # this is the total number of variables correctedly selected and zeros correctly retained
  prop_correct <- total_correct/total_vnumber
  return(prop_correct)
}

#################################################################################################
# Function to tune lasso, so that the chosen lasso leads to
# The number of zeroes in the generated data
#################################################################################################

findLasso <- function(dat, Ptrue, maxItr, lassou){
  
  percentageZeroesInData <- sum(Ptrue == 0)
  percentageInW <- 0
  i  <- 0
  lassol <- 0 
  converged <- FALSE
  conv0 <- 0
  lasso1 <- 1
  while(conv0 == 0){
    
    lasso <- (lassol + lassou) / 2
    fit <- MULTISTART_LSLVLASSO(dat, lambda = lasso, R = 3, MaxIter = 200, eps = 10^-6, nstarts = 20)
    percentageInW <- sum(round(fit$loadings,3) == 0)
    if( percentageZeroesInData > percentageInW){
      lassol  <- lasso
      
    } else {
      lassou  <- lasso
    }
    #print(lasso)
    
    if (abs(percentageZeroesInData - percentageInW) < 0){
      conv0 <- 1
    }
    
    
    if(i > maxItr){
      conv0 <- 1
    }
    
    
    if (abs(lasso1 - lasso) < .01){
      conv0 <- 1
    }
    lasso1 <- lasso
    i <- i + 1
  }
  
  if( i < maxItr ){
    converged <- TRUE
  } 
  
  return(list(lasso = lasso, converged = converged))
}

##### MODEL SELECTION
# index of sparseness 
IS_LSLVLASSO <- function(DATA, R, lambda, MaxIter, eps, nstarts){
  J <- dim(DATA)[2]
  
  VarSelect0 <- MULTISTART_LSLVLASSO(DATA, R, lambda = 0, MaxIter, eps, nstarts)
  P_hat0 <- VarSelect0$loadings
  T_hat0 <- VarSelect0$scores
  
  V_oo <- sum(DATA^2)
  V_s <- sum((T_hat0%*%t(P_hat0))^2) 
  
  VarSelect <- MULTISTART_LSLVLASSO(DATA, R, lambda, MaxIter, eps, nstarts)
  P_hat <- VarSelect$loadings
  T_hat <- VarSelect$scores
  
  V_a <- sum((T_hat %*% t(P_hat))^2)
  IS <- list()
  IS$value <- (V_a * V_s / V_oo^2) * (sum(round(P_hat,3) == 0) /(J*R))
  IS$vaf <- V_a/V_oo
  IS$nzero <- sum(round(P_hat,3) == 0)
  return(IS)
}



