
                               # Portfolio Optimization Machine #

# ------------------------------------------------------------------------------------------------------ #

# One-Over-N Optimization ------------------------------------------------------------------------------ #

OneOverN.Optimizer <- function(x, ...)
{
  # It derives the 1/N portfolio weights
  return(rep(1 / ncol(x), ncol(x)))
}

# Max Sharpe-Ratio ------------------------------------------------------------------------------------- #

Max.Sharpe <- function(R, Views, mean.vect, cap.weight, TargetRisk, BL= FALSE, mN = 2000, Rfree=0)
{
  # It derives the vector of portfolio weights of the market or tangency portfolio
  # accoding to the Markowitz MPT. Long-only by hardcode. 
  #
  # R           - Returns data
  # Views       - Adjust equlibrium returns with views ? yes/ no (TRUE/FALSE)
  # mean.vect   - expected return: either the future mean (perfect view), random expected return (random view) or historical mean
  # cap.weight  - weight of the market portfolio
  # Target.Risk - Variable from Backtest function not relevant in this function
  # BL          - Estimate expected returns with Black-Litterman ? yes/ no (TRUE/FALSE)
  # mN          - Number of returns for which the efficient frontier is calculated 
  # Rfree       - Risk free rate
  
  # calculate expected return vector
  if(BL==FALSE){
    mean_vect <- apply(R, 2, mean) 
  }else{
    mean_vect <- BlackLitterman(R, Views, mean.vect, cap.weight)
  }
  # calculate covariance + target returns
  cov_mat   <- cov(R)
  N         <- ncol(R)
  muP       <- seq(min(mean_vect) + 0.001, max(mean_vect) - 0.001, length=mN) 
  
  # constraint matrix
  Amat      <- cbind(rep(1, N), mean_vect, diag(1, N), diag(-1, N)) 
  
  QuadProg  <- function(i)
  {
    bvec    <-c(1, muP[i], rep(wmin,N), rep(-wmax,N)) # constraint vector
    result  <- tryCatch(solve.QP(Dmat = 2*cov_mat, dvec = rep(0, N), Amat = Amat, bvec = bvec, meq=2),
                        error=function(e){return(0)}) # tryCatch : avoids error interruption if no solution is found
    if(length(result) == 1){
      sdP     <- NA
      weights <- rep(NA,N)
    }else{
      sdP     <- sqrt(result$value)
      weights <- result$solution 
    }
    return(list(sdP, weights))
  }
  
  result <- lapply(1:mN, QuadProg)
  
  sdP      <- sapply(result, `[[`, 1) # Get standard deviations
  weights  <- sapply(result, `[[`, 2) # Get weights
  sharpe   <- (muP-Rfree)/sdP # Compute Sharpe’s ratios
  
  weights <- weights[, which.max(sharpe)]
  weights[weights<0] <- 0 # Eliminate even super small negative weights
  return(weights)
}

## Target Risk Optimization ################################################################################

# Target VaR --------------------------------------------------------------------------------------------- #

Target.VaR <- function(R, Views, mean.vect, cap.weight, TargetRisk, BL= FALSE)
{
  # The function maximizes portfolio return subject to a VaR constraint ("Target Risk"). Long-only by hardcode. 
  #
  # R           - Returns data
  # Views       - Adjust equlibrium returns with views ? yes/ no (TRUE/FALSE)
  # mean.vect   - expected return: either the future mean (perfect view), random expected return (random view) or historical mean
  # cap.weight  - weight of the market portfolio
  # Target.Risk - Maximize portfolio return constraint to a target risk ? yes/ no (TRUE/FALSE) 
  # BL          - Estimate expected returns with Black-Litterman ? yes/ no (TRUE/FALSE)
  
  # calculate expected return vector
  if(BL==FALSE){
    mean_vect <- apply(R, 2, mean) 
  }else{
    mean_vect <- BlackLitterman(R, Views, mean.vect, cap.weight)
  }
  
  # number of assets
  N         <- ncol(R) 
  
  # startvalues for x
  NP          <- 100 * ncol(R)
  startvalues <- sobol(NP, dim = ncol(R), init=TRUE, scrambling=1) / ncol(R) # "/ncol(R)" to make maximal sum = 1
  tolerance   <- 1e-06 # weight constraint tolerance
  
  optVaR <- function(x)
  {
    tmpmat <- rbind(as.matrix(x - wmax, ncol=1),         # penalty for small weights
                    as.matrix(wmin - x, ncol=1),         # penalty for large weights
                    sum(x)-1,                            # penalty for leverage portfolios
                    1-sum(x))                            # penalty for low invested portfolios
    
    weight.penalty <- sum(pmax(tmpmat - tolerance,0)) * 100 
    risk.penalty   <- max(TargetRisk - VaR_FUN(R %*% x),0)*100 # penalty for risk exceeding target risk
    obj            <- -(x %*% mean_vect)
    return( obj + weight.penalty + risk.penalty) 
  }
  
  controlDE <- list(NP = NP, initialpop = startvalues, itermax=1000, 
                    reltol = 0.000001, steptol = 150, c = 0.4, trace=FALSE)
  
  # run the parallel version
  registerDoMC()
  set.seed(1234)
  
  result <- DEoptim(optVaR,lower=rep(0,N),upper=rep(1,N), control=controlDE)
  
  # stop parallel processing
  wait()
  
  # For windows user:
  # #run the parallel version
  # cl <- makeSOCKcluster(3)
  # clusterEvalQ(cl, library(randtoolbox)) # load any necessary libraries
  # clusterExport(cl, 
  #               list("R","N","NP","optVaR", "tolerance", "startvalues", "VaR_FUN","controlDE", 
  #                    "TargetRisk","wmax")) # copy any necessary objects
  # registerDoSNOW(cl) # register foreach backend
  # set.seed(1234)
  # 
  # 
  # result <- DEoptim(optMDD,lower=rep(0,N),upper =rep(1,N), control=controlDE)
  # 
  # # stop parallel processing
  # stopCluster(cl) # stop cluster
  
  w      <- result$optim$bestmem
  
  if(sum(w)>1) w<- w/sum(w)
  return(w)
}

# Target CVaR -------------------------------------------------------------------------------------------- #

Target.CVaR <- function(R, Views, mean.vect, cap.weight, TargetRisk, BL= FALSE) 
{
  # The function maximizes portfolio return subject to a CVaR constraint ("Target Risk"). Long-only by hardcode. 
  #
  # R           - Returns data
  # Views       - Adjust equlibrium returns with views ? yes/ no (TRUE/FALSE)
  # mean.vect   - expected return: either the future mean (perfect view), random expected return (random view) or historical mean
  # cap.weight  - weight of the market portfolio
  # Target.Risk - Maximize portfolio return constraint to a target risk ? yes/ no (TRUE/FALSE) 
  # BL          - Estimate expected returns with Black-Litterman ? yes/ no (TRUE/FALSE)
  
  # create mean_vect
  if(BL==FALSE){
    mean_vect <- apply(R, 2, mean) 
  }else{
    mean_vect <- BlackLitterman(R, Views, mean.vect, cap.weight)
  }
  
  N         <- ncol(R) # number of assets
  S         <- nrow(R) # number of scenarios i.e. periods
  
  # constraint matrix
  # weights # alpha # auxiliary variables
  Amat <- rbind(c(rep(1,N), 0, rep(0,S)),
                c(rep(0,N), 1, rep(1/(alpha*S), S)),
                cbind(R, 1, diag(S)))
  
  # create objective vector
  objL <- c(mean_vect, rep(0, S+1)) #c^T
  
  # direction vector (supports general equality/inequality constraints)
  dir.vec <- c("<=","<=",rep(">=",S))
  
  # bounds on weights
  bounds <- list(lower = list(ind = 1:N, val = rep(wmin,N)),
                 upper = list(ind = 1:N, val = rep(wmax,N)))
  
  bvec    <- c(1,-(TargetRisk),rep(0,S))
  result  <- Rglpk_solve_LP(obj=objL, mat=Amat, dir=dir.vec, rhs=bvec, 
                            types=rep("C",length(objL)), max=T, bounds=bounds)
  
  w       <- as.numeric(result$solution[1:N])
  return(w)
}

# Target MDD -------------------------------------------------------------------------------------------- #

Target.MDD <- function(R, Views, mean.vect, cap.weight, TargetRisk, BL= FALSE) 
{
  # function orientated on the function PMaxDD() of the FRAPO r-package
  # The function maximizes portfolio return subject to a MDD constraint ("Target Risk"). Long-only by hardcode. 
  #
  # R           - Returns data
  # Views       - Adjust equlibrium returns with views ? yes/ no (TRUE/FALSE)
  # mean.vect   - expected return: either the future mean (perfect view), random expected return (random view) or historical mean
  # cap.weight  - weight of the market portfolio
  # Target.Risk - Maximize portfolio return constraint to a target risk ? yes/ no (TRUE/FALSE) 
  # BL          - Estimate expected returns with Black-Litterman ? yes/ no (TRUE/FALSE)
  
  RC        <- apply(R,2,cumsum)
  
  N <- ncol(RC)
  J <- nrow(RC)
  w <- rep(0, N) ## weights
  u <- rep(0, J) ## high-watermark
  x <- c(w, u)
  
  # Defining objective 
  if(BL==FALSE){
    # objective = expected end-wealth
    obj <- c(as.numeric(RC[J, ]), rep(0, J))
  }else{
    # objective = max BL - expected return
    r_hat <- BlackLitterman(R, Views, mean.vect, cap.weight)*252 
    obj <- c(r_hat, rep(0, J))
  }
  
  ## a1: budget constraint
  a1 <- c(rep(1, N), rep(0, J))
  d1 <- "<="
  b1 <- 1
  ## a2: draw-down constraint (1)
  a2 <- cbind(-1 * RC, diag(J))
  d2 <- rep("<=", J)
  b2 <- rep(-TargetRisk, J)
  ## a3: draw-down constraint (2)
  a3 <- a2
  d3 <- rep(">=", J)
  b3 <- rep(0, J)
  ## a4: draw-down constraint (3)
  D1 <- -1.0 * diag(J)
  udiag <- embed(1:J, 2)[, c(2, 1)] 
  D1[udiag] <- 1
  a4 <- cbind(matrix(0, ncol = N, nrow = J), D1)
  a4 <- a4[-J, ]
  d4 <- rep(">=", J-1)
  b4 <- rep(0, J-1)
  ## a5: draw-down constraint (4)
  a5 <- c(rep(0, N), 1, rep(0, J - 1))
  d5 <- "=="
  b5 <- 0  
  ## Combining restrictions
  Amat <- rbind(a1, a2, a3, a4, a5)
  Dvec <- c(d1, d2, d3, d4, d5)
  Bvec <- c(b1, b2, b3, b4, b5)

  # bounds on weights
  bounds <- list(lower = list(ind = 1:N, val = rep(wmin,N)),
                 upper = list(ind = 1:N, val = rep(wmax,N)))
  
  ## Solving LP
  result <- Rglpk_solve_LP(obj = obj, mat = Amat, dir = Dvec, rhs = Bvec,
                           max = TRUE, bounds=bounds)
  if(result$status != 0){
    warning(paste("GLPK had exit status:", opt$status))
  }
  
  w       <- as.numeric(result$solution[1:N])
  return(w)
  

}

# Target AvDD -------------------------------------------------------------------------------------------- #
Target.AvDD <- function(R, Views, mean.vect, cap.weight, TargetRisk, BL= FALSE) 
{
  # function orientated on the function PAveDD() of the FRAPO r-package
  # The function maximizes portfolio return subject to a MDD constraint ("Target Risk"). Long-only by hardcode. 
  #
  # R           - Returns data
  # Views       - Adjust equlibrium returns with views ? yes/ no (TRUE/FALSE)
  # mean.vect   - expected return: either the future mean (perfect view), random expected return (random view) or historical mean
  # cap.weight  - weight of the market portfolio
  # Target.Risk - Maximize portfolio return constraint to a target risk ? yes/ no (TRUE/FALSE) 
  # BL          - Estimate expected returns with Black-Litterman ? yes/ no (TRUE/FALSE)
  
  RC        <- apply(R,2,cumsum)
  
  N <- ncol(RC)
  J <- nrow(RC)
  w <- rep(0, N) ## weights
  u <- rep(0, J) ## high-watermark
  v <- rep(0, J) ## draw downs
  x <- c(w, u, v)
  
  ## Defining objective (end-wealth)
  if(BL==FALSE){
    obj <- c(as.numeric(RC[J, ]), rep(0, J), rep(0, J))
  }else{
    r_hat <- BlackLitterman(R, Views, mean.vect, cap.weight)*252
    obj <- c(r_hat, rep(0, J), rep(0, J))
  }
  
  ## a1: budget constraint
  a1 <- c(rep(1, N), rep(0, 2 * J))
  d1 <- "<="
  b1 <- 1
  ## a2: draw-down constraint (1) assigning summands to v
  a2 <- cbind(-1 * RC, diag(J), -1 * diag(J))
  d2 <- rep("==", J)
  b2 <- rep(0, J)
  ## a3: defining average constraint
  a3 <- c(rep(0, N), rep(0, J), rep(1 / J, J))
  d3 <- "<="
  b3 <- -TargetRisk
  ## a4: draw-down constraint (2)
  a4 <- cbind(-1 * RC, diag(J), matrix(0, nrow = J, ncol = J))
  d4 <- rep(">=", J)
  b4 <- rep(0, J)
  ## a5: draw-down constraint (3)
  D1 <- -1.0 * diag(J)
  udiag <- embed(1:J, 2)[, c(2, 1)] 
  D1[udiag] <- 1
  a5 <- cbind(matrix(0, ncol = N, nrow = J), D1, matrix(0, ncol = J, nrow = J))
  a5 <- a5[-J, ]
  d5 <- rep(">=", J-1)
  b5 <- rep(0, J-1)
  ## Combining restrictions
  Amat <- rbind(a1, a2, a3, a4, a5)
  Dvec <- c(d1, d2, d3, d4, d5)
  Bvec <- c(b1, b2, b3, b4, b5)

  # bounds on weights
  bounds <- list(lower = list(ind = 1:N, val = rep(wmin,N)),
                 upper = list(ind = 1:N, val = rep(wmax,N)))
  
  ## Solving LP
  result <- Rglpk_solve_LP(obj = obj, mat = Amat, dir = Dvec, rhs = Bvec,
                           max = TRUE, bounds=bounds)
  if(result$status != 0){
    warning(paste("GLPK had exit status:", opt$status))
  }
  
  w       <- as.numeric(result$solution[1:N])
  return(w)
}

## Target Risk Optimization ###############################################################################

# Minimum Variance-----------------------------------------------------------------------------------------

GlobalMinVar <- function(R, mN = 2000, ...)
{
  # It derives the vector of portfolio weights of the Global Minimum Variance portfolio
  # accoding to the Markowitz MPT. Long-only by hardcode. Expected returns as historical 
  # means by hardcode.
  #
  # R           - Return data
  # mN          - Number of returns for which the efficient frontier is calculated 
  
  # create mean_vect
  if(BL==FALSE){
    mean_vect <- apply(R, 2, mean) 
  }else{
    mean_vect <- BlackLitterman(R, mean.vect)
  }
  
  cov_mat   <- cov(R)
  N         <- ncol(R)
  muP       <- seq(min(mean_vect) + 0.001, max(mean_vect) - 0.001, length=mN) # Target returns
  
  Amat      <- cbind(rep(1, N), mean_vect, diag(1, N), diag(-1, N))
  
  # QuadProg Function
  QuadProg  <- function(i)
  {
    bvec    <-c(1, muP[i], rep(wmin,N), rep(-wmax,N)) # constraint vector
    result  <- tryCatch(solve.QP(Dmat = 2*cov_mat, dvec = rep(0, N), Amat = Amat, bvec = bvec, meq=2),
                        error=function(e){return(0)})
    if(length(result) == 1){
      sdP     <- NA
      weights <- rep(NA,N)
    }else{
      sdP     <- sqrt(result$value)
      weights <- result$solution 
    }
    return(list(sdP, weights))
  }
  
  result <- lapply(1:mN, QuadProg)
  
  sdP      <- sapply(result, `[[`, 1) # Get standard deviations
  weights  <- sapply(result, `[[`, 2) # Get weights
  
  weights <- weights[, which.min(sdP)]
  weights[weights<0] <- 0 # eliminate even super small negative weights
  return(weights)
}

# Minimum CVaR ----------------------------------------------------------------------------------------------

Min.CVaR <- function(R,...) 
{
  # The function returns the weights of the minimum-CVaR portfolio. Long-only by hardcode. 
  #
  # R           - Returns data
  
  N         <- ncol(R)
  S         <- nrow(R) # number of scenarios i.e. periods
  
  # constraint matrix
  Amat = rbind(c(rep(1,N), 0, rep(0,S)),
               cbind(R, 1, diag(S)))
  
  # create objective vector
  objL <- c(rep(0,N), -1, rep(-1/(alpha*S), S)) #c^T
  
  # direction vector (supports general equality/inequality constraints)
  dir.vec <- c("==",rep(">=",S))
  
  # bounds on weights
  bounds <- list(lower = list(ind = 1:N, val = rep(wmin,N)),
                 upper = list(ind = 1:N, val = rep(wmax,N)))
  
  bvec    <- c(1,rep(0,S))
  result  <- Rglpk_solve_LP(obj=objL, mat=Amat, dir=dir.vec, rhs=bvec, 
                            types=rep("C",length(objL)), max=T, bounds=bounds)
  
  w       <- as.numeric(result$solution[1:N])
  return(w)
}


# Minimum VaR -------------------------------------------------------------------------------------------- #

Min.VaR <- function(R, ...)
{
  # The function returns the weights of the minimum-VaR portfolio. Long-only by hardcode. 
  #
  # R           - Returns data
  
  # number of assets
  N         <- ncol(R) 
  
  # startvalues for x
  NP          <- 100 * ncol(R)
  startvalues <- sobol(NP, dim = ncol(R), init=TRUE, scrambling=1) / ncol(R) # "/ncol(R)" to make maximal sum = 1
  tolerance   <- 1e-06 # weight constraint tolerance
  
  optVaR <- function(x)
  {
    tmpmat <- rbind(as.matrix(x - wmax, ncol=1),         # penalty for small weights
                    as.matrix(wmin - x, ncol=1),            # penalty for large weights
                    sum(x)-1,                            # penalty for leverage portfolios
                    1-sum(x))                            # penalty for low invested portfolios
    
    penalty <- sum(pmax(tmpmat - tolerance,0)) * 100
    obj     <- abs(VaR_FUN(R %*% x))
    return( obj + penalty) 
  }
  
  controlDE <- list(NP = NP, initialpop = startvalues, itermax=1000, 
                    reltol = 0.000001, steptol = 150, c = 0.4, trace=FALSE)
  
  #run the parallel version
  registerDoMC()
  set.seed(1234)
  
  result <- DEoptim(optVaR,lower=rep(0,N),upper =rep(1,N), control=controlDE)
  
  # stop parallel processing
  wait()
  
  # For windows user:
  # #run the parallel version
  # cl <- makeSOCKcluster(3)
  # clusterEvalQ(cl, library(randtoolbox)) # load any necessary libraries
  # clusterExport(cl, 
  #               list("R","N","NP","optVaR", "tolerance", "startvalues", "VaR_FUN","controlDE", 
  #                    "wmax")) # copy any necessary objects
  # registerDoSNOW(cl) # register foreach backend
  # set.seed(1234)
  # 
  # 
  # result <- DEoptim(optMDD,lower=rep(0,N),upper =rep(1,N), control=controlDE)
  # 
  # # stop parallel processing
  # stopCluster(cl) # stop cluster
  
  w      <- result$optim$bestmem
  
  if(sum(w)>1) w<- w/sum(w) # ensure long-only portfolio
  return(w)
}


# Minimum MDD -------------------------------------------------------------------------------------------- #

Min.MDD <- function(R, ...)
{
  # The function returns the weights of the minimum-MDD portfolio. Long-only by hardcode. 
  #
  # R           - Returns data
  
  # number of assets
  N         <- ncol(R) 
  
  # startvalues for x
  NP          <- 100 * ncol(R)
  startvalues <- sobol(NP, dim = ncol(R), init=TRUE, scrambling=1) / ncol(R) # "/ncol(R)" to make maximal sum = 1
  tolerance   <- 1e-06 # weight constraint tolerance
  
  optMDD <- function(x)
  {
    tmpmat <- rbind(as.matrix(x - wmax, ncol=1),         # penalty for small weights
                    as.matrix(wmin - x, ncol=1),         # penalty for large weights
                    sum(x)-1,                            # penalty for leverage portfolios
                    1-sum(x))                            # penalty for low invested portfolios
    
    penalty <- sum(pmax(tmpmat - tolerance,0)) * 100
    obj     <- abs(MDD_FUN(R %*% x))
    return( obj + penalty) 
  }
  
  controlDE <- list(NP = NP, initialpop = startvalues, itermax=1000, 
                    reltol = 0.000001, steptol = 150, c = 0.4, trace=FALSE)
  
  #run the parallel version
  registerDoMC()
  set.seed(1234)
  
  result <- DEoptim(optVaR,lower=rep(0,N),upper =rep(1,N), control=controlDE)
  
  # stop parallel processing
  wait()
  
  # For windows user:
  # #run the parallel version
  # cl <- makeSOCKcluster(3)
  # clusterEvalQ(cl, library(randtoolbox)) # load any necessary libraries
  # clusterExport(cl, 
  #               list("R","N","NP","optMDD", "tolerance", "startvalues", "MDD_FUN","controlDE", 
  #                    "wmax")) # copy any necessary objects
  # registerDoSNOW(cl) # register foreach backend
  # set.seed(1234)
  # 
  # 
  # result <- DEoptim(optMDD,lower=rep(0,N),upper =rep(1,N), control=controlDE)
  # 
  # # stop parallel processing
  # stopCluster(cl) # stop cluster
  
  w      <- result$optim$bestmem
  
  if(sum(w)>1) w<- w/sum(w) # ensure long-only portfolio
  return(w)
}

# Minimum AvDD ------------------------------------------------------------------------------------------- #

Min.AvDD <- function(R,...) 
{
  # The function returns the weights of the minimum-AvDD portfolio. Long-only by hardcode. 
  #
  # R           - Returns data
  
  RC        <- apply(R,2,cumsum)
  
  N <- ncol(RC)
  J <- nrow(RC)
  w <- rep(0, N) ## weights
  u <- rep(0, J) ## high-watermark
  v <- rep(0, J) ## draw downs
  x <- c(w, u, v)
  
  obj <- c(rep(0, N), rep(0, J), rep(1 / J, J))
  
  ## Define constraints
  ## a1: budget constraint
  a1 <- c(rep(1, N), rep(0, 2 * J))
  d1 <- "=="
  b1 <- 1
  ## a2: draw-down constraint (1) assigning summands to v
  a2 <- cbind(-1 * RC, diag(J), -1 * diag(J))
  d2 <- rep("==", J)
  b2 <- rep(0, J)
  ## a3: draw-down constraint (2)
  a3 <- cbind(-1 * RC, diag(J), matrix(0, nrow = J, ncol = J))
  d3 <- rep(">=", J)
  b3 <- rep(0, J)
  ## a4: draw-down constraint (3)
  D1 <- -1.0 * diag(J)
  udiag <- embed(1:J, 2)[, c(2, 1)] 
  D1[udiag] <- 1
  a4 <- cbind(matrix(0, ncol = N, nrow = J), D1, matrix(0, ncol = J, nrow = J))
  a4 <- a4[-J, ]
  d4 <- rep(">=", J-1)
  b4 <- rep(0, J-1)
  ## Combining restrictions
  Amat <- rbind(a1, a2, a3, a4)
  Dvec <- c(d1, d2, d3, d4)
  Bvec <- c(b1, b2, b3, b4)
  
  # bounds on weights
  bounds <- list(lower = list(ind = 1:N, val = rep(wmin,N)),
                 upper = list(ind = 1:N, val = rep(wmax,N)))
  
  ## Solving LP
  result <- Rglpk_solve_LP(obj = obj, mat = Amat, dir = Dvec, rhs = Bvec,
                           max = F, bounds=bounds)
  if(result$status != 0){
    warning(paste("GLPK had exit status:", opt$status))
  }
  
  w       <- as.numeric(result$solution[1:N])
  return(w)
}

## Black-Litterman #########################################################################

# Implemented after Idzorek (2005)

BlackLitterman <- function(R, Views=F, mean.vect, cap.weight, lambda=3.07, tscalar=0.25)
{
  
  N <- ncol(R)
  
  # calculate covariance matrix (Sigma)
  Sigma <- cov(R)
  
  # derive market weights
  if(is.null(cap.weight)) cap.weight <- Max.Sharpe(R)
  
  # derive equlibrium returns
  eq_ret  <- lambda * Sigma %*% cap.weight
  
  if(Views){
    
    # Create View matrix
    P <- diag(N)
    
    Q <- mean.vect # real means
    
    # K = number of views
    K <- nrow(P)
    
    # create storage for variance of views
    var <- 0
    
    for(k in 1:K){
      p <- t(P[k,])
      var[k] <- p %*% Sigma %*% t(p)
    }
    
    Omega <- diag(var, K, K)*tscalar
    
    # calculate the new (posterior) combined return vector
    v1 <- ginv(tscalar * Sigma)
    v2 <- t(P) %*% ginv(Omega) %*% P
    
    v12 <- ginv(v1+v2)
    
    v3 <- (t(P) %*% ginv(Omega) %*% Q) + (v1 %*% eq_ret) 
    
    rhat <- v12 %*% v3
    
    return(as.numeric(rhat))
    
  }else{
    return(as.numeric(eq_ret))
  }
  
}
