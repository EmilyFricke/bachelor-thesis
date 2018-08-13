
                        # Functions for Performance Analysis #

# --------------------------------------------------------------------------------------------- #

# Load risk functions ------------------------------------------------------------------------- #

wd <- "~/Desktop/../" # <---- paste file path here
setwd(wd)
source(paste(wd,"Risk_FUN.R",sep="")) 

# Set risk function parameter ----------------------------------------------------------------- #

alpha <- 0.05 # (1-alpha)= confidence interval ---> VaR function
beta <- 0.95 # beta = confidence interval ---> VaR function

# Cummulated return - function ---------------------------------------------------------------- #

CumReturn <- function(R){
  R <- R[,14]
  CumReturn <- tail(cumsum(R), 1)
  return(round(CumReturn,2))
}

# Annualized mean return - function ----------------------------------------------------------- #

MeanFUN <- function(R){
  R <- R[,14]
  return(round(mean(R)*252, 2))
}

# Aannualized standard deviation - function --------------------------------------------------- #

SdFUN <- function(R){
  R <- R[,14]
  SD <- sd(R)*sqrt(252)
  return(round(SD, 2))
}

# Annualized Sharpe Ratio - function ----------------------------------------------------------- #

sharpeRatio<-function(R, Rfree = Cash)
{
  R <- R[,14] 
  Rfree <- mean(Cash[-(1:252)])
  
  Sharpe <- (mean(R)-Rfree)/sd(R)
  return(round(Sharpe*sqrt(252),2))
}

# Annualized Sortino Ratio - function ---------------------------------------------------------- #

sortinoRatio<-function(R, Rfree = Cash)
{
  R <- R[,14] 
  Rfree <- mean(Cash[-(1:252)])
  
  Sortino <- (mean(R)-Rfree)/sd(pmin(R,0))
  return(round(Sortino*sqrt(252),2))
}

# Annualized turnover - function ---------------------------------------------------------------- #

Turnover <- function(R)
{
  w_matrix <- R[,2:13]
  
  Jdifference <- 0
  Tdifference <- 0
  
  # rebalancing all 21 days
  for(t in seq(21,nrow(R)-21, by=21)){
    
    for(j in 1:ncol(w_matrix)){ Jdifference[j] <- abs(w_matrix[t+1,j]-w_matrix[t,j]) }
    
    Tdifference[t] <- sum(Jdifference)
  }
  
  difference <- sum( na.omit(Tdifference))
  
  turnover <- 1/nrow(w_matrix) * difference * 252 * 100
  
  return(round(turnover,2))
  
}

# Sum of squared portfolio weights (SSPW) - function -------------------------------------------- #

Concentration <- function(R){
  
  w_matrix <- R[,2:13]
  
  sum_ws_squared <- 0
  
  for(i in 1:nrow(R)){
    ws                <- w_matrix[i,]
    ws_squared        <- ws^2
    sum_ws_squared[i] <- sum(ws_squared)
  }
  
  concentration <- mean(sum_ws_squared)*100
  
  return(round(concentration,2))
}

# Creates LaTex table output (1/2) --------------------------------------------------------------- #

Statistics <- function(R, NameR){
  
  CompR <- CumReturn(R)
  MeanR <- MeanFUN(R)
  SRR   <- sharpeRatio(R)
  SR    <- sortinoRatio(R)
  TR    <- Turnover(R)
  CR    <- Concentration(R)
  
  results <- data.frame(CompR, MeanR, SRR, SR, TR, CR)
  colnames(results) <- c("CR_P", "mu_P", "SR_P", "S_P", "Turnover", "Concentration")
  rownames(results) <- NameR

  return(results)
}

# Creates LaTex table output (2/2) ---------------------------------------------------------------- #

RiskStats <- function(R, NameR){
  
  SdR   <- SdFUN(R)
  
  R <- R[,14]
  VAR   <- VaR_FUN(R)
  CVaR    <- CVaR_FUN(R)
  MDD   <- MDD_FUN(R)
  AvDD  <- AvDD_FUN(R)
  
  results <- data.frame(SdR, round(VAR,2), round(CVaR,2), round(MDD,2), round(AvDD,2))
  colnames(results) <- c("SD", "VaR_95%", "CVaR_95%", "MDD_10y", "AvDD_10y")
  rownames(results) <- NameR
  
  return(results)

} 

