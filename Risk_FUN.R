
                         # Functions for Performance Analysis #

# --------------------------------------------------------------------------------------------- #

# Value at Risk (VaR) ------------------------------------------------------------------------- #

VaR_FUN <- function(R)
{
  # Function that calculates empiric VaR(alpha)
  #
  # R = return vector
  # type 2 = Inverse of empirical distribution function + averaging at discontinuities
  # see ?quantile() for details
  
  return(as.numeric(quantile(R, p=alpha, type=2))) 
}

# Conditional Value at Risk (CVaR) ------------------------------------------------------------- #

CVaR_FUN <- function(R)
{
  # Function that calculates Expected Shortfall - CVaR(alpha)
  #
  # R = return vector
  
  R <- as.numeric(R)
  
  VaR <- quantile(R,p=alpha)
  n <- length(R)
  
  # average excess loss of all scenarios
  Excess_L <- -(1/n) * sum(pmax(VaR-R,0))
  
  # propability of excess loss
  p <- (1-beta)
  
  # average loss if loss occurs
  ES <- VaR + Excess_L/p
  
  return(as.numeric(ES))  
}

# Maximum Drawdown (MDD) ------------------------------------------------------------------------ #

MDD_FUN    <- function(R) {
  
  # Function that calculates Maximum Drawdown
  #
  # R = return vector
  
  cum.pnl  <- c(0, cumsum(R))
  drawdown <- cum.pnl - cummax(cum.pnl)
  drawdown <- tail(drawdown, -1)
  return(min(drawdown))
}

# Average Drawdown (CDD) ------------------------------------------------------------------------ #

AvDD_FUN <- function(R) {
  # Function that calculates Average Drawdown
  #
  # R = return vector
  
  cum.pnl  <- c(0, cumsum(R))
  dd <- cum.pnl - cummax(cum.pnl)
  dd <- tail(dd, -1)
  
  return(mean(dd))
}
