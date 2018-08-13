
                                         # Backtesting Machine #

# ------------------------------------------------------------------------------------------------------------ #

# Portfolio Constructor -------------------------------------------------------------------------------------- ##

Portfolio.Constructor <- function(InvestmentP, weights)
{
  # The function calculates and returns the weight changes and returns of the portfolio within a holding perio.d
  #
  # InvestmentP - Returns data
  # weights     - Asset weights
  
  Initial_Investment  <- 100 * weights
  
  New_Investment      <- matrix(NA, nrow=nrow(InvestmentP), ncol=ncol(InvestmentP))
  Investment          <- Initial_Investment * exp(InvestmentP[1, ]/100)
  New_Investment[1,]  <- Investment
  
  for(i in 2: nrow(InvestmentP)){
    New_Investment[i,]  <- Investment * exp(InvestmentP[i, ]/100)
    Investment          <- New_Investment[i,] 
  }
  
  Investment    <- rbind(Initial_Investment, New_Investment) 
  Investment    <- Investment[-nrow(Investment),]
  rownames(Investment) <- c()
  
  Weights      <- Investment/100
  Returns      <- InvestmentP * Weights 
  PortRet      <- apply(Returns, 1, sum) 

  return(list(Weights, PortRet))
}


# Backtest function ------------------------------------------------------------------------------------------ #

BackTester <- function(x, PortOpt, Name, BL=F, Views = F, N.false = 0, Target.Risk=F, Risk_FUN=NULL, ...)
{
  # Basic backtesting function that re-balances portfolio after a specified 
  # holding period.
  #
  # x           - Returns data
  # PortOpt     - Portfolio optimizer function to be used
  # Name        - Strategy name
  # BL          - Estimate expected returns with Black-Litterman ? yes/ no (TRUE/FALSE)
  # Views       - Adjust equlibrium returns with views ? yes/ no (TRUE/FALSE)
  # N.false     - Number of false views (if "0", all views are perfect and vice versa if "12" all views are random)
  # Target.Risk - Maximize portfolio return constraint to a target risk ? yes/ no (TRUE/FALSE)
  # Risk_FUN    - Risk function to be used
  
  rollingWindow<- function(i){
  
    print(i)
    EstimatingP <- x[(StartN - LookBack + i):(StartN - 1 + i),] # Lookback days
    InvestmentP <- x[(StartN + i):(StartN - 1 + HoldPer + i),] # HoldPer days
    Date        <- as.character(Dates[(StartN + i):(StartN - 1 + HoldPer + i)])
    market.w    <- MarketWeights[(StartN - LookBack + i):(StartN - 1 + i),]
    
    if(BL){
      
      cap.weight <- apply(market.w,2,mean)
      
      if(Views){
        mean.vect  <-  apply(InvestmentP, 2, mean)
        
        if(N.false != 0){
          set.seed(1234)
          mean.vect[1:N.false] <- runif(N.false, min= min(mean.vect)-1, max = max(mean.vect)+1)
        } 
      }else{ mean.vect <- NULL}
      
    }else{cap.weight <- NULL}
    
    if(Target.Risk){
      
      w           <- Max.Sharpe(EstimatingP, BL=BL, Views=Views, mean.vect=mean.vect, cap.weight = cap.weight)
      Target      <- Risk_FUN(EstimatingP%*%w) 
      
    }else{Target <- NULL}
    
    w <- PortOpt(EstimatingP, BL=BL, Views=Views, mean.vect=mean.vect, cap.weight = cap.weight, TargetRisk = Target, ...)

    results <- Portfolio.Constructor(InvestmentP, w) # calculates portfolio returns
    weights <- results[[1]]
    PortRet <- results[[2]]
    
    colnames(weights) <- etf.list
    
    if(is.null(Risk_FUN)){
      
      return(list(Date,weights,PortRet))
      
    } else{
      
      OptRisk     <- Risk_FUN(EstimatingP %*% w)
      if(Target.Risk){
        return(list(Date, weights,PortRet, rep(OptRisk,nrow(weights)), rep(Target,nrow(weights)) ) )
      } else{
        return(list(Date, weights,PortRet, rep(OptRisk,nrow(weights)) ) ) 
      }
      
    } 
    
  }
  
  i       <- seq(0,nrow(x)-StartN-HoldPer+1, by=HoldPer)
  result  <- lapply(i, rollingWindow)
  
  Date    <- lapply(result, `[[`, 1)
  weights <- lapply(result, `[[`, 2)
  Returns <- lapply(result, `[[`, 3)

  Date    <- do.call(c,Date)
  weights <- do.call(rbind,weights)
  Returns <- do.call(c,Returns)
  
  if(is.null(Risk_FUN)){
    
    results <- data.frame(Date,weights, Returns)
    colnames(results) <- c("Date",colnames(weights),Name)
    return(results)
    
  }else{
    
    if(Target.Risk){
      Risk_Opt    <- lapply(result, `[[`, 4)
      Risk_Real   <- lapply(result, `[[`, 5)
      
      Risk_Opt    <- do.call(c,Risk_Opt)
      Risk_Real   <- do.call(c,Risk_Real)
      
      results <- data.frame(Date, weights, Returns, Risk_Opt, Risk_Real)
      colnames(results) <- c("Date",colnames(weights),Name, "Optimized_Risk", "Target_Risk")
    } else{
      Risk_Opt    <- lapply(result, `[[`, 4)
      Risk_Opt    <- do.call(c,Risk_Opt)
      
      results <- data.frame(Date, weights, Returns, Risk_Opt)
      colnames(results) <- c("Date",colnames(weights),Name, "Optimized_Risk")
    }
    
    return(results)
  }
  
}

# ------------------------------------------------------------------------------------------------------------ #
