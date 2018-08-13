
                         # Cockpit for Portfolio Optimization #

# ------------------------------------------------------------------------------------------------------ #

# Packages and sourcing files -------------------------------------------------------------------------- #

library(Rglpk)       # solves linear program (Expected Shortfall Optimization)
library(DEoptim)     # solves non-linear program (VaR, Maximum Drawdown Optimization)
library(quadprog)    # solves quadratic program (Mean-Variance Optimization)
library(doMC)        # parallel computing -> speed up non-linear optimization
#library(doSNOW)     # for windows user: parallel computing -> speed up non-linear optimization
library(inline)      # kill parallel processing 
library(randtoolbox) # create sobol series -> starting values for non-linear optimization
library(MASS)        # for matrix algebra
library(ggplot2)     # for illustration
library(reshape2)    # for illustration
library(xtable)      # creates LaTex tables

wd <- "~/Desktop/.../" # <---- paste file path here
setwd(wd)

source(paste(wd,"PortfolioOpt_Machine.R",sep="")) # Optimization functions
source(paste(wd,"Backtesting_Machine.R",sep=""))  # Backtest functions
source(paste(wd,"Risk_FUN.R",sep=""))  # Risk functions
source(paste(wd,"PerformanceA_FUN.R",sep="")) # Performance Analysis functions

# Prep to do parallel computing ------------------------------------------------------------------------ #

# set this appropriate to your system
options(cores=3) # number of cores to use 

# function calls "waitpid" - stops parallel processing
# Source: https://stackoverflow.com/questions/25388139/r-parallel-computing-and-zombie-processes
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

# Data & parameter ------------------------------------------------------------------------------------- #

load("Data.rdata")
load("MarketWeights.rdata")

# Date and Assets
Dates      <- Data[,1]
Cash       <- Data[,2]*100
etf.logrtn <- as.matrix(Data[,3:14]*100)
etf.list   <- colnames(etf.logrtn)

# Market weights for B-L
MarketWeights <- MarketWeights[,2:13]

#  Parameter for rolling-window Backtest 
LookBack    <- 252    # Number of days used for asset allocation estimation (252 trading days ~ 1 year)
HoldPer     <- 21     # Number of days asset allocation is hold (21 trading days ~ 1 month)
StartN      <- LookBack + 1 # Start row

# Make Backtest reproducable
set.seed(1234) # needed because of randomness in non-linear optimization 

# Risk optimization parameter
alpha       <- 0.05       # VaR alpha
beta        <- 0.95       # CVaR confidence interval
wmin        <- 0          # Min. allowed weight per instrument

# Maximum weight constraint
wmax        <- 1       # Max. allowed weight per instrument
#wmax        <- 0.35    # Max. allowed weight per instrument

# Start Backtest  -------------------------------------------------------------------------------------- #

### No optimization
OneOverN <- BackTester(etf.logrtn, OneOverN.Optimizer,"OneOverN P")

### Mean-Variance
MaxSharpe <- BackTester(etf.logrtn, Max.Sharpe ,"MaxSharpe P")

## Minimum Risk Optimization
# Quadratic programming
GlobalMin.Var   <- BackTester(etf.logrtn, GlobalMinVar,"GlobalMinVar P")

# Linear programming
MinCVaR          <- BackTester(etf.logrtn, Min.CVaR, "MinCVaR", Risk_FUN = CVaR_FUN)
MinAvDD          <- BackTester(etf.logrtn, Min.AvDD, "MinAvDD", Risk_FUN = AvDD_FUN) 

# Non-linear programming
# **** NOTE!!! execution of "MinVaR" and "MinMDD" may take ~ 40 min (each) ***** 
Min.VaR         <- BackTester(etf.logrtn, Min.VaR, "MinVaR P", Risk_FUN = VaR_FUN)
MinMDD         <- BackTester(etf.logrtn, Min.MDD, "MinMDD", Risk_FUN = MDD_FUN)

### Target Risk Optimization
TargetCVaR        <- BackTester(etf.logrtn, Target.CVaR, "TargetCVaR", Target.Risk=T, Risk_FUN = CVaR_FUN)
# **** NOTE!!! execution of "TargetVaR" may take ~ 40 min ***** 
TargetVaR       <- BackTester(etf.logrtn, Target.VaR, "TargetVaR", Target.Risk=T, Risk_FUN = VaR_FUN)
TargetMDD       <- BackTester(etf.logrtn, Target.MDD, "TargetMDD", Target.Risk=T, Risk_FUN = MDD_FUN)
TargetAvDD      <- BackTester(etf.logrtn, Target.AvDD, "TargetAvDD", Target.Risk=T, Risk_FUN = AvDD_FUN)

### Target Risk + BL
# Without views
MaxSharpe       <- BackTester(etf.logrtn, Max.Sharpe ,"MaxSharpe P", BL=T)
TargetCVaR      <- BackTester(etf.logrtn, Target.CVaR, "TargetCVaR", BL=T, Target.Risk=T, Risk_FUN = CVaR_FUN)
# **** NOTE!!! execution of "TargetVaR" may take ~ 40 min ***** 
TargetVaR       <- BackTester(etf.logrtn, Target.VaR, "TargetVaR", BL=T, Target.Risk=T, Risk_FUN = VaR_FUN)
TargetMDD       <- BackTester(etf.logrtn, Target.MDD, "TargetMDD", BL=T, Target.Risk=T, Risk_FUN = MDD_FUN)
TargetAvDD      <- BackTester(etf.logrtn, Target.AvDD, "TargetAvDD", BL=T, Target.Risk=T, Risk_FUN = AvDD_FUN) 

# With perfect views
TargetCVaR        <- BackTester(etf.logrtn, Target.CVaR, "TargetCVaR", BL=T, Views=T, Target.Risk=T, Risk_FUN = CVaR_FUN)
# **** NOTE!!! execution of "TargetVaR" may take ~ 40 min ***** 
TargetVaR       <- BackTester(etf.logrtn, Target.VaR, "TargetVaR", BL=T, Views=T, Target.Risk=T, Risk_FUN = VaR_FUN)
TargetMDD       <- BackTester(etf.logrtn, Target.MDD, "TargetMDD", BL=T, Views=T, Target.Risk=T, Risk_FUN = MDD_FUN)
TargetAvDD      <- BackTester(etf.logrtn, Target.AvDD, "TargetAvDD", BL=T, Views=T, Target.Risk=T, Risk_FUN = AvDD_FUN) 

# With random views
TargetCVaR        <- BackTester(etf.logrtn, Target.CVaR, "TargetCVaR", BL=T, Views=T, N.false=12, Target.Risk=T, Risk_FUN = CVaR_FUN)
# **** NOTE!!! execution of "TargetVaR" may take ~ 40 min ***** 
TargetVaR       <- BackTester(etf.logrtn, Target.VaR, "TargetVaR", BL=T, Views=T, N.false=12, Target.Risk=T, Risk_FUN = VaR_FUN)
TargetMDD       <- BackTester(etf.logrtn, Target.MDD, "TargetMDD", BL=T, Views=T, N.false=12, Target.Risk=T, Risk_FUN = MDD_FUN)
TargetAvDD      <- BackTester(etf.logrtn, Target.AvDD, "TargetAvDD", BL=T, Views=T, N.false=12, Target.Risk=T, Risk_FUN = AvDD_FUN) 

# Save Results  -------------------------------------------------------------------------------------- #

## Save results as "rdata"
save(OneOverN, file="OneOverN.rdata")
save(MaxSharpe, file="MaxSharpe.rdata")
save(GlobalMin.Var, file="GlobalMin_Var.rdata")
save(MinCVaR, file="Min_CVaR.rdata")
save(Min.VaR , file="Min_VaR.rdata")
save(MinAvDD , file="Min_AvDD.rdata")
save(TargetCVaR , file="Target_CVaR.rdata")
save(TargetVaR , file="Target_VaR.rdata")
save(TargetMDD , file="Target_MDD.rdata")
save(TargetAvDD , file="Target_AvDD.rdata")
# load results using: load("FileNameExample.rdata")

## Save results as tables
# Traditional MVO
results1 <- data.frame(OneOverN,MaxSharpe)
write.table(results1, "MVO.txt", sep="\t", row.names = FALSE)

# Minimum Risk
results2 <- data.frame(OneOverN, GlobalMin.Var,MinCVaR, Min.VaR,MinAvDD)
write.table(results2, "MinRisk.txt", sep="\t", row.names = FALSE)

# Target Risk
results3 <- data.frame(OneOverN, MaxSharpe, TargetCVaR, TargetVaR, TargetMDD, TargetAvDD)
write.table(results3, "TargetRisk.txt", sep="\t", row.names = FALSE)

# Performance Analysis  ----------------------------------------------------------------------------------- #

## Traditional MVO
# Performance
Naive <- Statistics(OneOverN, "Naive")
MVO <- Statistics(MaxSharpe, "MVO")
MVO_con <- Statistics(MaxSharpe, "Constrained MVO")
results <- rbind(Naive, MVO, MVO_con)
xtable(results)

# Risk
Naive <- RiskStats(OneOverN, "Naive")
MVO <- RiskStats(MaxSharpe, "MVO")
MVO_con <- RiskStats(MaxSharpe, "Constrained MVO")
results <- rbind(Naive, MVO, MVO_con)
xtable(results)

## Maximum-return with risk limit strategies
# Performance
Naive <- Statistics(OneOverN, "Naive")
MVO <- Statistics(MaxSharpe, "MVO")
VaR <- Statistics(TargetVaR, NameR="VaR")
CVaR <- Statistics(TargetCVaR, NameR="CVaR")
MDD <- Statistics(TargetMDD, NameR="MDD")
AvDD <- Statistics(TargetAvDD, NameR="AvDD") 
results <- rbind(Naive, MVO, VaR, CVaR, MDD, AvDD)
xtable(results)

# Risk
Naive <- RiskStats(OneOverN, "Naive")
MVO <- RiskStats(MaxSharpe, "MVO")
VaR <- RiskStats(TargetVaR, NameR="VaR")
CVaR <- RiskStats(TargetCVaR, NameR="CVaR")
MDD <- RiskStats(TargetMDD, NameR="MDD")
AvDD <- RiskStats(TargetAvDD, NameR="AvDD") 
results <- rbind(Naive, MVO, VaR, CVaR, MDD, AvDD)
xtable(results)

## Minimum-risk strategies
# Performance
Naive <- Statistics(OneOverN, "Naive")
MVO <- Statistics(GlobalMin.Var, "Variance")
VaR <- Statistics(Min.VaR, NameR="VaR")
CVaR <- Statistics(MinCVaR, NameR="CVaR")
MDD <- Statistics(MinMDD, NameR="MDD")
AvDD <- Statistics(MinAvDD, NameR="AvDD") 
results <- rbind(Naive, MVO, VaR, CVaR, MDD, AvDD)
xtable(results)

# Risk
Naive <- RiskStats(OneOverN, "Naive")
MVO <- RiskStats(GlobalMin.Var, "Variance")
VaR <- RiskStats(Min.VaR, NameR="VaR")
CVaR <- RiskStats(MinCVaR, NameR="CVaR")
MDD <- RiskStats(Min.MDD, NameR="MDD")
AvDD <- RiskStats(MinAvDD, NameR="AvDD") 
results <- rbind(Naive, MVO, VaR, CVaR, MDD, AvDD)
xtable(results)

# ------------------------------------------------------------------------------------------------------ #
