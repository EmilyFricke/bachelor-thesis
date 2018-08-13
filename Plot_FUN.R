
                                    # Plot  Functions #

# ------------------------------------------------------------------------------------------------------- #

# Packages and sourcing files --------------------------------------------------------------------------- #

library(ggcorrplot)
library(extrafont)

wd <- "~/Desktop/.../"  # <---- paste file path here
setwd(wd)

source(paste(wd,"Risk_FUN.R",sep="")) 

# Data -------------------------------------------------------------------------------------------------- #

load('Data.rdata')

Date       <- Data[,1]
Cash       <- Data[,2]
etf.logrtn <- Data[,3:14]*100
etf.list   <- colnames(etf.logrtn)
N          <- ncol(etf.logrtn)


# Correlation matrix ------------------------------------------------------------------------------------ #

# Compute a correlation matrix
corr <- round(cor(etf.logrtn), 1)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(etf.logrtn)

# plot correlation matrix
ggcorrplot(corr, hc.order = TRUE, type = "lower", outline.col = "white",
           ggtheme = ggplot2::theme_gray, colors = c("black", "white", "darkred"),
           p.mat = p.mat) +
  theme(axis.text.x = element_text(family="Times", size=10, colour = "black"),
        axis.text.y = element_text(family="Times", size=10, colour = "black"),
        legend.title = element_text(family="Times", size=10, colour = "black"),
        legend.text = element_text(family="Times", size=10, colour = "black"))

# Create dataframe for ggplot function
Data[,3:14] <- Data[,3:14]*100
ggplot_data <- melt(Data, id = "Date", value.name = "Return", variable.name ="Asset")

ggplot(ggplot_data, aes(x = Asset, y = Return)) +
  geom_boxplot()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(family="Times", size=10, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(family="Times", size=10, colour = "black"))


# Boxplots ---------------------------------------------------------------------------------------------- #

# Create dataframe for ggplot function
Data[,3:14] <- Data[,3:14]*100
ggplot_data <- melt(Data, id = "Date", value.name = "Return", variable.name ="Asset")

# Plot Boxplots
ggplot(ggplot_data, aes(x = Asset, y = Return)) +
  geom_boxplot()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(family="Times", size=10, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(family="Times", size=10, colour = "black"))


# Return Distributions ---------------------------------------------------------------------------------- #

# Create dataframe for ggplot function
N <- nrow(Data)
Data[,3:14] <- Data[,3:14]*100
ggplot_data <- melt(Data, id = "Date", value.name = "Return", variable.name ="Asset")

HistPlots <- function(i){
  
  ggplot_data2 <- ggplot_data[(N+1+i):(N+N+1+i),]
  
  p <- ggplot(ggplot_data2, aes(x=ggplot_data2$Return)) + 
    geom_histogram(bins = 100, aes(y=..density..)) + 
    stat_function(fun = dnorm, color="red", args = list(mean = mean(ggplot_data2$Return, na.rm = TRUE) , sd = sd(ggplot_data2$Return, na.rm = TRUE)))+
    scale_x_continuous(limits = c(-6,6))+ 
    theme(axis.title.y = element_blank()) + xlab("(%)") + ggtitle(ggplot_data2$Asset)+
    theme(axis.text.x = element_text(family="Times", size=11, colour = "black"),
          axis.text.y = element_text(family="Times", size=11, colour = "black"),
          axis.title.x= element_text(family="Times", size=11, colour = "black"),
          plot.title = element_text(family="Times", size=14, colour = "black"))
  
  return(p)
  
}

i <- seq(0,(nrow(ggplot_data)-2*N), by=N)
x <- lapply(i, HistPlots)
ml <- marrangeGrob(x, nrow=4, ncol=3)
ggsave("ReturnDistributions.pdf", ml, width=8.27, height=11.69)

# ------------------------------------------------------------------------------------------------------- #
