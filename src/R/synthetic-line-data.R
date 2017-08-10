# =============================================================================
# Clear workspace
rm(list=ls())

set.seed(123)
# =============================================================================
setwd('/Users/aaron/playground/proj-ML-gallery/data/')
# =============================================================================
# libraries
library(ggplot2)
# =============================================================================
# sample from uniform distribution
# x ~ U(0,1)
nSamples <- 10
inputs <- runif( nSamples )
# =============================================================================
# Targets
# t = w0 + w1*x
w0 = 2.71
w1 = 0.43
targets <- w0 + w1 * inputs + rnorm( nSamples, 0, 0.02 )
# =============================================================================
# linear regression in R
linReg <- lm(targets ~ inputs)
# =============================================================================
# create dataframe for plots
df <- cbind.data.frame( inputs, targets )
# =============================================================================
# define functions to plot with data
f1 <- function(x) {w0 + w1*x}
f2 <- function(x) {2.72 + 0.403*x} # found from running c++ file
# =============================================================================
p <- ggplot(df, aes(x = inputs, y = targets)) +
  geom_point()

p + stat_function(fun = f1, aes( colour = "true")) +
  stat_function(fun = f2, aes(colour = "prediction")) +
  scale_colour_manual("Functions",values = c("red", "blue"), labels = c("true", "prediction"))
# =============================================================================
# write to file
write(inputs, file = "../data/linearRegression/inputs.txt", ncolumns = 1)
write(targets, file = "../data/linearRegression/targets.txt", ncolumns = 1)

