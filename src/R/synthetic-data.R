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
nSamples <- 20
inputs <- runif( nSamples )
# =============================================================================
# Targets
# t = sin(2*pi*x) + noise
targets <- sin( 2 * pi * inputs) + rnorm( nSamples, 0, 0.10 )

# t = w0 + w1*x
#w0 = 2.71
#w1 = 0.43
#targets <- w0 + w1 * inputs + rnorm( nSamples, 0, 0.02 )
# =============================================================================
# plot dataset
qplot( inputs, targets ) + geom_smooth(method = "lm")
# =============================================================================
# linear regression in R
linReg <- lm(targets ~ inputs)
# =============================================================================
# write to file
write(inputs, file = "../data/sineData/inputs.txt", ncolumns = 1)
write(targets, file = "../data/sineData/targets.txt", ncolumns = 1)

