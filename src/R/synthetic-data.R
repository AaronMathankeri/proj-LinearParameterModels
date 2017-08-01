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
# t = sin(2*pi*x) + noise
#targets <- sin( 2 * pi * inputs) + rnorm( nSamples, 0, 0.30 )

# t = w0 + w1*x
w0 = 2.71
w1 = 0.43
targets <- w0 + w1 * inputs + rnorm( nSamples, 0, 0.02 )
# =============================================================================
# plot dataset
qplot( inputs, targets )
# =============================================================================
# linear regression in R
linReg <- lm(targets ~ inputs)
# =============================================================================
# learn these parameters manually
# 1. specify parametric functional form of the model
# 2. determine the values of the parameters using max lnL/error function
# =============================================================================
# 1. y(x,w) = w0 + w1*x
# 2. E(w) = 0.5 * sum_i (y_i - t_i)^2

leastSquaresError <- function( x , w, t){
  error = 0.0;
  for( i in 1:nSamples ){
    y <- w[1] + w[2]*x[i]  
    error = error + ( y - t[i])^2
  }
  error = 0.5 * error
  return( error )
}

# random initialize the weights
weights <- c( 0.2 , 3.9)

myError <- leastSquaresError( inputs, weights, targets)
print( myError )
# =============================================================================
# learn best parameters
