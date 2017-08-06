/**
 *   \file linearRegression.hpp
 *   fit models of the form:
 *    y = w0 + w1x
 *
 */

#ifndef LINEARREGRESSION_HPP
#define LINEARREGRESSION_HPP

#include <vector>

double computeAverage(const double *x);

double computeSlope(const double *x , const double *t);

double computeIntercept(const double *x , const double *t);

#endif /* LINEARREGRESSION_HPP */
