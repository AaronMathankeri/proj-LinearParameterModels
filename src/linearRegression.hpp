/**
 *   \file linearRegression.hpp
 *   fit models of the form:
 *    y = w0 + w1x
 *
 */

#ifndef LINEARREGRESSION_HPP
#define LINEARREGRESSION_HPP

#include <vector>

double computeAverage(const std::vector<double> x);

double computeSlope(const std::vector<double> x , const std::vector<double> t);

double computeIntercept(const std::vector<double> x , const std::vector<double> t);

#endif /* LINEARREGRESSION_HPP */
