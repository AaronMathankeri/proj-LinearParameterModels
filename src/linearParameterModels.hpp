#ifndef LINEARPARAMETERMODELS_H
#define LINEARPARAMETERMODELS_H

#include <vector>
#include <math.h>
#include <stdlib.h>
#include "mkl.h"
#include "parameters.hpp"

void computeDesignMatrix( const double *x, double *Phi );
void computePseudoInverse( const double* Phi , double *phiPsuedoInverse );
void solveNormalEquations( const double *inversePhi, const double *t, double *w );
void computeOutputs( const double *x, const double *w , double *y );
double computeLeastSquaresError( const double *t, const double *y );
void setRandomWeights( double *weights );
double fRand( const double fMin, const double fMax);
//---------------------------
// older functions
/*
void computeOutputs( const std::vector<double> x, const std::vector<double> w , std::vector<double> &y);
double leastSquaresError( const std::vector<double> x , const std::vector<double> t, const std::vector<double> w);
float fRand(float fMin, float fMax);
void setRandomWeights( std::vector<double> &weights );
*/
#endif /* LINEARPARAMETERMODELS_H */
