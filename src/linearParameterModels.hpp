#ifndef LINEARPARAMETERMODELS_H
#define LINEARPARAMETERMODELS_H

#include <vector>
#include <math.h>
#include <stdlib.h>

void computeOutputs( const std::vector<double> x, const std::vector<double> w , std::vector<double> &y);

double leastSquaresError( const std::vector<double> x , const std::vector<double> t, const std::vector<double> w);

float fRand(float fMin, float fMax);

void setRandomWeights( std::vector<double> &weights );

#endif /* LINEARPARAMETERMODELS_H */
