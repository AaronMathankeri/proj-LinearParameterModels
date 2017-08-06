
#include "linearParameterModels.hpp"

void computeOutputs( const std::vector<double> x, const std::vector<double> w , std::vector<double> &y){

      int order = w.size() + 1;
      for (int i = 0; i < x.size( ) ; ++i) {
	    double temp = 0.0;
	    for (int j = 0; j < order; ++j) {
		  temp += w[j]*pow(x[i], j);
	    }
	    y[i] = temp;
      }
}

double leastSquaresError( const std::vector<double> x , const std::vector<double> t, const std::vector<double> w){
      // w only goes up to M--the order of polynomial!
      std::vector<double> y( x.size( ) , 0.0 );
      computeOutputs( x, w, y);

      double error = 0.0;
      for (int i = 0; i < x.size(); ++i) {
	    error += (y[i] * x[i] - t[i]) * (y[i] * x[i] - t[i]);
      }

      error /= 0.5;

      return error;
}

float fRand(float fMin, float fMax){
      float f = (float)rand() / RAND_MAX;
      return fMin + f * (fMax - fMin);
}
void setRandomWeights( std::vector<double> &weights ){
      for (int i = 0; i < weights.size(); ++i) {
	    float temp = fRand( -10.0, 10.0);
	    weights[i] = temp;
      }
}
