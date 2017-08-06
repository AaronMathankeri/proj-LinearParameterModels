#include "linearRegression.hpp"
#define NUM_PATTERNS 10

double computeAverage(const double *x){
      double average = 0.0;

      for (int i = 0; i < NUM_PATTERNS; ++i) {
	    average += x[i];
      }

      average /= NUM_PATTERNS;
      return average;
}

double computeSlope(const double *x , const double *t){
      double slope = 0.0;
      double xbar = computeAverage( x );
      double tbar = computeAverage( t );

      double top = 0.0;
      double bottom = 0.0;

      for (int i = 0; i < NUM_PATTERNS; ++i) {
	    top += ( x[i] - xbar) * ( t[i] - tbar) ;
	    bottom += ( x[i] - xbar) * ( x[i] - xbar);
      }

      slope = top/bottom;
      return slope;
}

double computeIntercept(const double *x , const double *t){
      double xbar = computeAverage( x );
      double tbar = computeAverage( t );
      double m = computeSlope( x , t);

      double intercept = 0.0;

      intercept = tbar - m * xbar;

      return intercept;
}
