#include "linearRegression.hpp"

double computeAverage(const std::vector<double> x){
      double average = 0.0;

      for (int i = 0; i < x.size(); ++i) {
	    average += x[i];
      }
      average /= x.size();

      return average;
}

double computeSlope(const std::vector<double> x , const std::vector<double> t){
      double slope = 0.0;
      double xbar = computeAverage( x );
      double tbar = computeAverage( t );

      double top = 0.0;
      double bottom = 0.0;

      for (int i = 0; i < x.size(); ++i) {
	    top += ( x[i] - xbar) * ( t[i] - tbar) ;
	    bottom += ( x[i] - xbar) * ( x[i] - xbar);
      }

      slope = top/bottom;
      return slope;
}

double computeIntercept(const std::vector<double> x , const std::vector<double> t){
      double xbar = computeAverage( x );
      double tbar = computeAverage( t );
      double m = computeSlope( x , t);

      double intercept = 0.0;

      intercept = tbar - m * xbar;

      return intercept;
}
