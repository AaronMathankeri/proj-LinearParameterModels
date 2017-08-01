//Copyright (C) 2017 Aaron Mathankeri (aaron@quantifiedag.com)
//License: QuantifiedAG Software License.  See LICENSE.txt for full license.

/*
 *   \file example.cpp
 *   \brief A Documented file.
 *
 *  Detailed description
 *
 */
#include <iostream>
#include <vector>
using namespace std;

void computeOutputs( double* x , double * w){
      cout << "Computing Outputs" << "\n";
}

void computeleastSquaresError( double *x , double *w, double *t){
      cout <<"Calculating Error" << endl;;

      computeOutputs( x , w );
}

void printVector( vector<double> x ){
      for ( auto i : x){
	    cout << i << " ";
      }
      cout << endl;
}

double computeAverage( vector<double> x){
      double average = 0.0;

      for (int i = 0; i < x.size(); ++i) {
	    average += x[i];
      }
      average /= x.size();

      return average;
}

double computeSlope( vector<double> x , vector<double> t){
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

double computeIntercept( vector<double> x , vector<double> t){
      double xbar = computeAverage( x );
      double tbar = computeAverage( t );
      double m = computeSlope( x , t);

      double intercept = 0.0;

      intercept = tbar - m * xbar;

      return intercept;
}

int main(int argc, char *argv[])
{
      cout << "Aaron's Back." << endl;

      //      double *x,*w,*t;

      //      computeleastSquaresError( x , w , t);
      const int nPatterns = 10;
      vector<double> x{0.7883051, 0.4089769, 0.8830174, 0.9404673, 0.0455565,
		  0.5281055, 0.8924190, 0.5514350, 0.4566147, 0.9568333}; 

      vector<double> t{3.046626, 2.889522, 3.115309, 3.079856, 2.763393,
		  2.947162, 3.144307, 2.958099, 2.911109, 3.100460};

      vector<double> weights{0.2 , 3.9}; 

      cout << "Inputs" << endl;
      printVector( x );

      cout << "Targets" << endl;
      printVector( t );

      cout << "Initial Weights" << endl;
      printVector( weights );

      cout << "\nComputing Line of Best Fit ..." << endl;
      weights[1] = computeSlope( x , t);
      weights[0] = computeIntercept( x , t);
      
      cout << "w0 = " << weights[0] << endl;
      cout << "w1 = " << weights[1] << endl;

      
      return 0;
}

