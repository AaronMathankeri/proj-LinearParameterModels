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
#include <fstream>
#include <sstream>
#include <math.h>
#include "linearRegression.hpp"

using namespace std;

void printVector( const vector<double> x ){
      for ( auto i : x){
	    cout << i << endl;
      }
      cout << endl;
}

void loadData( vector<double> &x , string fileName ){
      ifstream file  ( fileName );
      if(file.is_open()) {
	    for (int i = 0; i < x.size(); ++i) {
		  file >> x[i];
	    }
      }
}

void computeOutputs( const vector<double> x, const vector<double> w , vector<double> &y){

      int order = w.size() + 1;
      for (int i = 0; i < x.size( ) ; ++i) {
	    double temp = 0.0;
	    for (int j = 0; j < order; ++j) {
		  temp += w[j]*pow(x[i], j);
	    }
	    y[i] = temp;
      }
}

double leastSquaresError( const vector<double> x , const vector<double> t, const vector<double> w){
      // w only goes up to M--the order of polynomial!
      vector<double> y( x.size( ) , 0.0 );
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
void setRandomWeights( vector<double> &weights ){
      for (int i = 0; i < weights.size(); ++i) {
	    float temp = fRand( -10.0, 10.0);
	    weights[i] = temp;
      }
}


int main(int argc, char *argv[])
{
      cout << "Aaron's Back." << endl;

      const int nPatterns = 20;
      vector<double> x( nPatterns, 0.0 );
      vector<double> t( nPatterns, 0.0 );

      string inputsFile = "./data/sineData/inputs.txt";
      string targetsFile = "./data/sineData/targets.txt";

      loadData( x , inputsFile );
      loadData( t , targetsFile );

      cout << "Inputs" << endl;
      printVector( x );

      cout << "Targets" << endl;
      printVector( t );

      const int M = 3; //polynomial order
      vector<double> weights( M + 1, 0.0 );
      setRandomWeights( weights );

      cout<< "Weights" << endl;
      printVector( weights );

      cout << "Error in estimation is = " << leastSquaresError( x, t, weights ) << endl;

      //perform optimization

      return 0;
}

/*
  cout << "\nComputing Line of Best Fit ..." << endl;
  weights[1] = computeSlope( x , t);
  weights[0] = computeIntercept( x , t);
      
  cout << "w0 = " << weights[0] << endl;
  cout << "w1 = " << weights[1] << endl;
*/
