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
#include "linearRegression.hpp"
#include "mkl.h"

using namespace std;

#define NUM_PATTERNS 10

void printVector( const double *x , const int length ){
      for (int i = 0; i < length; i++) {
	    printf ("%12.3f", x[i]);
	    printf ("\n");
      }
}

void loadData( double *x , string fileName ){
      ifstream file  ( fileName );
      if(file.is_open()) {
	    for (int i = 0; i < NUM_PATTERNS; ++i) {
		  file >> x[i];
	    }
      }
}

int main(int argc, char *argv[])
{
      cout << " Aaron's Back." << endl;

      string inputsFile = "./data/linearRegression/inputs.txt";
      string targetsFile = "./data/linearRegression/targets.txt";
      
      double *x, *t, *w;
      int order = 1;
      
      x = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      t = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      w = (double *)mkl_malloc( (order + 1)*sizeof( double ), 64 );

      memset( x, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( t, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( w, 0.0,  (order + 1) * sizeof(double));      

      loadData( x , inputsFile );
      loadData( t , targetsFile );
      
      cout << "Inputs" << endl;
      printVector( x , NUM_PATTERNS );

      cout << "Targets" << endl;
      printVector( t , NUM_PATTERNS );

      cout << "\nComputing Line of Best Fit ..." << endl;
      w[1] = computeSlope( x , t);
      w[0] = computeIntercept( x , t);
      cout << " w0 = " << w[0] << endl;
      cout << " w1 = " << w[1] << endl;

      printf ("\n Deallocating memory \n\n");
      mkl_free( x );
      mkl_free( t );
      mkl_free( w );
      printf (" Example completed. \n\n");

      return 0;
}

/*
//LPModels -- polynomial
  const int M = 3; //polynomial order
  vector<double> weights( M + 1, 0.0 );
  setRandomWeights( weights );

  cout<< "Weights" << endl;
  printVector( weights );

  cout << "Error in estimation is = " << leastSquaresError( x, t, weights ) << endl;
*/

