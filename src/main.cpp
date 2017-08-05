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

int main(int argc, char *argv[])
{
      cout << "Aaron's Back." << endl;

      const int nPatterns = 20;
      vector<double> x( nPatterns, 0.0 );
      vector<double> t( nPatterns, 0.0 );
      vector<double> weights( 2, 0.0 );
      string inputsFile = "./data/sineData/inputs.txt";
      string targetsFile = "./data/sineData/targets.txt";

      loadData( x , inputsFile );
      loadData( t , targetsFile );

      cout << "Inputs" << endl;
      printVector( x );

      cout << "Targets" << endl;
      printVector( t );
      /*
      cout << "\nComputing Line of Best Fit ..." << endl;
      weights[1] = computeSlope( x , t);
      weights[0] = computeIntercept( x , t);
      
      cout << "w0 = " << weights[0] << endl;
      cout << "w1 = " << weights[1] << endl;
      */
      return 0;
}

