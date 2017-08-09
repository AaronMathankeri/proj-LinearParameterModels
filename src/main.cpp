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
#include "mkl.h"
#include "ioFunctions.hpp"
#include "parameters.hpp"
#include "linearParameterModels.hpp"

using namespace std;

int main(int argc, char *argv[])
{
      cout << " Aaron's Back." << endl;

      //--------------------------------------------------------------------------------
      // declare variables for calculations
      double *x, *t, *w, *y, *designMatrix, *designPseudoInverse;

      x = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      t = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      w = (double *)mkl_malloc( (ORDER)*sizeof( double ), 64 );
      y = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      designMatrix = (double *)mkl_malloc( NUM_PATTERNS*ORDER*sizeof( double ), 64 );
      designPseudoInverse = (double *)mkl_malloc( ORDER*NUM_PATTERNS*sizeof( double ), 64 );
      
      memset( x, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( t, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( w, 0.0,  ORDER* sizeof(double));
      memset( y, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( designMatrix, 0.0, NUM_PATTERNS * ORDER* sizeof(double));     
      memset( designPseudoInverse, 0.0,  ORDER*NUM_PATTERNS * sizeof(double));      
      //--------------------------------------------------------------------------------
      //read data
      string inputsFile = "./data/linearRegression/inputs.txt";
      string targetsFile = "./data/linearRegression/targets.txt";

      loadData( x , inputsFile );
      loadData( t , targetsFile );
      
      cout << "Inputs" << endl;
      printVector( x, NUM_PATTERNS );

      cout << "Targets" << endl;
      printVector( t, NUM_PATTERNS );
      //--------------------------------------------------------------------------------
      cout << "Randomly initialize weights:" << endl;
      setRandomWeights( w );
      printVector( w, ORDER );
      /*
      computeOutputs( x, w , y);
      cout << "\nOutputs" << endl;
      printVector( y , NUM_PATTERNS );
      cout << "Least Squares Error :" << computeLeastSquaresError( t , y ) << endl;
      */
      //--------------------------------------------------------------------------------
      printf ("\n Design Matrix: \n");
      computeDesignMatrix( x, designMatrix );
      printMatrix( designMatrix, NUM_PATTERNS, ORDER );
      //--------------------------------------------------------------------------------
      // compute moore-penrose psuedo inverse of design matrix
      computePseudoInverse( designMatrix , designPseudoInverse );
      //--------------------------------------------------------------------------------
      cout << "\nSolving Normal Equations ..." << endl;
      solveNormalEquations( designPseudoInverse, t , w);

      cout << " w0 = " << w[0] << endl;
      cout << " w1 = " << w[1] << endl;
      //--------------------------------------------------------------------------------
      // compute outputs
      computeOutputs( x, w , y);
      cout << "\nOutputs" << endl;
      printVector( y , NUM_PATTERNS );
      //--------------------------------------------------------------------------------
      cout << "Least Squares Error :" << computeLeastSquaresError( t , y ) << endl;
      //--------------------------------------------------------------------------------
      printf ("\n Deallocating memory \n\n");
      mkl_free( x );
      mkl_free( t );
      mkl_free( w );
      mkl_free( y );
      mkl_free( designMatrix );
      mkl_free( designPseudoInverse );
      printf (" Example completed. \n\n");

      return 0;
}

/*
//LPModels -- polynomial
  const int M = 3; //polynomial ORDER
  vector<double> weights( M + 1, 0.0 );
  setRandomWeights( weights );

  cout<< "Weights" << endl;
  printVector( weights );

  cout << "Error in estimation is = " << leastSquaresError( x, t, weights ) << endl;
*/

