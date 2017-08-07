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

void printMatrix( const double *x, const int nRows, const int nCols){
      for (int i=0; i < nRows; i++) {
	    for (int j=0; j < nCols; j++) {
		  printf ("%12.3f", x[i*nCols +j]);
	    }
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
      const int ORDER = 2;
      
      x = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      t = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      w = (double *)mkl_malloc( (ORDER)*sizeof( double ), 64 );

      memset( x, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( t, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( w, 0.0,  (ORDER) * sizeof(double));      

      loadData( x , inputsFile );
      loadData( t , targetsFile );
      
      cout << "Inputs" << endl;
      printVector( x , NUM_PATTERNS );

      cout << "Targets" << endl;
      printVector( t , NUM_PATTERNS );


      cout << "\nComputing Design Matrix ..." << endl;
      double *designMatrix = (double *)mkl_malloc( NUM_PATTERNS*ORDER*sizeof( double ), 64 );
      memset( designMatrix, 0.0, NUM_PATTERNS * ORDER* sizeof(double));

      //set first column to 1--dummy index to calculate w0
      // maybe in an opportunity for blas routines?
      for (int i = 0; i < NUM_PATTERNS*ORDER; ++i) {
	    designMatrix[i] = 1.0;
	    i++;
      }

      // phi(x) = x--> identity basis function for linear regression
      for (int i = 0; i < NUM_PATTERNS; ++i) {
	    for (int j = 0; j < ORDER; ++j) {
		  if (( j % 2) != 0) {
		  designMatrix[i*ORDER + j] = x[i];			
		  }
	    }
      }

      printf ("\n Design Matrix: \n");
      printMatrix( designMatrix, NUM_PATTERNS, ORDER );

      //printf (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
      double *designTranspose = (double *)mkl_malloc( ORDER*NUM_PATTERNS*sizeof( double ), 64 );
      memset( designTranspose, 0.0, ORDER* NUM_PATTERNS*sizeof(double));

      //perform transpose
      double alpha, beta;
      alpha = 1.0;
      beta = 0.0;

      mkl_domatcopy('R' , 'T' ,  NUM_PATTERNS, ORDER, alpha, designMatrix, ORDER ,designTranspose, NUM_PATTERNS);

      printf ("\n Design Matrix Transpose: \n");
      printMatrix( designTranspose, ORDER, NUM_PATTERNS );

      //TODO
      //2. Inversion
      //3. multiplication

      // A = (Phi'*Phi)
      // B = Phi'*t
      // w = A^-1 * B
      double *A = (double *)mkl_malloc( ORDER*ORDER*sizeof( double ), 64 );
      double *B = (double *)mkl_malloc( ORDER*sizeof( double ), 64 );

      memset( A, 0.0, ORDER* ORDER*sizeof(double));
      memset( B, 0.0, ORDER*sizeof(double));
      
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  ORDER, ORDER, NUM_PATTERNS, alpha, designTranspose, NUM_PATTERNS, designMatrix, ORDER, beta, A, ORDER);

      printf ("\n Matrix Product A : \n");
      printMatrix( A, ORDER, ORDER );

      //calculate inverse
      int *IPIV = new int[ORDER + 1];
      int LWORK = ORDER*ORDER;
      double *WORK = new double[LWORK];
      int INFO;

      dgetrf( &ORDER, &ORDER, A, &ORDER, IPIV, &INFO );
      dgetri( &ORDER, A, &ORDER, IPIV, WORK, &LWORK, &INFO );

      delete IPIV;
      delete WORK;

      printf ("\n Inverse A : \n");
      printMatrix( A, ORDER, ORDER );

      //get B
      cblas_dgemv( CblasRowMajor, CblasNoTrans, ORDER, NUM_PATTERNS,
		   alpha, designTranspose, NUM_PATTERNS, t, 1, beta, B, 1);

      printf(" Vector B\n");
      printVector( B , ORDER );      

      //final solution should just be A * b
      cblas_dgemv( CblasRowMajor, CblasNoTrans, ORDER, ORDER,
		   alpha, A, ORDER, B, 1, beta, w, 1);
      
      cout << "\nComputing Line of Best Fit ..." << endl;
      //w[1] = computeSlope( x , t);
      //w[0] = computeIntercept( x , t);
      cout << " w0 = " << w[0] << endl;
      cout << " w1 = " << w[1] << endl;
      
      printf ("\n Deallocating memory \n\n");
      mkl_free( x );
      mkl_free( t );
      mkl_free( w );
      mkl_free( designMatrix );
      mkl_free( designTranspose );
      
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

