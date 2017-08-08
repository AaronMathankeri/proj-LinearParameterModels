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

const int NUM_PATTERNS = 10;
const int ORDER = 2;

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

void computeDesignMatrix( double *x, double *Phi ){
      //set first column to 1--dummy index to calculate w0
      // maybe in an opportunity for blas routines?
      for (int i = 0; i < NUM_PATTERNS*ORDER; ++i) {
	    Phi[i] = 1.0;
	    i++;
      }
      // phi(x) = x--> identity basis function for linear regression
      for (int i = 0; i < NUM_PATTERNS; ++i) {
	    for (int j = 0; j < ORDER; ++j) {
		  if (( j % 2) != 0) {
			Phi[i*ORDER + j] = x[i];			
		  }
	    }
      }
}

void solveNormalEquations( double *moorePenrosePhi, double *t, double *w ){
      // compute normal equations...
      //final solution should just be Moore-Penrose * t
      double alpha = 1.0;
      double beta = 0.0;
      cblas_dgemv( CblasRowMajor, CblasNoTrans, ORDER, NUM_PATTERNS,
      		   alpha, moorePenrosePhi, NUM_PATTERNS, t, 1, beta, w, 1);
}

void computeMoorePenroseInverse( double* Phi , double *MoorePenrose ){
      double *PhiTranspose = (double *)mkl_malloc( ORDER*NUM_PATTERNS*sizeof( double ), 64 );
      double *A = (double *)mkl_malloc( ORDER*ORDER*sizeof( double ), 64 );
      double alpha = 1.0;
      double beta = 0.0;
      //declare MKL variables for inverse calculation
      int LWORK = ORDER*ORDER;
      int INFO;
      int *IPIV = (int *)mkl_malloc( (ORDER+1)*sizeof( int ), 64 );
      double *WORK = (double *)mkl_malloc( LWORK*sizeof( double ), 64 );

      memset( PhiTranspose, 0.0, ORDER* NUM_PATTERNS*sizeof(double));
      memset( A, 0.0, ORDER* ORDER*sizeof(double));
      //memset( B, 0.0, ORDER*NUM_PATTERNS*sizeof(double));

      // calculate transpose
      mkl_domatcopy('R' , 'T' ,  NUM_PATTERNS, ORDER, alpha,
		    Phi, ORDER, PhiTranspose, NUM_PATTERNS);

      //calculate product
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  ORDER, ORDER, NUM_PATTERNS, alpha, PhiTranspose,
		  NUM_PATTERNS, Phi, ORDER, beta, A, ORDER);
      //calculate inverse
            dgetrf( &ORDER, &ORDER, A, &ORDER, IPIV, &INFO );
      dgetri( &ORDER, A, &ORDER, IPIV, WORK, &LWORK, &INFO );

      // A * phi' = moore penrose!
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  ORDER, NUM_PATTERNS, ORDER, alpha, A,
		  ORDER, PhiTranspose, NUM_PATTERNS, beta, MoorePenrose, NUM_PATTERNS);

      mkl_free( PhiTranspose );
      mkl_free( A );
      mkl_free( IPIV );
      mkl_free( WORK );
}

int main(int argc, char *argv[])
{
      cout << " Aaron's Back." << endl;

      string inputsFile = "./data/linearRegression/inputs.txt";
      string targetsFile = "./data/linearRegression/targets.txt";

      // declare variables for calculations
      double *x, *t, *w, *designMatrix, *moorePenroseDesignInverse;

      x = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      t = (double *)mkl_malloc( NUM_PATTERNS*sizeof( double ), 64 );
      w = (double *)mkl_malloc( (ORDER)*sizeof( double ), 64 );
      designMatrix = (double *)mkl_malloc( NUM_PATTERNS*ORDER*sizeof( double ), 64 );
      moorePenroseDesignInverse = (double *)mkl_malloc( ORDER*NUM_PATTERNS*sizeof( double ), 64 );
      
      memset( x, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( t, 0.0,  NUM_PATTERNS * sizeof(double));
      memset( w, 0.0,  ORDER* sizeof(double));
      memset( designMatrix, 0.0, NUM_PATTERNS * ORDER* sizeof(double));     
      memset( moorePenroseDesignInverse, 0.0,  ORDER*NUM_PATTERNS * sizeof(double));      

      loadData( x , inputsFile );
      loadData( t , targetsFile );
      
      cout << "Inputs" << endl;
      printVector( x , NUM_PATTERNS );

      cout << "Targets" << endl;
      printVector( t , NUM_PATTERNS );
      //--------------------------------------------------------------------------------
      printf ("\n Design Matrix: \n");
      computeDesignMatrix( x, designMatrix );
      printMatrix( designMatrix, NUM_PATTERNS, ORDER );
      //--------------------------------------------------------------------------------
      // compute moore-penrose psuedo inverse of design matrix
      computeMoorePenroseInverse( designMatrix , moorePenroseDesignInverse );
      //--------------------------------------------------------------------------------
      cout << "\nSolving Normal Equations ..." << endl;
      solveNormalEquations( moorePenroseDesignInverse, t , w);

      cout << " w0 = " << w[0] << endl;
      cout << " w1 = " << w[1] << endl;
      //--------------------------------------------------------------------------------
      printf ("\n Deallocating memory \n\n");
      mkl_free( x );
      mkl_free( t );
      mkl_free( w );
      mkl_free( designMatrix );
      mkl_free( moorePenroseDesignInverse );
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

