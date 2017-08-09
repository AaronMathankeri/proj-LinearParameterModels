#include "linearParameterModels.hpp"

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

void computePseudoInverse( double* Phi , double *phiPsuedoInverse ){
      double *A = (double *)mkl_malloc( ORDER*ORDER*sizeof( double ), 64 );
      double alpha = 1.0;
      double beta = 0.0;
      //declare MKL variables for inverse calculation
      int LWORK = ORDER*ORDER;
      int INFO;
      int *IPIV = (int *)mkl_malloc( (ORDER+1)*sizeof( int ), 64 );
      double *WORK = (double *)mkl_malloc( LWORK*sizeof( double ), 64 );

      memset( A, 0.0, ORDER* ORDER*sizeof(double));

      //calculate product = Phi' * Phi = A
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
		  ORDER, ORDER, NUM_PATTERNS, alpha, Phi,
		  ORDER, Phi, ORDER, beta, A, ORDER);

      //calculate inverse A = (A)^-1
      dgetrf( &ORDER, &ORDER, A, &ORDER, IPIV, &INFO );
      dgetri( &ORDER, A, &ORDER, IPIV, WORK, &LWORK, &INFO );

      // A * phi' = moore penrose!
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
		  ORDER, NUM_PATTERNS, ORDER, alpha, A,
		  ORDER, Phi, ORDER, beta, phiPsuedoInverse, NUM_PATTERNS);

      mkl_free( A );
      mkl_free( IPIV );
      mkl_free( WORK );
}

void solveNormalEquations( double *inversePhi, double *t, double *w ){
      // compute normal equations...
      //final solution should just be Moore-Penrose * t
      double alpha = 1.0;
      double beta = 0.0;
      cblas_dgemv( CblasRowMajor, CblasNoTrans, ORDER, NUM_PATTERNS,
      		   alpha, inversePhi, NUM_PATTERNS, t, 1, beta, w, 1);
}

//------------------------------------------------------------------------------------------------------
//older functions
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
