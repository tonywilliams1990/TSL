/// Function definitions for LAPACK routines
#include <complex>

namespace TSL
{
  /// LU factorization of a general M-by-N matrix A
  extern "C" {
    void dgetrf_(int* M, int* N, double* A, int* lda, int* ipiv, int* info);
   }

  /// Solves a system of linear equations A * X = B or  A**T * X = B
  extern "C" {
    void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );
  }

  /// Solve the eigenvalue problem A * x = lambda * x for a complex matrix A
  extern "C" {
    void zgeev_( char* jobvl, char* jobvr, int* n, std::complex<double>* a,
                int* lda, std::complex<double>* w, std::complex<double>* vl,
                int* ldvl, std::complex<double>* vr, int* ldvr,
                std::complex<double>* work, int* lwork, double* rwork, int* info );
  }

  /// Solve the generalised eigenvalue problem A * x = lambda * B * x for complex matrices A & B
  extern "C"{
    void zggev_(char *jobvl, char *jobvr, int *n, std::complex<double> *a,
                int *lda, std::complex<double> *b, int *ldb, std::complex<double> *alpha,
                std::complex<double> *beta, std::complex<double> *vl,
                int *ldvl, std::complex<double> *vr, int *ldvr,
                std::complex<double> *work, int *lwork, double *rwork, int *info);
  }

} // End of namespace TSL
