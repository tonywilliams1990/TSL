// LAPACK test code

#include <iostream>
#include <vector>

#include "LAPACK.h"
#include "Core"

using namespace std;
using namespace TSL;

//extern "C" { void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv, int* info); }
//extern "C" { void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO ); }

int main()
{
    /*char trans = 'N';
    int dim = 2;
    int nrhs = 1;
    int LDA = dim;
    int LDB = dim;
    int info;

    vector<double> a, b;

    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(-1);

    b.push_back(2);
    b.push_back(0);

    int ipiv[3];

    TSL::dgetrf_(&dim, &dim, &*a.begin(), &LDA, ipiv, &info);
    TSL::dgetrs_(&trans, &dim, &nrhs, & *a.begin(), &LDA, ipiv, & *b.begin(), &LDB, &info);


    std::cout << "solution is:";
    std::cout << "[" << b[0] << ", " << b[1] << ", " << "]" << std::endl;
    std::cout << "Info = " << info << std::endl;*/

   /*
   ZGEEV Example.
   ==============

   Program computes the eigenvalues and left and right eigenvectors of a general
   rectangular matrix A:

   ( -3.84,  2.25) ( -8.94, -4.75) (  8.95, -6.53) ( -9.87,  4.82)
   ( -0.66,  0.83) ( -4.40, -3.82) ( -3.50, -4.26) ( -3.15,  7.36)
   ( -3.99, -4.73) ( -5.88, -6.60) ( -3.36, -0.40) ( -0.75,  5.23)
   (  7.74,  4.18) (  3.66, -7.53) (  2.58,  3.60) (  4.59,  5.41)

   Eigenvalues
   ( -9.43,-12.98) ( -3.44, 12.69) (  0.11, -3.40) (  5.76,  7.13)

   Left eigenvectors
   (  0.24, -0.18) (  0.61,  0.00) ( -0.18, -0.33) (  0.28,  0.09)
   (  0.79,  0.00) ( -0.05, -0.27) (  0.82,  0.00) ( -0.55,  0.16)
   (  0.22, -0.27) ( -0.21,  0.53) ( -0.37,  0.15) (  0.45,  0.09)
   ( -0.02,  0.41) (  0.40, -0.24) (  0.06,  0.12) (  0.62,  0.00)

   Right eigenvectors
   (  0.43,  0.33) (  0.83,  0.00) (  0.60,  0.00) ( -0.31,  0.03)
   (  0.51, -0.03) (  0.08, -0.25) ( -0.40, -0.20) (  0.04,  0.34)
   (  0.62,  0.00) ( -0.25,  0.28) ( -0.09, -0.48) (  0.36,  0.06)
   ( -0.23,  0.11) ( -0.10, -0.32) ( -0.43,  0.13) (  0.81,  0.00)
   */

    // Test zgeev
    /*int n = 4;
    int lda = 4;
    int ldvl = 4;
    int ldvr = 4;
    int info, lwork;

    std::complex<double> wkopt;
    std::complex<double>* work;

    double rwork[2*n];
    //std::complex<double> w[n];
    std::complex<double> vl[ldvl*n], vr[ldvr*n];
    std::complex<double> a[lda*n] = {
      {-3.84,  2.25}, {-0.66,  0.83}, {-3.99, -4.73}, { 7.74,  4.18},
      {-8.94, -4.75}, {-4.40, -3.82}, {-5.88, -6.60}, { 3.66, -7.53},
      { 8.95, -6.53}, {-3.50, -4.26}, {-3.36, -0.40}, { 2.58,  3.60},
      {-9.87,  4.82}, {-3.15,  7.36}, {-0.75,  5.23}, { 4.59,  5.41}
    };

    // Executable statements
    cout << endl << " ZGEEV Example Program Results" << endl;
    // Query and allocate the optimal workspace
    lwork = -1;
    zgeev_( ( char* ) "V", ( char* ) "V", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
            &wkopt, &lwork, rwork, &info );
    lwork = (int)wkopt.real();
    work = (std::complex<double>*)malloc( lwork*sizeof(std::complex<double>) );
    // Solve eigenproblem
    zgeev_( ( char* ) "V", ( char* ) "V", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
    work, &lwork, rwork, &info );
    // Check for convergence
    if( info > 0 ) {
          cout << "The algorithm failed to compute eigenvalues." << endl;
    }

    cout << "Eigenvalues = ";
    for ( int i=0; i<n; ++i)
    {
      cout << w[i] << " ";
    }
    cout << endl << "Eigenvectors (right) = " << endl;
    for ( int i=0; i<n; ++i)
    {
      for ( int j=0; j<n; ++j)
      {
        cout << vr[ j * ldvr + i] << " ";
      }
      cout << endl;
    }*/

    // Test zggev
    char jobvl = 'N';             // Compute the left eigenvectors N = no, V = yes
    char jobvr = 'V';             // Compute the right eigenvectors
    int n = 4;
    int lda = 4;
    int ldb = 4;
    int ldvl = 4;
    int ldvr = 4;
    int info, lwork;
    std::complex<double> wkopt;
    std::complex<double>* work;
    double rwork[2*n];
    std::complex<double> vl[ldvl*n], vr[ldvr*n];

    TSL::Vector< std::complex<double> > alphas( n, 0.0 );
    TSL::Vector< std::complex<double> > betas( n, 0.0 );

    // A matrix
    Matrix< std::complex<double> > A( 4, 4, 0.0 );

    A( 0, 0 ) = std::complex<double>( -21.10, -22.50 );
    A( 0, 1 ) = std::complex<double>( 53.50, -50.50 );
    A( 0, 2 ) = std::complex<double>( -34.50, 127.50 );
    A( 0, 3 ) = std::complex<double>( 7.50, 0.50 );

    A( 1, 0 ) = std::complex<double>( -0.46, -7.78 );
    A( 1, 1 ) = std::complex<double>( -3.50, -37.50 );
    A( 1, 2 ) = std::complex<double>( -15.50, 58.50 );
    A( 1, 3 ) = std::complex<double>( -10.50, -1.50 );

    A( 2, 0 ) = std::complex<double>( 4.30, -5.50 );
    A( 2, 1 ) = std::complex<double>( 39.70, -17.10 );
    A( 2, 2 ) = std::complex<double>( -68.50, 12.50 );
    A( 2, 3 ) = std::complex<double>( -7.50, -3.50 );

    A( 3, 0 ) = std::complex<double>( 5.50, 4.40 );
    A( 3, 1 ) = std::complex<double>( 14.40, 43.30 );
    A( 3, 2 ) = std::complex<double>( -32.50, -46.00 );
    A( 3, 3 ) = std::complex<double>( -19.00, -32.50 );

    cout << "A = " << A << endl;

    // Convert to vector
    std::vector< std::complex<double> > a_vec( n * n, 0.0);
    for ( int i=0; i<n; ++i)
    {
      for ( int j=0; j<n; ++j)
      {
        a_vec[ i * n + j ] = A(i,j);
      }
    }

    /*std::vector< std::complex<double> > a_vec = {
      {-21.10, -22.50}, { 53.50, -50.50}, {-34.50, 127.50}, { 7.50,  0.50  },
      { -0.46, -7.78 }, { -3.50, -37.50}, {-15.50, 58.50 }, {-10.50, -1.50 },
      {  4.30, -5.50 }, { 39.70, -17.10}, {-68.50, 12.50 }, { -7.50, -3.50 },
      {  5.50, 4.40  }, { 14.40, 43.30 }, {-32.50, -46.00}, {-19.00, -32.50}
    };*/

    // B matrix
    Matrix< std::complex<double> > B( 4, 4, 0.0 );

    B( 0, 0 ) = std::complex<double>( 1.00, -5.00 );
    B( 0, 1 ) = std::complex<double>( 1.60, 1.20 );
    B( 0, 2 ) = std::complex<double>( -3.00, 0.00 );
    B( 0, 3 ) = std::complex<double>( 0.00, -1.00 );

    B( 1, 0 ) = std::complex<double>( 0.80, -0.60 );
    B( 1, 1 ) = std::complex<double>( 3.00, -5.00 );
    B( 1, 2 ) = std::complex<double>( -4.00, 3.00 );
    B( 1, 3 ) = std::complex<double>( -2.40, -3.20 );

    B( 2, 0 ) = std::complex<double>( 1.00, 0.00 );
    B( 2, 1 ) = std::complex<double>( 2.40, 1.80 );
    B( 2, 2 ) = std::complex<double>( -4.00, -5.00 );
    B( 2, 3 ) = std::complex<double>( 0.00, -3.00 );

    B( 3, 0 ) = std::complex<double>( 0.00, 1.00 );
    B( 3, 1 ) = std::complex<double>( -1.80, 2.40 );
    B( 3, 2 ) = std::complex<double>( 0.00, -4.00 );
    B( 3, 3 ) = std::complex<double>( 4.00, -5.00 );

    cout << "B = " << B << endl;
    // Convert to vector
    std::vector< std::complex<double> > b_vec( n * n, 0.0);
    for ( int i=0; i<n; ++i)
    {
      for ( int j=0; j<n; ++j)
      {
        b_vec[ i * n + j ] = B(i,j);
      }
    }

    /*std::vector< std::complex<double> > b_vec = {
      { 1.00, -5.00}, { 1.60,  1.20}, {-3.00,  0.00}, { 0.00, -1.00},
      { 0.80, -0.60}, { 3.00, -5.00}, {-4.00,  3.00}, {-2.40, -3.20},
      { 1.00,  0.00}, { 2.40,  1.80}, {-4.00, -5.00}, { 0.00, -3.00},
      { 0.00,  1.00}, {-1.80,  2.40}, { 0.00, -4.00}, { 4.00, -5.00}
    };*/
    /* Executable statements */
    cout << endl << " ZGGEV Example Program Results" << endl;
    /* Query and allocate the optimal workspace */
    lwork = -1;
    zggev_( &jobvl, &jobvr, &n, &a_vec[0], &lda, &b_vec[0], &ldb, &alphas[0], &betas[0], vl, &ldvl, vr, &ldvr,
            &wkopt, &lwork, rwork, &info );
    lwork = (int)wkopt.real();
    work = (std::complex<double>*)malloc( lwork*sizeof(std::complex<double>) );
    /* Solve eigenproblem */
    zggev_( &jobvl, &jobvr, &n, &a_vec[0], &lda, &b_vec[0], &ldb, &alphas[0], &betas[0], vl, &ldvl, vr, &ldvr,
    work, &lwork, rwork, &info );

    if( info > 0 ) {
          cout << "The algorithm failed to compute eigenvalues." << endl;
    }
    cout << "Alphas = ";
    for ( int i=0; i<n; ++i)
    {
      cout << alphas[i] << " ";
    }
    cout << endl << "Betas = ";
    for ( int i=0; i<n; ++i)
    {
      cout << betas[i] << " ";
    }
    cout << endl << "Eigenvalues = ";
    for ( int i=0; i<n; ++i)
    {
      cout << alphas[i] / betas[i] << " ";
    }
    cout << endl << "Eigenvectors (right) = " << endl;
    for ( int i=0; i<n; ++i)
    {
      for ( int j=0; j<n; ++j)
      {
        cout << vr[ j * ldvr + i] << " ";
      }
      cout << endl;
    }

    /* Free workspace */
    free( (void*)work );

    cout << endl << "FINISHED" << endl;

    return(0);
}
