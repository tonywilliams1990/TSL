/* Eigensystem class - Solves generalised eigenvalue problems
        of the form A*v=lambda*B*v where A and B are matrices,
        lambda is an eigenvalue and v is an eigenvector
*/

#include "Eigensystem.h"
#include "LAPACK.h"

namespace TSL
{

	/* ----- Methods ----- */

  /// Compute the eigenvalues ( and optionally the eigenvectors ) for real matrices
  template <>
  void Eigensystem<double>::compute( const Matrix<double>& A, const Matrix<double>& B,
                                bool compute_evecs )
  {
    // Check matrix dimensions agree
    if ( A.ROWS != B.ROWS || A.COLS != B.COLS )
    {
      std::string problem = "Eigensystem: Matrix dimensions do not agree.";
      throw Error( problem );
    }
    // Check matrix is square
    if ( A.ROWS != A.COLS )
    {
      std::string problem = "Eigensystem: Matrix must be square.";
      throw Error( problem );
    }
    N = A.ROWS;
    // Create dynamically sized generalized Eigen solver
    /*Eigen::GeneralizedEigenSolver< Eigen::Matrix<T,-1,-1> > ges;*/
		// Try memory preallocation construction
		Eigen::GeneralizedEigenSolver< Eigen::Matrix<double, -1, -1> > ges( N );
    ges.setMaxIterations( 100000 );
    // Compute the eigenvalues
    if ( compute_evecs )
    {
      ges.compute( A.MATRIX, B.MATRIX, true );
    }
    else
    {
      ges.compute( A.MATRIX, B.MATRIX );
    }
    EIGENVALUES_COMPUTED = true;

    // Put the eigenvalues into the storage vector
    for ( std::size_t i=0; i<N; ++i )
    {
      ALPHAS.push_back( ges.alphas()[i] );
      BETAS.push_back( ges.betas()[i] );
      // If beta[i] = 0 the eigenvalue is at infinity so set to very large number
      if ( ges.betas()[i] == 0 )
      {
        EIGENVALUES.push_back( ges.alphas()[i] / 1e-100 );
      }
      else
      {
        EIGENVALUES.push_back( ges.eigenvalues()[i] );
      }
    }

    // Store eigenvectors if required
    if ( compute_evecs )
    {
      EIGENVECTORS_COMPUTED = true;

      // Put the eigenvectors into the storage matrix
      EIGENVECTORS.resize( A.ROWS, A.COLS );
      for ( std::size_t i=0; i<A.ROWS; ++i )
      {
        for ( std::size_t j=0; j<A.COLS; ++j )
        {
          EIGENVECTORS( i, j ) = ges.eigenvectors()( i, j );
        }
      }
    }

  }

  /// Compute the eigenvalues ( and optionally the eigenvectors ) for complex matrices
  template <>
  void Eigensystem<std::complex<double>>::compute( const Matrix<std::complex<double>>& A,
                                                   const Matrix<std::complex<double>>& B,
                                                   bool compute_evecs )
  {
    // Check matrix dimensions agree
    if ( A.ROWS != B.ROWS || A.COLS != B.COLS )
    {
      std::string problem = "Eigensystem: Matrix dimensions do not agree.";
      throw Error( problem );
    }
    // Check matrix is square
    if ( A.ROWS != A.COLS )
    {
      std::string problem = "Eigensystem: Matrix must be square.";
      throw Error( problem );
    }
    N = A.ROWS;

    // Setup LAPACK routine zggev_
    char jobvl = 'N';             // Compute the left eigenvectors N = no, V = yes
    char jobvr = 'N';             // Compute the right eigenvectors
    if ( compute_evecs ) { jobvr = 'V'; }
    int n = static_cast<int>(N);
    int lda = n;
    int ldb = n;
    int ldvl = n;
    int ldvr = n;
    int info, lwork;
    std::complex<double> wkopt;
    std::complex<double>* work;
    double rwork[2*n];
    std::complex<double> vl[ldvl*n], vr[ldvr*n];

    Vector< std::complex<double> > alphas( n, 0.0 );
    Vector< std::complex<double> > betas( n, 0.0 );

    // Convert A and B matrices to vectors
    std::vector< std::complex<double> > a_vec( n * n, 0.0);
    std::vector< std::complex<double> > b_vec( n * n, 0.0);
    for ( int i=0; i<n; ++i)
    {
      for ( int j=0; j<n; ++j)
      {
        a_vec[ i * n + j ] = A(i,j);
        b_vec[ i * n + j ] = B(i,j);
      }
    }

    // Query and allocate the optimal workspace
    lwork = -1;
    zggev_( &jobvl, &jobvr, &n, &a_vec[0], &lda, &b_vec[0], &ldb, &alphas[0], &betas[0], vl, &ldvl, vr, &ldvr,
            &wkopt, &lwork, rwork, &info );
    lwork = (int)wkopt.real();
    work = (std::complex<double>*)malloc( lwork*sizeof(std::complex<double>) );
    // Solve eigenproblem
    zggev_( &jobvl, &jobvr, &n, &a_vec[0], &lda, &b_vec[0], &ldb, &alphas[0], &betas[0], vl, &ldvl, vr, &ldvr,
    work, &lwork, rwork, &info );

    if( info > 0 ) {
      std::string problem = "zggev_ failed to compute eigenvalues.";
      throw Error( problem );
    }

    EIGENVALUES_COMPUTED = true;
    ALPHAS = alphas;
    BETAS = betas.real();
    // Put the eigenvalues into the storage vector
    for ( int i=0; i<n; ++i )
    {
      // If betas[i] = 0 the eigenvalue is at infinity so set to very large number
      if ( betas[i] == 0.0 )
      {
        EIGENVALUES.push_back( alphas[i] / 1e-100 );
      }
      else
      {
        EIGENVALUES.push_back( alphas[i] / betas[i] );
      }
    }

    // Eigenvectors
    if ( compute_evecs )
    {
      EIGENVECTORS_COMPUTED = true;

      // Put the eigenvectors into the storage matrix
      EIGENVECTORS.resize( n, n );
      for ( std::size_t i=0; i<N; ++i )
      {
        for ( std::size_t j=0; j<N; ++j )
        {
          EIGENVECTORS( i, j ) = vr[ i * n + j ];
        }
      }
    }

  }


  // Templated versions
  template class Eigensystem<double>;
	template class Eigensystem< std::complex<double> >;

} // End of namespace TSL
