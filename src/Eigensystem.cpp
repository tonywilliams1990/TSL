/* Eigensystem class - Solves generalised eigenvalue problems
        of the form A*v=lambda*B*v where A and B are matrices,
        lambda is an eigenvalue and v is an eigenvector 
*/

#include "Eigensystem.h"

namespace TSL
{

	/* ----- Methods ----- */
  
  /// Compute the eigenvalues ( and optionally the eigenvectors )
  template <typename T>
  void Eigensystem<T>::compute( const Matrix<T>& A, const Matrix<T>& B, 
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
    // Create dynamically sized generalized Eigen solver
    Eigen::GeneralizedEigenSolver< Eigen::Matrix<T,-1,-1> > ges;
    // Compute the eigenvalues
    ges.compute( A.MATRIX, B.MATRIX, true );
    EIGENVALUES_COMPUTED = true; 

    N = ges.eigenvalues().size();

    // Put the eigenvalues into the storage vector
    for ( std::size_t i=0; i<N; ++i )
    {
      EIGENVALUES.push_back( ges.eigenvalues()[i] );
      ALPHAS.push_back( ges.alphas()[i] );
      BETAS.push_back( ges.betas()[i] );
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


  // Templated versions
  template class Eigensystem<double>;
	//template class Eigensystem< std::complex<double> >;

} // End of namespace TSL
