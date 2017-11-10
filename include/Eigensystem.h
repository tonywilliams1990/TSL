/* Eigensystem class - Solves generalised eigenvalue problems
        of the form A*v=lambda*B*v where A and B are matrices,
        lambda is an eigenvalue and v is an eigenvector
*/

#ifndef EIGENSYSTEM_H
#define EIGENSYSTEM_H

#include <complex>
#include <Eigen/Eigenvalues>

#include "Vector.h"
#include "Matrix.h"
#include "Error.h"

namespace TSL
{

  // TODO maybe the class doesn't need to be templated as it will always be complex<double>?

  template <class T>
  class Eigensystem
  {
    protected:
      std::size_t N;                                  // Number of eigenvalues
      Vector< std::complex<double> > EIGENVALUES;     // Vector for eigenvalues
      Vector< std::complex<double> > ALPHAS;          // Vector for complex numerators
      Vector< double > BETAS;                         // Vector for real denominators
      Matrix< std::complex<double> > EIGENVECTORS;    // Matrix for eigenvectors
      bool EIGENVALUES_COMPUTED;                      // True if eigenvalues are calculated
      bool EIGENVECTORS_COMPUTED;                     // True if eigenvectors are calculated

    public:

      /* ----- Constructors and Destructor ----- */

      /// Constructor
      Eigensystem() : EIGENVALUES_COMPUTED( false ), EIGENVECTORS_COMPUTED( false ) { }

      /// Destructor
      ~Eigensystem() { }

      /* ----- Methods ----- */

      /// Return true if the eigenvalues have been computed
      bool eigenvalues_computed()  { return EIGENVALUES_COMPUTED; }

      /// Return true if the eigenvectors have been computed
      bool eigenvectors_computed()  { return EIGENVECTORS_COMPUTED; }

      /// Compute the eigenvalues ( and optionally the eigenvectors ) for real matrices
      void compute( const Matrix<double>& A, const Matrix<double>& B, bool compute_evecs = true );

      /// Compute the eigenvalues ( and optionally the eigenvectors ) for complex matrices
      void compute( const Matrix<std::complex<double>>& A, const Matrix<std::complex<double>>& B, bool compute_evecs = true );

      //TODO non-generalised A*x = lambda*x (no B matrix) -> zgeev

      /// Return the computed eigenvalues in a vector
      Vector< std::complex<double> > eigenvalues() const
      {
        if ( EIGENVALUES_COMPUTED ) { return EIGENVALUES; }
        else { throw Error( "Eigensystem: eigenvalues not computed." ); }
      }

      /// Return the complex numerators of the eigenvalues
      Vector< std::complex<double> > alphas() const
      {
        if ( EIGENVALUES_COMPUTED ) { return ALPHAS; }
        else { throw Error( "Eigensystem: eigenvalues not computed." ); }
      }

      /// Return the real denominators of the eigenvalues
      Vector<double> betas() const
      {
        if ( EIGENVALUES_COMPUTED ) { return BETAS; }
        else { throw Error( "Eigensystem: eigenvalues not computed." ); }
      }

      /// Return the matrix of eigenvectors ( each column is an eigenvector )
      Matrix< std::complex<double> > eigenvector_matrix() const
      {
        if ( EIGENVECTORS_COMPUTED ) { return EIGENVECTORS; }
        else { throw Error( "Eigensystem: eigenvectors not computed." ); }
      }

      /// Return an std::vector of eigenvectors
      std::vector< Vector< std::complex<double> > > eigenvectors() const
      {
        if ( EIGENVECTORS_COMPUTED )
        {
          std::vector< Vector< std::complex<double> > > evecs;
          std::size_t rows = EIGENVECTORS.rows();
          std::size_t cols = EIGENVECTORS.cols();

          for (std::size_t j=0; j<cols; ++j)
          {
              Vector< std::complex<double> > evec;
              for (std::size_t i=0; i<rows; ++i)
              {
                  evec.push_back( EIGENVECTORS( i, j ) );
              }
              evecs.push_back( evec );
          }
          return evecs;
        }
        else { throw Error( "Eigensystem: eigenvectors not computed." ); }
      }

  }; // End of class Eigensystem

} // End of namespace TSL

#endif
