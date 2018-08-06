/* SparseEigenSystem class - Specification of the sparse linear eigensystem
  class. This class links to complex SLEPc to perform the solver phase.
*/

#ifndef SPARSEEIGENSYSTEM_H
#define SPARSEEIGENSYSTEM_H

#include <set>

#include "SparseMatrix.h"
#include "Matrix.h"
#include "Vector.h"

#ifdef SLEPC
#include <slepceps.h>

#include <slepc.h>
#include <petscksp.h>

namespace TSL
{
  /// A linear Nth-order generalised eigensystem class.
  template <typename T>
  class SparseEigenSystem
  {

  public:

    /// Constructor for a linear system object.
    SparseEigenSystem( SparseMatrix<T>* Aptr, SparseMatrix<T>* Bptr );

    /// Destructor for a linear system object.
    ~SparseEigenSystem();

    /// Access the (actual) number of (converged) eigenvalues found.
    unsigned get_nconv() const;

    /// Request a certain number of eigenvalues
    void set_nev( unsigned n );

    /// Set target for the shift-invert algorithm.
    /// \param target The target eigenvalue
    void set_target( std::complex<double> target );

    /// Defines the ordering of returned set of eigenvalues/vectors
    /// \param order_string A string that defines the SLEPc enum options for ordering (see EPSWhich)
    void set_order( std::string order_string );

    /// Gives a handle to the boolean that sets a region in the complex plane
    bool& region_defined();

    /// Return a handle to the CALC_EIGENVECTORS variables
    bool& calc_eigenvectors();

    /// Set a rectangular region of the complex plane in which to look for eigenvalues.
    void set_region( const double& left, const double& right, const double& bottom, const double& top );

    /// Gives a handle to the boolean that sets if an initial guess has been used
    bool& guess_defined();

    /// Set the initial guess
    /// \param guess The vector of the guess
    void set_initial_guess( const Vector<T>& guess );

    /// Solve the matrix linear eigensystem
    void eigensolve();

    /// Return all the eigenvalues
    Vector< std::complex<double> > eigenvalues() const
    {
      return ALL_EIGENVALUES;
    }

    /// Return all the eigenvectors in a matrix
    Matrix< std::complex<double> > eigenvectors() const
    {
      return ALL_EIGENVECTORS;
    }

    /// Return the tagged eigenvalues TODO
    Vector< std::complex<double> > get_tagged_eigenvalues() const;

    /// Return the tagged eigenvectors TODO
    Matrix< std::complex<double> > get_tagged_eigenvectors();

  private:

    //static PetscErrorCode monitor(EPS eps,int its,int nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal* errest,int nest,void *mctx);

    /// Solve the generalised eigenproblem and compute eigenvectors
    void eigensolve_slepc();

    /// number of eigenvalues requested
    unsigned NEV;
    /// number of (converged) eigenvalues located by the solver -- will be changed by the solver.
    unsigned NCONV;

    /// defines how to seek and order eigenvalues (e.g. smallest real part?)
    EPSWhich ORDER;

    /// defines if ev's are to be sought in a specific region
    bool REGION_DEFINED;
    /// defines if an initial guess is to be used
    bool GUESS_DEFINED;

    /// four numbers that define a rectangular region in the complex plane
    double REAL_L,REAL_R,IMAG_B,IMAG_T;

    /// stores an initial guess to work with
    Vector<T> INITIAL_GUESS;

    /// pointer to the LHS matrix
    SparseMatrix<T>* p_A;
    /// pointer to the RHS matrix
    SparseMatrix<T>* p_B;

    /// storage for eigenvectors and eigenvalues
    Vector< std::complex<double> > ALL_EIGENVALUES;
    Matrix< std::complex<double> > ALL_EIGENVECTORS;

    bool CALC_EIGENVECTORS;
    std::complex<double> SHIFT;

  };

} //end namepsace
#endif

#endif
