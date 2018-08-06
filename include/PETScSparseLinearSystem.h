/* PETScSparseLinearSystem class - Specification of the sparse linear system
  class. This class links to complex PETSc to perform the solver phase.
*/

#ifndef PETSCSPARSELINEARSYSTEM_H
#define PETSCSPARSELINEARSYSTEM_H

#include <set>

#include "SparseMatrix.h"
#include "Matrix.h"
#include "Vector.h"

#if defined(PETSC_Z) || defined(PETSC_D)
   #include "petscksp.h"
#endif

namespace TSL
{
  template <typename T>
  class PETScSparseLinearSystem
  {

  public:

    /// Constructor for a PETSc linear system object.
    PETScSparseLinearSystem( SparseMatrix<T>* Aptr, Vector<T>* Bptr );

    /// Destructor for a linear system object.
    ~PETScSparseLinearSystem();

    /// deallocates some objects
    void cleanup();

    /// Solve the sparse system
    void solve();

    /// Factorise the Ax=B system
    void factorise();

    /// Resolve the same system using the same factorisation
    void solve_using_factorisation();


  private:
    /// Solve by using PETSc
    void solve_petsc();

    /// Factorise by linking to PETSc
    void factorise_petsc();  

    /// Pointer to a sparse LHS matrix
    SparseMatrix<T>* p_A;
    /// Pointer to the RHS vector
    Vector<T>* p_B;

    /// indicates that the matrix has been factorised
    bool factorised_;

    #if defined(PETSC_Z) || defined(PETSC_D)
      Vec            x_,B_;      /* B = RHS and x = soln */
      Mat            F_;
      KSP            ksp_;       /* linear solver context */
      PC             pc_;        /* preconditioner -- though hard wired for MUMPS direct method */
      PetscMPIInt    rank_, size_;
    #endif

  };

} //end namepsace
#endif
