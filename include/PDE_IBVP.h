/// A specification of a class for an \f$ n^{th} \f$-order IBVP of the form
/// \f[ M_1( {\underline f}(y,t), y, t )\cdot {\underline f}_t (y,t)+ M_0( {\underline f}(y,t), y, t ) \cdot {\underline f}_y (y,t) = {\underline R}( {\underline f}(y,t), y, t )\,, \f]
/// subject to \f$ n \f$ conditions defined at \f$ y = y_{left} \f$ and
/// \f$ y_{right} \f$ for some components of \f$ {\underline f}(y) \f$.
/// Here \f$ M_{0,1} \f$ are matrices.
/// The solution at the new time step \f$ t+\Delta t \f$ is
/// \f[ {\underline f}^{new} = {\underline F} + {\underline g} \f]
/// where \f$ {\underline F} \f$ is the current guess at the solution
/// and \f$ {\underline g} \f$ is the linearised correction.
/// The solution at the previous time \f$ t \f$ is
/// \f[ {\underline f}^{old} = {\underline O} \f]
/// A Crank-Nicolson method is employed with the linearised problem at the mid-time point
/// \f$ t + \Delta t /2 \f$ being:
/// \f[ \frac{2}{\Delta t} M_1 \cdot {\underline g } + 2 M_0 \cdot {\underline g}_y - J \cdot {\underline g} + J_2 \cdot \frac{\underline F - \underline O}{\Delta t} \cdot {\underline g} + J_1 \cdot \frac{\underline F_y + \underline O_y}{2} \cdot {\underline g} = 2 {\underline R} - \frac{2}{\Delta t} M_2 \cdot ( {\underline F} - {\underline O} ) -  M_1 \cdot ( {\underline F}_y + {\underline O}_y )\f]
/// Where \f$ M_{0,1}, J, J_{1,2}, R \f$ are evaluated at the mid-time step with arguments \f$ \left ( \frac{\underline F + \underline O}{2}, y, t + \frac{\Delta t}{2} \right ) \f$,
/// with \f$ J_{1,2} \f$ denoting the Jacobian of the matrices \f$ \partial {(M_{0,1})}_{ij} / \partial f_k \f$.
/// This problem is solved by second-order central differencing at the spatial (\f$ y \f$) inter-node mid points.

#ifndef PDE_IBVP_H
#define PDE_IBVP_H

#include <Vector.h>
#include <Matrix.h>
#include <SparseMatrix.h>
#include <Equation_2matrix.h>
#include <Residual_with_coords.h>
#include <OneD_node_mesh.h>

namespace TSL
{

  /// A templated object for real/complex vector system
  /// of unsteady equations.

  template <typename _Type>
  class PDE_IBVP
  {
  public:

    /// The class is defined by a vector function for the system.
    /// \param equation_ptr A pointer to an inherited Equation object.
    /// \param nodes A vector that defines the nodal positions.
    /// \param ptr_to_left_residual A pointer to a residual object that defines the LHS boundary conditions.
    /// \param ptr_to_right_residual A pointer to a residual object that defines the RHS boundary conditions.
    PDE_IBVP( Equation_2matrix<_Type > *equation_ptr,
              const Vector<double>& nodes,
              Residual_with_coords<_Type>* ptr_to_left_residual,
              Residual_with_coords<_Type>* ptr_to_right_residual );

    /// Destructor
    ~PDE_IBVP();

    /// A Crank-Nicolson 'time' stepper.
    void step2( const double& dt );

    /// Assembles the matrix problem for a BVP solve at the
    /// current time level.
    /// \param a The LHS (banded) matrix.
    /// \param b The RHS (dense) vector.
    /// \param dt The 'time step' to be taken.
    void assemble_matrix_problem( SparseMatrix<_Type>& a, Vector<_Type>& b, 
                                  const double& dt );

    /// Return a reference to the current value of the 'timelike/parabolic' coordinate
    /// \return A handle to the current time stored in the object
    double& coord()
    {
      return T;
    }

    /// \return A handle to the solution mesh
    OneD_node_mesh<_Type>& solution();

    /// Access method to the tolerance
    /// \return A handle to the private member data TOL
    double& tolerance()
    {
      return TOL;
    }

    /// Access method to the maximum number of iterations
    /// \return A handle to the private member data MAX_ITERATIONS
    int& max_itns()
    {
      return MAX_ITERATIONS;
    }

  private:
    /// The solution at the current time level
    OneD_node_mesh<_Type> SOLN;
    /// The solution at the previous time step
    OneD_node_mesh<_Type> PREV_SOLN;
    /// tolerance
    double TOL;
    /// The current value of the timelike variable
    double T;
    /// maximum number of iterations
    int MAX_ITERATIONS;
    /// The function associated with this instance.
    Equation_2matrix<_Type > *p_EQUATION;
    /// Pointer to the residual defining the LHS BC
    Residual_with_coords<_Type > *p_LEFT_RESIDUAL;
    /// Pointer to the residual defining the RHS BC
    Residual_with_coords<_Type > *p_RIGHT_RESIDUAL;

  }
  ; // end class

  template <typename _Type>
  inline OneD_node_mesh<_Type>& PDE_IBVP<_Type>::solution()
  {
    return SOLN;
  }


} // end namespace

#endif

