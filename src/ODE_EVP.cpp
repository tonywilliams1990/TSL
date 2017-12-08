/* ODE_EVP - Here we define the ODE_EVP class used for solving ODE linear
             eigenvalue problems.
*/

#include <string>

#include "Equation.h"
#include "ODE_EVP.h"
#include "Error.h"
#include "Eigensystem.h"
#include "Equation_2matrix.h"

namespace TSL
{
  template <typename T>
  ODE_EVP<T>::ODE_EVP( Equation_2matrix<T >* ptr_to_equation,
                           const Vector<double> &nodes,
                           Residual<T>* ptr_to_left_residual,
                           Residual<T>* ptr_to_right_residual ) :
      p_EQUATION( ptr_to_equation ),
      p_LEFT_RESIDUAL( ptr_to_left_residual ),
      p_RIGHT_RESIDUAL( ptr_to_right_residual ),
      NODES( nodes )
  {
    unsigned n( nodes.size() );
    unsigned order( p_EQUATION -> get_order() );
    A_DENSE = Matrix<T>( n * order, n * order, 0.0 );
    B_DENSE = Matrix<T>( n * order, n * order, 0.0 );
  }

  template <typename T>
  ODE_EVP<T>::~ODE_EVP()
  {}

  template <typename T>
  Eigensystem<T>& ODE_EVP<T>::eigensystem()
  {
    return SYSTEM;
  }

  template <typename T>
  void ODE_EVP<T>::eigensolve( bool compute_evecs )
  {
    // Construct the eigenvalue problem
    assemble_dense_problem();
    // solve the system
    if ( compute_evecs )
    {
      SYSTEM.compute( A_DENSE, B_DENSE, true );
      EIGENVECTORS_COMPUTED = true;
    }
    else
    {
      SYSTEM.compute( A_DENSE, B_DENSE, false );
      EIGENVECTORS_COMPUTED = false;
    }
    EIGENVALUES_COMPUTED = true;
  }


  template <typename T>
  void ODE_EVP<T>::assemble_dense_problem()
  {
    // clear the A & B matrices, as they could have been filled with
    // pivoting if this is a second solve.
    A_DENSE.fill( 0.0 );
    B_DENSE.fill( 0.0 );
    // eqn order
    unsigned order( p_EQUATION -> get_order() );
    // Jacobian matrix for the equation
    Matrix<T> jac_midpt( order, order, 0.0 );
    // local state variable and functions
    // the problem has to be linear, so this is a dummy vector.
    Vector<T> temp_dummy( order, 0.0 );
    // number of nodes in the mesh
    std::size_t num_of_nodes( NODES.size() );
    // row counter
    std::size_t row( 0 );

    // update the BC residuals for the current iteration
    p_LEFT_RESIDUAL -> update( temp_dummy );
    // add the (linearised) LHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_LEFT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        A_DENSE( row, var ) = p_LEFT_RESIDUAL -> jacobian()( i, var );
        B_DENSE( row, var ) = 0.0;
      }
      ++row;
    }

    // inner nodes of the mesh, node = 0,1,2,...,num_of_nodes-2
    for ( std::size_t node = 0; node <= num_of_nodes - 2; ++node )
    {
      const std::size_t lnode = node;
      const std::size_t rnode = node + 1;
      // get the current step length
      double h = NODES[ rnode ] - NODES[ lnode ];
      // reciprocal of the spatial step
      const double invh( 1. / h );
      // mid point of the independent variable
      double x_midpt = 0.5 * ( NODES[ lnode ] + NODES[ rnode ] );
      p_EQUATION -> coord(0) = x_midpt;
      // Update the equation to the mid point position
      p_EQUATION -> update( temp_dummy );
      // loop over all the variables and fill the matrix
      for ( unsigned var = 0; var < order; ++var )
      {
        std::size_t placement_row( row + var );
        // deriv at the MID POINT between nodes
        A_DENSE( placement_row, order * rnode + var ) = invh;
        A_DENSE( placement_row, order * lnode + var ) = -invh;
        // add the Jacobian terms to the linearised problem
        for ( unsigned i = 0; i < order; ++i )
        {
          A_DENSE( placement_row, order * lnode + i ) -= 0.5 * p_EQUATION -> jacobian()( var, i );
          A_DENSE( placement_row, order * rnode + i ) -= 0.5 * p_EQUATION -> jacobian()( var, i );
        }
        // RHS
        for ( unsigned i = 0; i < order; ++i )
        {
          B_DENSE( placement_row, order * lnode + i ) -= 0.5 * p_EQUATION -> matrix1()( var, i );
          B_DENSE( placement_row, order * rnode + i ) -= 0.5 * p_EQUATION -> matrix1()( var, i );
        }
      }
      // increment the row
      row += order;
    }

    // update the BC residuals for the current iteration
    p_RIGHT_RESIDUAL -> update( temp_dummy );
    // add the (linearised) LHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_RIGHT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        A_DENSE( row, order * ( num_of_nodes - 1 ) + var ) = p_RIGHT_RESIDUAL -> jacobian()( i, var );
        B_DENSE( row, order * ( num_of_nodes - 1 ) + var ) = 0.0;
      }
      ++row;
    }

    if ( row != num_of_nodes * order )
    {
      std::string problem( "\n The ODE_BVP has an incorrect number of boundary conditions. \n" );
      throw Error( problem );
    }

  }


  // Templated versions
  template class ODE_EVP<double>;
  //template class ODE_EVP<std::complex<double> >;

} // End of namespace TSL
