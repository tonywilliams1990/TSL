/// Implementation of a class for an \f$ n^{th} \f$-order IBVP of the form
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

#include <string>
#include <cassert>

#include "PDE_IBVP.h"
#include "Utility.h"

namespace TSL
{

  template <class T>
  PDE_IBVP<T>::PDE_IBVP( Equation_2matrix<T >* ptr_to_equation,
                             const Vector<double> &nodes,
                             Residual_with_coords<T>* ptr_to_left_residual,
                             Residual_with_coords<T>* ptr_to_right_residual ) :
      TOL( 1.e-8 ),
      TIME( 0.0 ),
      MAX_ITERATIONS( 12 ),
      p_EQUATION( ptr_to_equation ),
      p_LEFT_RESIDUAL( ptr_to_left_residual ),
      p_RIGHT_RESIDUAL( ptr_to_right_residual )
  {
    SOLN = OneD_node_mesh<T>( nodes, p_EQUATION -> get_order() );
    PREV_SOLN = SOLN;
    p_EQUATION -> coord( 1 ) = TIME;
    if ( p_EQUATION -> get_order() - p_LEFT_RESIDUAL -> get_order() - p_RIGHT_RESIDUAL -> get_order() != 0 )
    {
      std::string problem( "\n The PDE_IBVP class has been constructed, but the order of the \n");
      problem += "system does not match the number of boundary conditions.\n";
     throw Error( problem );      
    }

  }

  template <class T>
  PDE_IBVP<T>::~PDE_IBVP()
  {

  }

  template <class T>
  void PDE_IBVP<T>::step2( const double& dt )
  {
    // the order of the problem
    unsigned order( p_EQUATION -> get_order() );
    // get the number of nodes in the mesh
    // -- this may have been refined by the user since the last call.
    unsigned ny( SOLN.get_nnodes() );
    // measure of maximum residual
    double max_residual( 1.0 );
    // iteration counter
    int counter = 0;
    // store this soln as the 'previous SOLUTION'
    PREV_SOLN = SOLN;
    // ANY LARGE STORAGE USED IN THE MAIN LOOP IS
    // DEFINED HERE TO AVOID REPEATED CONSTRUCTION.
    // Note we blank the A matrix after every iteration.
    //
    // Banded LHS matrix - max obove diagonal band width is
    // from first variable at node i to last variable at node i+1
    SparseMatrix<T> a( ny * order, ny * order );
    // RHS
    Vector<T> b( ny * order, 0.0 );
    // loop until converged or too many iterations
    do
    {
      // iteration counter
      ++counter;

      assemble_matrix_problem( a, b, dt );
      max_residual = b.norm_inf();
#ifdef DEBUG
      std::cout << " PDE_IBVP.solve : Residual_max = " << max_residual << " tol = " << TOL << "\n";
#endif

      // linear solver
      b = a.solve( b );
      //system.solve();
      // keep the solution in a OneD_GenMesh object
      for ( std::size_t var = 0; var < order; ++var )
      {
        for ( std::size_t i = 0; i < ny; ++i )
        {
          SOLN( i, var ) += b[ i * order + var ];
        }
      }

    }
    while ( ( max_residual > TOL ) && ( counter < MAX_ITERATIONS ) );
    if ( max_residual > TOL ) 
    {
      // restore to previous state because this step failed
      SOLN = PREV_SOLN;
      std::string problem( "\n The PDE_IBVP.step2 method took too many iterations.\n");
      problem += "Solution has been restored to the previous accurate state.\n";
      throw Error( problem );
    }
    // set the time to the updated level
    TIME += dt;
  }


  template <class T>
  void PDE_IBVP<T>::assemble_matrix_problem( SparseMatrix<T>& a, 
                                                 Vector<T>& b, const double& dt )
  {
    // clear the Jacobian matrix
    a.clear( );
    // inverse of the time step
    const double inv_dt( 1. / dt );
    // the order of the problem
    const unsigned order( p_EQUATION -> get_order() );
    // number of spatial nodes
    const unsigned ny( SOLN.get_nnodes() );
    // row counter
    std::size_t row( 0 );
    // a matrix that is used in the Jacobian of the mass matrix terms
    Matrix<T> h0( order, order, 0.0 );
    Matrix<T> h1( order, order, 0.0 );
    // local state variable and functions
    Vector<T> F_midpt( order, 0.0 );
    Vector<T> O_midpt( order, 0.0 );
    Vector<T> state( order, 0.0 );
    Vector<T> state_dt( order, 0.0 );
    Vector<T> state_dy( order, 0.0 );
    // BCn equation is evaluated at the next time step
    p_LEFT_RESIDUAL -> coord( 0 ) = TIME + dt;
    // update the BC residuals for the current iteration
    p_LEFT_RESIDUAL -> update( SOLN.get_nodes_vars( 0 ) );
    // add the (linearised) LHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_LEFT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        a( row, var ) = p_LEFT_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - p_LEFT_RESIDUAL -> residual()[ i ];
      ++row;
    }
    // inner nodes of the mesh, node = 0,1,2,...,N-2
    for ( std::size_t node = 0; node <= ny - 2; ++node )
    {
      const std::size_t l_node = node;
      const std::size_t r_node = node + 1;
      // inverse of step size
      const double inv_dy = 1. / ( SOLN.coord( r_node ) - SOLN.coord( l_node ) );
      // set the current solution at this node by 2nd order evaluation at mid point
      for ( unsigned var = 0; var < order; ++var )
      {
        const T F_midpt = ( SOLN( l_node, var ) + SOLN( r_node, var ) ) / 2.;
        const T O_midpt = ( PREV_SOLN( l_node, var ) + PREV_SOLN( r_node, var ) ) / 2.;
        state_dy[ var ] = ( SOLN( r_node, var ) - SOLN( l_node, var )
                            + PREV_SOLN( r_node, var ) - PREV_SOLN( l_node, var ) ) * inv_dy / 2.;
        state[ var ] = ( F_midpt + O_midpt ) / 2.;
        state_dt[ var ] = ( F_midpt - O_midpt ) * inv_dt;
      }
      // set the equation's y & t values to be mid points
      // y
      p_EQUATION -> coord(0) = 0.5 * ( SOLN.coord( l_node ) + SOLN.coord( r_node ) );
      // t
      p_EQUATION -> coord(1) = TIME + dt / 2;
      // Update the equation to the mid point position
      p_EQUATION -> update( state );
      // evaluate the Jacobian of mass contribution multiplied by state_dy
      p_EQUATION -> get_jacobian_of_matrix0_mult_vector( state, state_dy, h0 );
      // evaluate the Jacobian of mass contribution multiplied by state_dt
      p_EQUATION -> get_jacobian_of_matrix1_mult_vector( state, state_dt, h1 );
      // loop over all the variables
      //
      for ( unsigned var = 0; var < order; ++var )
      {
       // add the matrix mult terms to the linearised problem
        for ( unsigned i = 0; i < order; ++i )  // dummy index
        {
          // add the Jacobian terms
          // a( row, order * l_node + i ) -= jac_midpt( var, i ) * 0.5;
          a( row, order * l_node + i ) -= p_EQUATION -> jacobian()( var, i ) * 0.5;
          // a( row, order * r_node + i ) -= jac_midpt( var, i ) * 0.5;
          a( row, order * r_node + i ) -= p_EQUATION -> jacobian()( var, i ) * 0.5;
          // add the Jacobian of mass terms for dt terms
          // a( row, order * l_node + i ) += h1( var, i ) * 0.5;
          a( row, order * l_node + i ) += h0( var, i ) * 0.5;
          // a( row, order * r_node + i ) += h1( var, i ) * 0.5;
          a( row, order * r_node + i ) += h0( var, i ) * 0.5;
          // add the Jacobian of mass terms for dt terms
          // a( row, order * l_node + i ) += h2( var, i ) * 0.5;
          a( row, order * l_node + i ) += h1( var, i ) * 0.5;
          // a( row, order * r_node + i ) += h2( var, i ) * 0.5;
          a( row, order * r_node + i ) += h1( var, i ) * 0.5;
          // add the mass matrix terms
          // a( row, order * l_node + i ) -= mass_midpt( var, i ) * inv_dy;
          a( row, order * l_node + i ) -= p_EQUATION -> matrix0()( var, i ) * inv_dy;
          // a( row, order * r_node + i ) += mass_midpt( var, i ) * inv_dy;
          a( row, order * r_node + i ) += p_EQUATION -> matrix0()( var, i ) * inv_dy;
          // a( row, order * l_node + i ) += mass_midpt( var, i ) * inv_dt;
          a( row, order * l_node + i ) += p_EQUATION -> matrix1()( var, i ) * inv_dt;
          // a( row, order * r_node + i ) += mass_midpt( var, i ) * inv_dt;
          a( row, order * r_node + i ) += p_EQUATION -> matrix1()( var, i ) * inv_dt;
        }
        // RHS
        b[ row ] = p_EQUATION -> residual()[ var ];

        Vector<T> mat0vec;
        Matrix<T> mat0;
        mat0 = p_EQUATION -> matrix0();
        mat0vec = mat0[ var ];
        b[ row ] -= Utility::dot( mat0vec, state_dy );

        Vector<T> mat1vec;
        Matrix<T> mat1;
        mat1 = p_EQUATION -> matrix1();
        mat1vec = mat1[ var ];
        b[ row ] -= Utility::dot( mat1vec, state_dt );

        b[ row ] *= 2;
        // increment the row
        row += 1;
      }
    }
    // BCn equation is evaluated at the next step time point as
    // they cannot depend on d/dt terms
    p_RIGHT_RESIDUAL -> coord( 0 ) = TIME + dt;
    // update the BC residuals for the current iteration
    p_RIGHT_RESIDUAL -> update( SOLN.get_nodes_vars( ny - 1 ) );
    // add the (linearised) RHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_RIGHT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at RHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        a( row, order * ( ny - 1 ) + var ) = p_RIGHT_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - p_RIGHT_RESIDUAL -> residual()[ i ];
      ++row;
    }
#ifdef PARANOID
    if ( row != ny * order )
    {
      std::string problem( "\n The ODE_BVP has an incorrect number of boundary conditions. \n" );
      throw Error( problem );
    }
#endif
  }

  // the templated versions that we require are:
  template class PDE_IBVP<double>
  ;
  template class PDE_IBVP< std::complex<double> >
  ;

} // end namespace

