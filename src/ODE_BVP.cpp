/* ODE_BVP - Implementation for solving ODE boundary value problems
*/
#include "ODE_BVP.h"
#include "Error.h"
#include <complex>

namespace TSL
{
	/* ----- Constructors and destructors ----- */

	// Constructor 
	template <typename T, typename X>
    ODE_BVP<T, X>::ODE_BVP( Equation<T, X>* ptr_equation, const Vector<X>& nodes, Residual<T>* ptr_left_residual, 
                      Residual<T>* ptr_right_residual ) :
        MAX_ITER( 20 ), TOL( 1.0e-8 ), ptr_EQUATION( ptr_equation ), NODES( nodes ), DELTA( 1.0e-8 ),
        ptr_LEFT_RESIDUAL( ptr_left_residual ), ptr_RIGHT_RESIDUAL( ptr_right_residual )
    {
        // Set up the solution mesh  
        SOLUTION = OneD_node_mesh<T, X>( NODES, ptr_EQUATION -> get_order() );

        // Check the order of the system and the boundary conditions are consistent
        if ( ( ptr_LEFT_RESIDUAL -> get_number_of_vars() != ptr_EQUATION -> get_order() ) ||
         ( ptr_RIGHT_RESIDUAL -> get_number_of_vars() != ptr_EQUATION -> get_order() ) ||
         ( ptr_LEFT_RESIDUAL -> get_order() + ptr_RIGHT_RESIDUAL -> get_order() != ptr_EQUATION -> get_order() ) )
        {
         std::cout << "order " << ptr_EQUATION -> get_order() << "\n";
         std::cout << "left nvars " << ptr_LEFT_RESIDUAL -> get_number_of_vars() << "\n";
         std::cout << "right nvars " << ptr_RIGHT_RESIDUAL -> get_number_of_vars() << "\n";
         std::cout << "left order " << ptr_LEFT_RESIDUAL -> get_order() << "\n";
         std::cout << "right order " << ptr_RIGHT_RESIDUAL -> get_order() << "\n";
         std::string problem;
         problem  = "It looks like the ODE_BVP equation and boundary conditions are\n";
         problem += "not well posed. The number of variables for each boundary condition\n";
         problem += "has to be the same as the order of the equation. Also the order of\n";
         problem += "both boundary conditions has to sum to the order of the equation.\n";
         throw Error( problem );
        }
    }

	// Copy constructor
	template <typename T, typename X>
	ODE_BVP<T, X>::ODE_BVP(  const ODE_BVP& source )
	{
		*this = source;
	}

	// Destructor
	template <typename T, typename X>
  ODE_BVP<T, X>::~ODE_BVP()
	{}

	/* ----- Methods ----- */

    // Return a pointer to DELTA
    template <typename T, typename X>
    double& ODE_BVP<T, X>::delta()
    {
        return DELTA;
    }

    // Return a pointer to MAX_ITER
    template <typename T, typename X>
    std::size_t& ODE_BVP<T, X>::max_iterations()
    {
        return MAX_ITER;
    }

    // Return a pointer to TOL
    template <typename T, typename X>
    double& ODE_BVP<T, X>::tolerance()
    {
        return TOL;
    }

    // Return the current solution mesh
	  template <typename T, typename X>
	  OneD_node_mesh<T, X>& ODE_BVP<T, X>::solution() 
	  {
		  return SOLUTION;
	  }

    // Solve the system using Newton iteration (points to solve_bvp)
    template <class T, class X>
    void ODE_BVP<T, X>::solve( Vector<T>& state )
    {
      SOLUTION.set_vars_from_vector( state );
      solve_bvp();
      state = SOLUTION.vars_as_vector();
    }

    template <class T, class X>
    void ODE_BVP<T, X>::assemble_matrix_problem( SparseMatrix<T>& A, Vector<T>& B )
    {
      std::size_t order( ptr_EQUATION -> get_order() );// Get the order of the problem 
      std::size_t N( SOLUTION.get_nnodes() );          // Number of nodes in the mesh
      std::size_t Size( order * N );
      Matrix<T> J( order, order, 0.0 ); 			          // Jacobian matrix
		  Vector<T> R( order, 0.0 );						            // RHS of equation at each node
		  ODE_BVP<T,X> temp( *this );						            // Temporary copy of *this
		  Vector<T> F_fun( order, 0.0 );					          // Vector for function evaluation
		  Vector<T> F_fun_star( order, 0.0 );			          // Vector for perturbed function
		  Matrix<T> Identity( order, order, 0.0 );		      // Identity matrix
		  Matrix<T> K( order, order, 0.0 );				          // Coefficient matrix K
		  Matrix<T> L( order, order, 0.0 );				          // Coefficient matrix L

      A.clear();										                    // Reset the A matrix
			B.assign( Size, 0.0 );							              // Reset B vector
			std::size_t row( 0 );							                // Row counter

      // Fill in the LHS boundary conditions in the matrix

        // Update the BC residuals for the current iteration
        ptr_LEFT_RESIDUAL -> update( SOLUTION.get_nodes_vars( 0 ) );
			  // Add the (linearised) LHS BCs to the matrix problem
        for ( unsigned i = 0; i < ptr_LEFT_RESIDUAL -> get_order(); ++i )
        {
          for ( unsigned var = 0; var < order; ++var )
          {
            A( row, var ) = ptr_LEFT_RESIDUAL -> jacobian()( i, var );
          }
          B[ row ] = - ptr_LEFT_RESIDUAL -> residual()[ i ];
          ++row;
        }

			  // Calculate the equations at i + 1/2
			  J.fill( 0.0 );								                // Reset Jacobian
			  R.assign( order, 0.0 );						            // Reset RHS of equation

			for (std::size_t i=0; i < N-1; ++i)
			{
				T x_i = NODES[ i ];						                    // Get the nodal positions x_i
				T x_i_1 = NODES[ i + 1 ];				                  // and x_i+1
				T Dx_i = x_i_1 - x_i;					                    // Step between x_i and x_i+1 
				T x_i_half = x_i + Dx_i / 2.0; 			              // Mid-node location
			
				Vector<T> u_G_i, u_G_i_1, u_G_i_half, u_G_i_half_star;
				u_G_i = SOLUTION.get_nodes_vars( i );			        // Solution stored at node i
				u_G_i_1 = SOLUTION.get_nodes_vars( i+1 );		      // Solution stored at node i+1
				u_G_i_half = ( u_G_i + u_G_i_1 ) / 2; 	   // Linearly interpolate at i + 1/2

        // Function evaluation at i + 1/2
				F_fun.assign( order, 0.0 );				                    // Clear F_fun vector
				ptr_EQUATION->residual_fn( u_G_i_half, F_fun );

                // Jacobian at i + 1/2
				for ( std::size_t alpha = 0; alpha < order; ++alpha )	  // Row index
				{				
				  	for ( std::size_t beta = 0; beta < order; ++beta )	// Column index
					{
						Vector<T> delta( order ,0.0); 					            // Perturbation vector
						delta[ beta ] += DELTA;
						u_G_i_half_star = u_G_i_half + delta;			          // Perturb the known value

						// Perurbed function evaluation at i + 1/2
						F_fun_star.assign( order, 0.0 );				            // Clear vector	
						ptr_EQUATION->residual_fn( u_G_i_half_star, F_fun_star );

						// Approximate the entry in the Jacobian
						J(alpha,beta) = ( F_fun_star[ alpha ] - F_fun[ alpha ] ) / DELTA;
					}
				}
				R = F_fun - (u_G_i_1 - u_G_i)/Dx_i;						          // RHS vector	

        // Coefficient matrices
				Identity.fill_diag(1.0/Dx_i);
				K.fill( 0.0 );											                    // Reset matrices
				L.fill( 0.0 );
				K = Identity - J * 0.5;							                    // Fill matrices
				L = -Identity - J * 0.5;	

				// Fill the rows in the A matrix and B vector
				for (std::size_t n=0; n<order; ++n)
				{
					for (std::size_t j=0; j<order; ++j)
					{
						A(row, i*order + j ) = L(n,j);
						A(row, (i+1)*order + j ) = K(n,j);
					}
					B[ row ] = R[ n ];
					++row;
				}
			  }

        // Fill in the RHS boundary conditions in the matrix

        // Update the BC residuals for the current iteration
        ptr_RIGHT_RESIDUAL -> update( SOLUTION.get_nodes_vars( N - 1 ) );
        // Add the (linearised) RHS BCs to the matrix problem
        for ( unsigned i = 0; i < ptr_RIGHT_RESIDUAL -> get_order(); ++i )
        {
          // loop thru variables at RHS of the domain
          for ( unsigned var = 0; var < order; ++var )
          {
            A( row, order * ( N - 1 ) + var ) = ptr_RIGHT_RESIDUAL -> jacobian()( i, var );
          }
          B[ row ] = - ptr_RIGHT_RESIDUAL -> residual()[ i ];
          ++row;
        }

    }
    

    // Solve the BVP - do iterate while residual > tol and iter < max_iter
    template <typename T, typename X>
	  void ODE_BVP<T, X>::solve_bvp()
    {
      std::size_t order( ptr_EQUATION -> get_order() );// Get the order of the problem 
      std::size_t counter( 0 );                        // Initialise an iteration counter
      double max_res( 1.0 );                           // Measure of the maximum residual
      std::size_t N( SOLUTION.get_nnodes() );          // Number of nodes in the mesh

      // Declare A matrix and vector RHS
      std::size_t Size( order * N );
		  SparseMatrix<T> A( Size, Size );				          // Jacobian matrix
		  Vector<T> B( Size, 0.0 );						              // RHS vector
		  Vector<T> sol( Size, 0.0 );					              // Vector for storing solution     

      do
		  {	
			  ++counter;										                  // Increment counter 				
			  /* ----- Iterate to a solution ----- */

			  assemble_matrix_problem( A, B );                // Assemble the matrix problem
			  sol.assign( Size, 0.0 );									      // Reset solution vector
			  sol = A.solve( B ); 	                          // Solve the system of equations

			  // Update the solution in the OneD_node_mesh (SOLUTION)
			  for ( std::size_t i=0; i < N; ++i )
			  {	
				  for (std::size_t j=0; j < order; j++ )
				  {
					  SOLUTION(i,j) += sol[ i * order + j ];
				  }
			  }
			  max_res = B.norm_inf();
			
		  }while( counter < MAX_ITER && max_res > TOL );
    }	

    /// Arc-length solve the system (init_arc must be called first)
    template <class T, class X>
    double ODE_BVP<T, X>::arclength_solve( const double& step )
    { 
      this -> ds() = step;
      std::size_t order( ptr_EQUATION -> get_order() );         // Order of the equation
      std::size_t N( SOLUTION.get_nnodes() );                   // Number of nodes
      std::size_t Size( order * N );
      
      SparseMatrix<T> Jac( Size, Size );                        // Jacobian matrix
      // Residuals over all nodes
      Vector<T> Res1( Size, 0.0 );                              
      Vector<T> Res2( Size, 0.0 );
      Vector<T> dRes_dp( Size, 0.0 );
      // RHS vectors for the linear solvers
      Vector<T> y( Size, 0.0 );
      Vector<T> z( Size, 0.0 );
      Vector<T> Jac_E;
      // Make backups in case we can't find a converged solution
      Vector<T> backup_state( SOLUTION.vars_as_vector() );
      T backup_parameter( *( this -> ptr_PARAM ) );
      // Generate a 1st order guess for the next state and parameter
      Vector<T> x( this -> LAST_X + this -> X_DERIV_S * this -> DS );
      *( this -> ptr_PARAM ) = this -> LAST_PARAM + this -> PARAM_DERIV_S * this -> DS;
      
      SOLUTION.set_vars_from_vector( x );                       // Update the solution mesh
      int det_sign( 0 );                                        // Determinant monitor
      bool step_succeeded( false );                             // Check for success
      std::size_t itn( 0 );                                     // Iteration counter
      
      do
      { 
        ++itn;
        double E1 = this -> arclength_residual( x );            // Arclength residual
        assemble_matrix_problem( Jac, Res1 );                   // Assemble matrix problem
        y = Res1; //TODO do we need this?
        if ( Res1.norm_inf() < TOL && itn > 1 )
        {
          step_succeeded = true;
          break;
        }
        y = Jac.solve( Res1 );
        // Derivatives wrt parameter
        *( this -> ptr_PARAM ) += DELTA;
        assemble_matrix_problem( Jac, Res2 );
        double E2 = this -> arclength_residual( x );
        *( this -> ptr_PARAM ) -= DELTA;
        dRes_dp = ( Res2 - Res1 ) / DELTA;
        double dE_dp = ( E2 - E1 ) / DELTA;
        z = dRes_dp; //TODO again do we need this?
        z = Jac.solve( dRes_dp );
        Jac_E = this -> Jac_arclength_residual( x );
        T delta_p = -( E1 + Jac_E.dot( y ) ) / ( dE_dp + Jac_E.dot( z ) );
        Vector<T> delta_x = y + z * delta_p;
        // Update the state variables and the parameter with the corrections
        x += delta_x;
        *( this -> ptr_PARAM ) += delta_p;
        SOLUTION.set_vars_from_vector( x );
        // Check for convergence
        if ( delta_x.norm_inf() < TOL )
        {
          step_succeeded = true;
          break;
        }
        if ( itn > MAX_ITER )
        {
          step_succeeded = false;
          break;
        }
       }while( true );
       // If this isn't a successful step
       if ( !step_succeeded )
       {
          // Restore things using the backups we made earlier
          SOLUTION.set_vars_from_vector( backup_state );            // Restore state
          *( this -> ptr_PARAM ) = backup_parameter;                // Restore parameter
          this -> DS /= this -> ARCSTEP_MULTIPLIER;                 // Reduce step length
       }
       else
       {
          this -> update( SOLUTION.vars_as_vector() );
          // TODO Bifurcation detection using sign of determinant (ODE_BVP 318)
          if ( itn > 8 || std::abs( this -> DS ) > this -> MAX_DS )
          {
            // Converging too slowly, so decrease DS
            this -> DS /= this -> ARCSTEP_MULTIPLIER;
          }          
          if ( itn < 4 ) // Converging too quickly, so increase DS
          {
            if ( std::abs( this -> DS * this -> ARCSTEP_MULTIPLIER ) < this -> MAX_DS )
            {
              this -> DS *= this -> ARCSTEP_MULTIPLIER;
            }
          }

       }
       return this -> DS;
    }

  	// Templated versions
  	template class ODE_BVP<double>;
	  template class ODE_BVP< std::complex<double> >;
    template class ODE_BVP< std::complex<double>, std::complex<double> >;

} // End of namespace TSL
