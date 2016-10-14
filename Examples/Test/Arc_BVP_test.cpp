// Test the ODE_BVP class
#include "Core"

// ODE enumeration
enum{ f, fd, fdd };

namespace TSL
{
  namespace Example
  {
    double KB( 0.0 );

	  class test_equation : public Equation<double>
	  {
		  public:
        double beta;

			  // The test equation is 3rd order
			  test_equation() : Equation<double> ( 3 ) {} 

			  // Define the equation
			  void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			  {
				  F[ f ]   = u[ fd ];
				  F[ fd ]  = u[ fdd ];
				  F[ fdd ] = - u[ f ] * u[ fdd ] - beta * ( 1.0 - u[ fd ] * u[ fd ] );     
			  }
	  }; // End of class test_equation

    class plate_BC : public Residual<double>
    {
        public:
            plate_BC() : Residual<double> ( 2, 3 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ f ] + KB;
            B[ 1 ] = z[ fd ];
        }
    }; // End of class plate_BC

    class free_BC : public Residual<double>
    {
        public:
            free_BC() : Residual<double> ( 1, 3 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ fd ] - 1.0;
        }
    }; // End of class free_BC
  } // End of namespace Example 
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{
  cout << "----- TESTING ODE_BVP (arc-length continuation) -----" << endl;

  /* ----- TESTING ODE_BVP class ----- */	

	double Inf( 10.0 );									  // Infinite boundary 
	size_t N_nodes( 1000 );
	TSL::Vector<double> nodes;						// Declare vector of nodes ( uniformly spaced )
	nodes.linspace(0,Inf,N_nodes); 
	

  // Create instances of the equation and BCs
  Example::test_equation equation;
  Example::plate_BC left_BC;
  Example::free_BC right_BC;

	// Create boundary value problem
  ODE_BVP<double> ode( &equation, nodes, &left_BC, &right_BC ); 
	
    /* ----- Set the initial guess ----- */
	for (std::size_t j=0; j < N_nodes; ++j )
	{
		double eta = nodes[ j ];					// eta value at node j
		ode.solution()( j , f )  		= eta + exp( -eta );
    ode.solution()( j , fd ) 		= 1.0 - exp( -eta ); 
		ode.solution()( j , fdd )  	= exp( -eta );
	}

  ode.solve_bvp();                        // Solve the system numerically

  ode.solution().output( "./DATA/Solution_mesh_test.dat" );

  // Output initial solution
  cout << "beta = 0" << ", U'(0) = " << ode.solution()( 0, fdd ) << endl;  

  // Arc-length continuation

  // Instantiate the problem
  Example::test_equation arc_eqn;
  arc_eqn.beta = 0.0;

  // ODE boundary value problem
  TSL::ODE_BVP<double> ode_bvp_arc( &arc_eqn, nodes, &left_BC, &right_BC );

  // Initialise
  for (std::size_t j=0; j < N_nodes; ++j )
	{
		double eta = nodes[ j ];					// eta value at node j
		ode_bvp_arc.solution()( j , f )  		= eta + exp( -eta );
    ode_bvp_arc.solution()( j , fd ) 		= 1.0 - exp( -eta ); 
		ode_bvp_arc.solution()( j , fdd )  	= exp( -eta );
	}
  
  ode_bvp_arc.init_arc( &arc_eqn.beta, 0.01, 0.1 );

  // Arc-length solve the system
  double arc_step( 0.01 );
  // Output initial solution
  cout << "beta = " << arc_eqn.beta << ", U'(0) = " << ode_bvp_arc.solution()( 0, fdd );
  cout << endl; 
  do
  {
    arc_step = ode_bvp_arc.arclength_solve( arc_step );
    cout << "beta = " << arc_eqn.beta << ", U'(0) = " << ode_bvp_arc.solution()( 0, fdd );
    cout << endl; 
  }while( arc_eqn.beta < 1.0 );

	cout << "FINISHED" << endl;

}
