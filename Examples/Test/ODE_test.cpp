// Test the ODE_BVP class
#include "Core"

// ODE enumeration
enum{ f, fd, fdd };

namespace TSL
{
	class test_equation : public Equation<double>
	{
		public:
			// The test equation is 3rd order
			test_equation() : Equation<double> ( 3 ) {} 

			// Define the equation
			void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			{
				F[ f ]   = u[ fd ];
				F[ fd ]  = u[ fdd ];
				F[ fdd ] = - u[ f ] * u[ fdd ];    
			}
	};

    class plate_BC : public Residual<double>
    {
        public:
            plate_BC() : Residual<double> ( 2, 3 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ f ];
            B[ 1 ] = z[ fd ];
        }
    };

    class free_BC : public Residual<double>
    {
        public:
            free_BC() : Residual<double> ( 1, 3 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ fd ] - 1.0;
        }
    }; 
} 

using namespace std;
using namespace TSL;

int main()
{
  cout << "----- TESTING ODE_BVP -----" << endl;

  /* ----- TESTING ODE_BVP class ----- */	

  /* ----- Test Tode_bvp class ----- */
	double Inf( 10.0 );									  // Infinite boundary 
	size_t N_nodes( 1000 );
	TSL::Vector<double> nodes;						// Declare vector of nodes ( uniformly spaced )
	nodes.linspace(0,Inf,N_nodes); 
	

  // Create instances of the equation and BCs
  test_equation equation;
  plate_BC left_BC;
  free_BC right_BC;

	// Create boundary value problem
  ODE_BVP<double> ode( &equation, nodes, &left_BC, &right_BC ); 

  ODE_BVP<double> ode_copy( ode );      // Test copy construction
	
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

  // Test delta() method
  cout << "DELTA = " << ode.delta() << endl;
  ode.delta() = 1.0e-10;
  cout << "DELTA = " << ode.delta() << endl;

  // Test max_iterations() method
  cout << "MAX_ITER = " << ode.max_iterations() << endl;
  ode.max_iterations() = 50;
  cout << "MAX_ITER = " << ode.max_iterations() << endl;

  // Test tolerance() method
  cout << "TOL = " << ode.tolerance() << endl;
  ode.tolerance() = 1.0e-12;
  cout << "TOL = " << ode.tolerance() << endl;

	cout << "FINISHED" << endl;

}
