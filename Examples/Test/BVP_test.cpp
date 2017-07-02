// Test the ODE_BVP class
#include "Core"

// ODE enumeration
enum{ y, yd };

namespace TSL
{
	class test_equation : public Equation<double>
	{
    
		public:
			// The test equation is 2nd order ( Airy equation y'' - x*y = 0 )
			test_equation() : Equation<double> ( 2 ) {} 

			// Define the equation
			void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			{
        double x( coord( 0 ) );    
				F[ y ]   = u[ yd ];
				F[ yd ]  = x * u[ y ];   
			}
	};

    class plate_BC : public Residual<double>
    {
        public:
            plate_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ y ];
        }
    };

    class free_BC : public Residual<double>
    {
        public:
            free_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ y ] - 1.0;
        }
    }; 

    class test_equation_2 : public Equation<double>
	{
    
		public:
			// The test equation is 2nd order ( Airy equation y'' - sin(y) = 0 )
			test_equation_2() : Equation<double> ( 2 ) {} 

			// Define the equation
			void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			{   
				F[ y ]   = u[ yd ];
				F[ yd ]  = sin( u[ y ] );   
			}
	};

    class zero_BC : public Residual<double>
    {
        public:
            zero_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ y ];
        }
    };

    class one_BC : public Residual<double>
    {
        public:
            one_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ y ] - 3.14159;
        }
    }; 
} 

using namespace std;
using namespace TSL;

int main()
{
  cout << "----- TESTING ODE_BVP -----" << endl;

  /* ----- TESTING ODE_BVP class ----- */	

	double Inf( 10.0 );									  // Infinite boundary 
	size_t N_nodes( 20000 );
	TSL::Vector<double> nodes;						// Declare vector of nodes ( uniformly spaced )
	nodes.linspace(0,Inf,N_nodes); 
	
  // Create instances of the equation and BCs
  test_equation equation;
  plate_BC left_BC;
  free_BC right_BC;
  

	// Create boundary value problem
  ODE_BVP<double> ode( &equation, nodes, &left_BC, &right_BC ); 

  ode.max_iterations() = 50;
	
    /* ----- Set the initial guess ----- */
	for (std::size_t j=0; j < N_nodes; ++j )
	{
		double x = nodes[ j ];					// eta value at node j
		ode.solution()( j , y )  		= x;
    ode.solution()( j , yd ) 		= 1.; 
	}

  ode.solve_bvp();                        // Solve the system numerically

  ode.solution().output( "./DATA/Solution_mesh_test.dat" ); 

  double x = 8;
  Vector<double> soln( ode.solution().get_interpolated_vars( x ) );

  cout << "y(x=8) = " << soln[ y ] << endl;

  /* ------ y'' - sin(y) = 0 ----- */

  size_t N( 200 );

  nodes.linspace(0,1,N);
  test_equation_2 equation_2;
  zero_BC left_BC_2;
  one_BC right_BC_2;

  ODE_BVP<double> ode_2( &equation_2, nodes, &left_BC_2, &right_BC_2 );
  ode_2.max_iterations() = 50;

  for (std::size_t j=0; j < N; ++j )
	{
		double x = nodes[ j ];					// eta value at node j
		ode_2.solution()( j , y )  		= x;
    ode_2.solution()( j , yd ) 		= 1.; 
	}

  ode_2.solve_bvp();                        // Solve the system numerically

  ode_2.solution().output( "./DATA/Solution_test.dat" ); 

  double x_2 = 1;
  Vector<double> soln_2( ode_2.solution().get_interpolated_vars( x_2 ) );

  cout << "y(x=1) = " << soln_2[ y ] << endl;

    

	cout << "FINISHED" << endl;

}
