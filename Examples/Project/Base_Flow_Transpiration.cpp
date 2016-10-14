#include <cassert>
#include <cmath>

#include "Core"


// ODE enumeration
enum{ f, fd, fdd, g, gd, gdd };

namespace TSL
{
    namespace Base_Flow
    {
        double K( 0.0 );                                    // Transpiration parameter (+ve is blowing)
        double n( 0.0 );                                    // Pressure gradient parameter
        // 2D parameters
        double K_start_2D( -1.0 );
        double K_end_2D( 0.876 );
        std::size_t K_steps_2D( 30 );                      
        // 3D parameters
        double K_start_3D( -1.0 );
        double K_end_3D( 1.18 );
        std::size_t K_steps_3D( 30 );                    

	    class Blasius : public Equation<double>
	    {
		    public:
			    // The Blasius equation is 3rd order
			    Blasius() : Equation<double> ( 3 ) {} 

			    // Define the equation
			    void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			    {
				    F[ f ]   = u[ fd ];
				    F[ fd ]  = u[ fdd ];
				    F[ fdd ] = -(n + 1.0) * u[ f ] * u[ fdd ] + 2 * n * ( u[ fd ] * u[ fd ] - 1.0 );    
			    }
	    };

        class plate_BC : public Residual<double>
        {
            public:
                plate_BC() : Residual<double> ( 2, 3 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ f ] + K / ( n + 1.0 );
                B[ 1 ] = z[ fd ];
            }
        };

        class far_BC : public Residual<double>
        {
            public:
                far_BC() : Residual<double> ( 1, 3 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ fd ] - 1.0;
            }
        };

        class Alternative : public Equation<double>
	    {
		    public:
			    // The equation is 6th order
			    Alternative() : Equation<double> ( 6 ) {} 

			    // Define the equation
			    void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			    {
				    F[ f ]   = u[ fd ];
				    F[ fd ]  = u[ fdd ];
				    F[ fdd ] = - ( (n + 1.0) * u[ f ] + (1.0 - n) * u[ g ] ) * u[ fdd ] + 2.0 * n * ( u[ fd ] * u[ fd ] - 1.0 );
                    F[ g ]   = u[ gd ];
                    F[ gd ]  = u[ gdd ];
                    F[ gdd ] = - ( (n + 1.0) * u[ f ] + (1.0 - n) * u[ g ] ) * u[ gdd ] - 2.0 * n * ( 1.0 - u[ gd ] * u[ fd ] )
                               - 2.0 * u[ gd ] * u[ fd ] - (n - 1.0) * u[ gd ] * u[ gd ];     
			    }
	    };

        class alt_plate_BC : public Residual<double>
        {
            public:
                alt_plate_BC() : Residual<double> ( 4, 6 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ f ] + K;
                B[ 1 ] = z[ fd ];
                B[ 2 ] = z[ g ];
                B[ 3 ] = z[ gd ];
            }
        };

        class alt_far_BC : public Residual<double>
        {
            public:
                alt_far_BC() : Residual<double> ( 2, 6 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ fd ] - 1.0;
                B[ 1 ] = z[ gd ];
            }
        };
 
    } // End of namespace Base_Flow
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{
    cout << "----- Base flow similarity solutions for a flat plate boundary layer with transpiration -----" << endl;	
 	
	// Define the domain
	double Inf( 10.0 );											        // Infinite boundary 
	size_t N_nodes( 1000 );                         // Number of nodes
	TSL::Vector<double> nodes;								      // Declare vector of nodes (uniform)
	nodes.linspace(0,Inf,N_nodes); 

    cout << "*** Pressure gradient parameter n = " << Base_Flow::n << endl;
	
    /* ----- Solve the 2D system ----- */

    // Vector of K_values
    Vector<double> K_values_2D;
    K_values_2D.linspace( Base_Flow::K_start_2D, Base_Flow::K_end_2D, Base_Flow::K_steps_2D + 1 );

    // Mesh for storing the values of U'(eta=0) for different values of K
    OneD_node_mesh<double> K_mesh_2D( K_values_2D, 1 );

    // Create instances of the 2D equation and BCs
    Base_Flow::Blasius equation_2D;
    Base_Flow::plate_BC left_BC_2D;
    Base_Flow::far_BC right_BC_2D;

    // Create boundary value problem
    ODE_BVP<double> bvp_2D( &equation_2D, nodes, &left_BC_2D, &right_BC_2D );      

    // Set the initial guess 
	  for (std::size_t j=0; j < N_nodes; ++j )
	  {
		  double eta = nodes[ j ];					                      // eta value at node j
		  bvp_2D.solution()( j , f )  		= eta + exp( -eta );
      bvp_2D.solution()( j , fd ) 		= 1.0 - exp( -eta ); 
		  bvp_2D.solution()( j , fdd )    = exp( -eta );
	  }

    // Solve once for a good initial guess
    bvp_2D.solve_bvp();                                                                 

    cout << "*** For the 2D Blasius equation " << endl;

    // Solve the system for different values of K
    for (std::size_t i=0; i < K_values_2D.size(); ++i )
    {
      Base_Flow::K = K_values_2D[ i ];                    // Update the value of K
      bvp_2D.solve_bvp();                                     // Solve the system
      K_mesh_2D( i, 0 ) = bvp_2D.solution()( 0, fdd );    // Put the solution into the mesh
      cout << "K = " << Base_Flow::K << ", U'(eta=0) =" << bvp_2D.solution()( 0, fdd ) << endl;
    }

    K_mesh_2D.output( "./DATA/Transpiration_shear_values_2D.dat" );  // Output the solution
    

    /* ----- Solve the 3D system ----- */

    // Vector of K_values
    Vector<double> K_values_3D;
    K_values_3D.linspace( Base_Flow::K_start_3D, Base_Flow::K_end_3D, Base_Flow::K_steps_3D + 1 );

    // Mesh for storing the values of U'(eta=0) for different values of K
    OneD_node_mesh<double> K_mesh_3D( K_values_3D, 1 );

    // Create instances of the 3D equation and BCs
    Base_Flow::Alternative equation_3D;
    Base_Flow::alt_plate_BC left_BC_3D;
    Base_Flow::alt_far_BC right_BC_3D;

    // Create boundary value problem
    ODE_BVP<double> bvp_3D( &equation_3D, nodes, &left_BC_3D, &right_BC_3D );      
   
    // Set the initial guess 
	  for (std::size_t j=0; j < N_nodes; ++j )
	  {
		  double eta = nodes[ j ];					// eta value at node j
		  bvp_3D.solution()( j , f )  	= eta + exp( -eta );
      bvp_3D.solution()( j , fd ) 	= 1.0 - exp( -eta ); 
		  bvp_3D.solution()( j , fdd )  = exp( -eta );
      bvp_3D.solution()( j , g )  	= 0.35 * (1.0 - exp( -eta ));
      bvp_3D.solution()( j , gd ) 	= 1 - exp( -eta ) - exp( -1 / (eta * eta) ); 
		  bvp_3D.solution()( j , gdd )  = exp( -eta ) - 0.5 * tanh( eta ) + 0.5 * tanh( eta - 2.0 );
	  }

    // Solve once for a good initial guess
    Base_Flow::K = 0.0;
    bvp_3D.solve_bvp();

    cout << "*** For the 3D Alternative equation " << endl;

    // Solve the system for different values of K
    for (std::size_t i=0; i < K_values_3D.size(); ++i )
    {
        Base_Flow::K = K_values_3D[ i ];                 // Update the value of K
        bvp_3D.solve_bvp();                                  // Solve the system
        K_mesh_3D( i, 0 ) = bvp_3D.solution()( 0, fdd ); // Put the solution into the mesh
        cout << "K = " << Base_Flow::K << ", U'(eta=0) =" << bvp_3D.solution()( 0, fdd ) << endl;
    }

    K_mesh_3D.output( "./DATA/Transpiration_shear_values_3D.dat" ); // Output the solution

    cout << "FINISHED" << endl;
}
