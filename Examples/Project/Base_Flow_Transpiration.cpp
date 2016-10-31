#include <cassert>
#include <cmath>

#include "Core"


// ODE enumeration
enum{ f, fd, fdd, g, gd, gdd };

namespace TSL
{
  namespace Base_Flow
  {
    double KB( 2.0 );                      // Transpiration parameter (+ve is blowing)
                       
    class Falkner : public Equation<double>
    {
      public:
        double beta;

        // The Falkner-Skan equation is 3rd order
        Falkner() : Equation<double> ( 3 ) {} 

        // Define the equation
        void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
        {
          F[ f ]   = u[ fd ];
          F[ fd ]  = u[ fdd ];
          F[ fdd ] = - u[ f ] * u[ fdd ] - beta * ( 1.0 - u[ fd ] * u[ fd ] );    
        }
    };

    class plate_BC : public Residual<double>
    {
      public:
        plate_BC() : Residual<double> ( 2, 3 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
          B[ 0 ] = z[ f ] + KB;
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
          double beta;

			    // The equation is 6th order
			    Alternative() : Equation<double> ( 6 ) {} 

			    // Define the equation
			    void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			    {
				    F[ f ]    =  u[ fd ];
            F[ fd ]   =  u[ fdd ];
            F[ fdd ]  = -( u[ f ] + ( 2.0 - beta ) * u[ g ] ) * u[ fdd ] 
                        - beta * ( 1.0 - u[ fd ] * u[ fd ] );
            F[ g ]    =  u[ gd ];
            F[ gd ]   =  u[ gdd ];
            F[ gdd ]  = -( u[ f ] + ( 2.0 - beta ) * u[ g ] ) * u[ gdd ]
                        -( 2.0 * ( 1.0 - beta ) * u[ fd ] 
                        - ( 2.0 - beta) * u[ gd ] ) * u[ gd ];    
			    }
    };

    class alt_plate_BC : public Residual<double>
    {
      public:
        double K_alt;                // Set K=0 then iterate to K=KB

        alt_plate_BC() : Residual<double> ( 4, 6 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
          B[ 0 ] = z[ f ] + K_alt;
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
  cout << "----- Base flow transpiration arc-length continuation -----" << endl;	
 	
  double Inf( 20.0 );									  // Infinite boundary   
  size_t N_nodes( 1000 );
	TSL::Vector<double> nodes;						// Declare vector of nodes
	//nodes.linspace( 0.0, Inf, N_nodes ); 
  nodes.power( 0.0, Inf, N_nodes, 2.0 );

  /* ----- 2D problem ----- */
  
  // Create instances of the equation and BCs
  Base_Flow::Falkner  equation;
  Base_Flow::plate_BC left_BC;
  Base_Flow::far_BC   right_BC;
  equation.beta = 1.01; 
 
  // Create boundary value problem
  ODE_BVP<double> ode_2D( &equation, nodes, &left_BC, &right_BC );

  // Set the initial guess
	for (std::size_t j=0; j < N_nodes; ++j )
	{
		double eta = nodes[ j ];					// eta value at node j
		ode_2D.solution()( j , f )  		= eta + exp( -eta );
    ode_2D.solution()( j , fd ) 		= 1.0 - exp( -eta ); 
		ode_2D.solution()( j , fdd )  	= exp( -eta );
	}
  // Initialise the arc-length continuation 
  double arc_step( -0.01 );
  double max_step( 0.025 );
  ode_2D.init_arc( &equation.beta, arc_step, max_step ); 

  // Output initial solution
  double stress( ode_2D.solution()( 0, fdd ) );
  cout << "beta = " << equation.beta << ", U'(0) = "<< stress << endl;
  
  // Vectors for storing beta and U'(0) values
  Vector<double> betas, shear;
  betas.push_back( equation.beta );
  shear.push_back( stress );

  // Iterate 
  do
  {
    arc_step = ode_2D.arclength_solve( arc_step );
    stress = ode_2D.solution()( 0, fdd );
    cout << "beta = " << equation.beta << ", U'(0) = "<< stress << endl;
    betas.push_back( equation.beta );
    shear.push_back( stress );
  }while( std::abs( equation.beta ) > 0.005 || stress > 0.0 );

  // Create mesh for output
  OneD_node_mesh<double> Solution_2D( betas, 1 );
  Solution_2D.set_vars_from_vector( shear );
  Solution_2D.output( "./DATA/Solution_2D.dat" );
  
  /* ----- 3D problem ----- */
  
  // Create instances of the equation and BCs
  Base_Flow::Alternative  alt_eqn;
  Base_Flow::alt_plate_BC alt_left_BC;
  Base_Flow::alt_far_BC   alt_right_BC;
  alt_eqn.beta = 0.1;
  alt_left_BC.K_alt = 0.0;

  // Create the BVP
  ODE_BVP<double> ode_3D( &alt_eqn, nodes, &alt_left_BC, &alt_right_BC );

  // Set the initial guess
  for (std::size_t j=0; j < N_nodes; ++j )
	{
		double eta = nodes[ j ];					                   // eta value at node j
		ode_3D.solution()( j, f )  	= eta + exp( -eta ) - 1.0;
    ode_3D.solution()( j, fd ) 	= 1.0 - exp( -eta ); 
		ode_3D.solution()( j, fdd ) = exp( -eta );
    ode_3D.solution()( j, g )  	= 0.35 * (1.0 - exp( -eta ));
    ode_3D.solution()( j, gd ) 	= 1.0 - exp( -eta ) - exp( -1.0 / (eta * eta) ); 
		ode_3D.solution()( j, gdd ) = exp( -eta ) - 0.5 * tanh( eta ) + 0.5 * tanh( eta -2.0 );
	}

  // Solve the system for K=0 then arc-length continue until K=KB
  ode_3D.init_arc( &alt_left_BC.K_alt, 0.01, 0.1 );
  cout << "K_alt = " << alt_left_BC.K_alt << ", U'(0) = ";
  cout << ode_3D.solution()( 0, fdd ) << endl;
  arc_step = 0.01;
  do
  {
    arc_step = ode_3D.arclength_solve( arc_step );
    cout << "K_alt = " << alt_left_BC.K_alt << ", U'(0) = ";
    cout << ode_3D.solution()( 0, fdd ) << endl;

  }while( alt_left_BC.K_alt < Base_Flow::KB ); 

  alt_left_BC.K_alt = Base_Flow::KB;
  ode_3D.solve_bvp();                             // Solve once more with K=KB
  cout << "K_alt = " << alt_left_BC.K_alt << ", U'(0) = ";
  cout << ode_3D.solution()( 0, fdd ) << endl;

  cout << "KB = " << Base_Flow::KB << endl;
  cout << "K_alt = " << alt_left_BC.K_alt << endl;

  arc_step = -0.01;
  ode_3D.init_arc( &alt_eqn.beta, arc_step, max_step );

  // Output initial solution
  cout << " ----- Alternative equation ----- " << endl;
  double stress_3D( ode_3D.solution()( 0, fdd ) );
  cout << "beta = " << alt_eqn.beta << ", U'(0) = "<< stress_3D << endl;

  arc_step = 0.01;

  // Arc-length solve to beta > 1.0
  do
  {
    arc_step = ode_3D.arclength_solve( arc_step );
    stress_3D = ode_3D.solution()( 0, fdd );
    cout << "beta = " << alt_eqn.beta << ", U'(0) = "<< stress_3D << endl;

  }while( alt_eqn.beta < 1.0 );

  cout << " ----- beta = approx 1 ----- " << endl;

  // Vectors for storing beta and U'(0) values
  Vector<double> betas_3D, shear_3D;
  betas_3D.push_back( alt_eqn.beta );
  shear_3D.push_back( stress_3D );
  
  arc_step = -0.01;
  // Iterate 
  do
  {
    arc_step = ode_3D.arclength_solve( arc_step );
    stress_3D = ode_3D.solution()( 0, fdd );
    cout << "beta = " << alt_eqn.beta << ", U'(0) = "<< stress_3D << endl;
    betas_3D.push_back( alt_eqn.beta );
    shear_3D.push_back( stress_3D );
  }while( std::abs( alt_eqn.beta ) > 0.005 || stress_3D > 1e-5 );

  // Create mesh for output
  OneD_node_mesh<double> Solution_3D( betas_3D, 1 );
  Solution_3D.set_vars_from_vector( shear_3D );
  Solution_3D.output( "./DATA/Solution_3D.dat" );

  cout << "FINISHED" << endl;
}
