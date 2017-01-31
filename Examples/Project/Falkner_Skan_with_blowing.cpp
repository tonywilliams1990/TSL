#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <sstream>

#include "Core"
#include "Sparse"

// Enumerations
enum{ f, fd, fdd };                               // Base ODE
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };                 // Base ODE 

#define BASE_2D
// Either UNIFORM or NONUNIFORM for uniform of non-uniform mesh
#define NONUNIFORM

namespace TSL
{
    namespace Param
    {
      double eta_top( 128.0 );          // Size of the domain in the eta direction
      const std::size_t M( 400 );       // Number of intervals in the eta direction
      double beta( 0.1 );               // Hartree parameter
      double KB( 0.0 );                 // Base flow transpiration ( +ve = blowing )
      double KB_max( 5.0 );             // Maximum value of the transpiration    
      std::size_t KB_n( 51 );            // Number of KB values 

    } // End of namespace Param

    namespace Example
    {
      std::string output_path;          // Output path                                     

    } // End of namespace Example


    namespace Mesh
    {
#ifdef UNIFORM    
      double Y( const double& eta )
      {
        return eta;
      }
      double Yd( const double& eta )
      {
        return 1;
      }
      double Ydd( const double& eta )
      {
        return 0;
      }
#endif
#ifdef NONUNIFORM
      const double b1( 0.3 );
      const double b2( 0.3 );   // Y = (eta_hat + b1)^b2

      double Y( const double& eta )
      {
        return std::pow(eta + b1, b2) - pow(b1,b2);
      }
      double Yd( const double& eta )
      {
        return b2 * std::pow(eta + b1, b2 - 1);
      }
      double Ydd( const double& eta )
      {
        return b2 * (b2 - 1) * std::pow(eta + b1, b2 - 2);
      }
#endif

      class invert_eta : public Residual<double>
      {
        // USED TO INVERT THE NON-UNIFORM MESH MAPPING
        public:
        double Y0;

        invert_eta() : Residual<double>( 1 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &f ) const
        {
          f[ 0 ] = Y( z[0] ) - Y0;
        }
      };
 
    } // End of namespace Mesh

    namespace Base_Flow
    {  
#ifdef BASE_2D
      class equation : public Equation<double>
      {
        public:
          double beta;                            // Hartree parameter
          // Falkner-Skan equation is 3rd order
          equation() : Equation<double> ( 3 ) {}
          // Define the equation
          void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
          {
            F[ f ]   = u[ fd ];
            F[ fd ]  = u[ fdd ];
            F[ fdd ] = - u[ f ] * u[ fdd ] - beta * ( 1.0 - u[ fd ] * u[ fd ] ); 
          }
      }; // End Falkner-Skan equation class

      class plate_BC : public Residual<double>
      {
        public:
          double KB;                              // Transpiration parameter

          plate_BC() : Residual<double> ( 2, 3 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + KB;
            B[ 1 ] = z[ fd ];
          }
      }; // End Falkner-Skan plate_BC class

      class far_BC : public Residual<double>
      {
        public:
          far_BC() : Residual<double> ( 1, 3 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ fd ] - 1.0;
          }
      }; // End Falkner-Skan far_BC class
#endif
                 
    } // End of namespace Base_Flow
   
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{ 
  cout << "*** ---------- Falkner-Skan with blowing ---------- ***" << endl;
  cout << "  * We are solving using " << Param::M + 1 << " points in the mesh with " 
       << "eta_inf = " << Param::eta_top << "." << endl;
  /* ----- Make the output directory ----- */
  std::ostringstream ss;
  ss << "./DATA/Falkner_Skan_with_blowing" << "_beta_" << Param::beta << "/";
  Example::output_path = ss.str();               
  int status = mkdir( Example::output_path.c_str(), S_IRWXU );
  cout << "  * Output directory " + Example::output_path + 
          " has been made successfully." << endl;
  cout << "  * Hartree parameter beta = " << Param::beta << endl;

  /* ----- Setup the mesh ----- */

  // define the remapped (non-uniform mesh) domain
  double bottom = Mesh::Y(0.0);
  double top    = Mesh::Y( Param::eta_top );

  // number of points to solve for 
  std::size_t N_eta = Param::M + 1;
  std::size_t N_Y( N_eta );

  // nodal positions in the remapped domain (spanned by X,Y)
  Vector<double> Y_nodes;
  Y_nodes.linspace( bottom, top, N_Y );

  // Vectors for original coordinates for writing data on the original zeta-eta domain
  Vector<double> eta_nodes;
  eta_nodes.linspace( 0.0, Param::eta_top, N_eta );

  // to find eta=eta(Y) and zeta=zeta(X) we will use Newton iteration 
  Mesh::invert_eta find_eta;
  Newton<double> newton_eta( &find_eta );
  for ( unsigned j = 0; j < N_Y; ++j )
  {
    unsigned kmin(0); double min(99);
    for ( unsigned k = 0; k < N_Y; ++k )
    {
      if ( std::abs( Mesh::Y( eta_nodes[k] ) - Y_nodes[j] ) < min )
      {
        min = std::abs( Mesh::Y( eta_nodes[k] ) - Y_nodes[j] );
        kmin = k;
      }
    }
    find_eta.Y0 = Y_nodes[ j ];
    Vector<double> guess( 1, 1.0 );
    guess[ 0 ] = eta_nodes[ kmin ];
    newton_eta.iterate( guess );
    eta_nodes[j] = guess[ 0 ];
  }  
  //

  // step sizes in the remapped domain : these should be constants
  const double dY( Y_nodes[ 1 ] - Y_nodes[ 0 ] );

  cout << "*** Solving the ODE ***" << endl;

  // Setup the base flow ODE problem
  Base_Flow::equation equation;
  Base_Flow::plate_BC plate_BC;
  Base_Flow::far_BC far_BC;
  equation.beta = 0.0;
  plate_BC.KB = 0.0;
  ODE_BVP<double> base( &equation, eta_nodes, &plate_BC, &far_BC );
  

  Vector<double> KB_vals;
  KB_vals.linspace( Param::KB, Param::KB_max, Param::KB_n );
  Vector<double> eta_half_vals;

  cout << "  *   KB    |   eta_half   |   U'(eta=0) " << endl;

do{   // Iterate over values of K

  // Set the initial guess
#ifdef BASE_2D
  for (std::size_t j=0; j < N_eta; ++j )
	{
		double eta = eta_nodes[ j ];				                      // eta value at node j
		base.solution()( j, f )  	= eta + exp( -eta );
    base.solution()( j, fd ) 	= 1.0 - exp( -eta ); 
		base.solution()( j, fdd ) = exp( -eta );
	}
#endif


  // Solve the system with KB = 0 then arc-length continue until KB = Param::KB
  double arc_step( 0.01 );
  double max_arc_step( 0.1 );
  base.init_arc( &plate_BC.KB, arc_step, max_arc_step ); 
  do
  {
    arc_step = base.arclength_solve( arc_step );
  }while( plate_BC.KB < Param::KB );
  plate_BC.KB = Param::KB;
  base.solve_bvp();                               // Solve once more with KB = Param::KB
  
  // Solve the system with beta = 0.1 then arc-length continue until beta = Param::beta
  arc_step = -0.01;
  if ( Param::beta >= 0.0 ) { arc_step = 0.01; }

  base.init_arc( &equation.beta, arc_step, max_arc_step );
  do
  {
    arc_step = base.arclength_solve( arc_step ); 
  }while( equation.beta < Param::beta );
  equation.beta = Param::beta;
  base.solve_bvp();

  // Store the solution in a mesh
  OneD_node_mesh<double> Base_soln( eta_nodes, 6 );
#ifdef BASE_2D
  for (std::size_t j=0; j < N_eta; ++j )
	{
		Base_soln( j, UB )      =   base.solution()( j, fd );
    Base_soln( j, UBd )     =   base.solution()( j, fdd );
    Base_soln( j, PhiB )    =   base.solution()( j, f );
    Base_soln( j, ThetaB )  =   ( 1.0 - Param::beta ) * base.solution()( j, fdd );
    Base_soln( j, ThetaBd ) =   ( 1.0 - Param::beta ) * ( - base.solution()( j, f ) *
                                base.solution()( j, fdd ) - Param::beta * ( 1.0 -   
                                base.solution()( j, fd ) * base.solution()( j, fd ) ) );
    Base_soln( j, PsiB )    =   ( 1.0 - Param::beta ) * base.solution()( j, fd );
	}
#endif

  // Output the solution to a file
  Base_soln.output( Example::output_path + "Base_soln_KB_" 
                  + Utility::stringify( abs( Param::KB ), 3 ) + ".dat" ); 
  // Output data to the screen

  cout << fixed << setprecision(4) << "  * " << plate_BC.KB << ",\t";


  // Find value of eta on zeta=0 at which U=1/2
  double eta_half;
  std::size_t lower = 0;
  std::size_t upper = 1;
  for (std::size_t j=0; j < Param::M; ++j)
  {
    if ( Base_soln(j,UB) < 0.5 && Base_soln(j+1,UB) > 0.5 ) { lower = j; upper=j+1; } 
  }
  // linearly interpolate
  eta_half =  ( 0.5 - Base_soln(lower,UB) ) * ( eta_nodes[upper] - eta_nodes[lower] ) 
            / ( Base_soln(upper,UB) - Base_soln(lower,UB)  ) + eta_nodes[lower];
  cout << eta_half << ", \t";
  eta_half_vals.push_back( eta_half );
  cout << Base_soln( 0, UBd ) << endl;

  Param::KB += KB_vals[ 1 ] - KB_vals[ 0 ];

}while( Param::KB <= Param::KB_max + 0.0001 );

  cout << "KB_vals.size() = " <<  KB_vals.size() << endl;
  cout << "eta_half_vals.size() = " <<  eta_half_vals.size() << endl;

  OneD_node_mesh<double> eta_half_KB( KB_vals, 1 );
  eta_half_KB.set_vars_from_vector( eta_half_vals );
  eta_half_KB.output( Example::output_path + "eta_half_KB.dat" );

  cout << "FINISHED" << endl;
}
