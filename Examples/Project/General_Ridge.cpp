#include <cassert>
#include <cmath>

#include "Core"

// Enumerations
enum{ f, fd, fdd, g, gd, gdd };                         // Base ODE enumeration
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };           // Base flow enumeration

// Either BASE_2D or BASE_3D for 2D or 3D base flows
#define BASE_2D
// Either UNIFORM or NONUNIFORM for uniform of non-uniform mesh
#define NONUNIFORM

namespace TSL
{
    namespace Param
    {
      double hzeta_right( 20.0 );       // Size of the domain in the zeta_hat direction
      double eta_top( 30.0 );           // Size of the domain in the eta direction
      const std::size_t N( 400 );       // Number of intervals in the zeta_hat direction
      const std::size_t M( 400 );       // Number of intervals in the eta direction
      const std::size_t Nvar( 4 );      // Number of variables
      double n( 0.0 );                  // Pressure gradient parameter
      double beta = (2.0*n)/(n+1.0);    // Hartree parameter
      double KB( 0.0 );                 // Base flow transpiration

    } // End of namespace Param

    namespace Mesh
    {
#ifdef UNIFORM    
      double X( const double& zeta )
      {
        return zeta;
      }
      double Xd( const double& zeta )
      {
        return 1;
      }
      double Xdd( const double& zeta )
      {
        return 0.0; 
      }

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
      const double a1( 0.1 );
      const double a2( 0.5 );   // X = (zeta + a1)^a2 

      double X( const double& zeta )
      {
        return std::pow(zeta + a1, a2);
      }
      double Xd( const double& zeta )
      {
        return a2 * std::pow(zeta + a1, a2 - 1);
      }
      double Xdd( const double& zeta )
      {
        return a2 * (a2 - 1) * std::pow(zeta + a1, a2 - 2);
      }

      const double b1( 0.1 );
      const double b2( 0.5 );   // Y = (eta_hat + b1)^b2

      double Y( const double& eta )
      {
        return std::pow(eta + b1, b2);
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

      class invert_zeta : public Residual<double>
      {
        // USED TO INVERT THE NON-UNIFORM MESH MAPPING
        public:
        double X0;

        invert_zeta() : Residual<double>( 1 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &f ) const
        {
          f[ 0 ] = X( z[0] ) - X0;
        }
      };
 
    } // End of namespace Mesh

    namespace Base_Flow
    {
      double K = Param::KB;             // Base transpiration parameter (+ve = blowing)   
      double beta = Param::beta;        // Hartree parameter

#ifdef BASE_2D
      class equation : public Equation<double>
      {
        public:
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
          plate_BC() : Residual<double> ( 2, 3 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + K;
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
#ifdef BASE_3D
      class equation : public Equation<double>
      {
        public:
          // The 3D alternative equation is 6th order
          equation() : Equation<double> ( 6 ) {}
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
                        - beta * ( 1.0 - u[ gd ] * u[ gd ] ) 
                        - 2.0 * ( 1.0 - beta ) * ( u[ fd ] - u[ gd ] ) * u[ gd ]; 
          }
      }; // End 3D alternative equation class

      class plate_BC : public Residual<double>
      {
        public:
          plate_BC() : Residual<double> ( 4, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + K;
            B[ 1 ] = z[ fd ];
            B[ 2 ] = z[ g ];
            B[ 3 ] = z[ gd ];
          }
      }; // End 3D alternative plate_BC class

      class far_BC : public Residual<double>
      {
        public:
          far_BC() : Residual<double> ( 2, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ fd ] - 1.0;
            B[ 1 ] = z[ gd ];
          }
      }; // End 3D alternative far_BC class
#endif                  
    } // End of namespace Base_Flow
   
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{ 
  cout << "*** ---------- General Ridge Code ---------- ***" << endl;
  //TODO output information about the mesh and number of ODE points etc
  

  /* ----- Setup the mesh ----- */

  // define the remapped (non-uniform mesh) domain
  double left   = Mesh::X( 0.0 );
  double right  = Mesh::X( Param::hzeta_right );
  double bottom = Mesh::Y(0.0);
  double top    = Mesh::Y( Param::eta_top );

  // number of points to solve for 
  std::size_t N_hzeta = Param::N + 1;
  std::size_t N_X( N_hzeta );
  std::size_t N_eta = Param::M + 1;
  std::size_t N_Y( N_eta );

  // nodal positions in the remapped domain (spanned by X,Y)
  Vector<double> X_nodes, Y_nodes;
  X_nodes.linspace( left, right, N_X );
  Y_nodes.linspace( bottom, top, N_Y );

  // Vectors for original coordinates for writing data on the original zeta-eta domain
  Vector<double> eta_nodes, hzeta_nodes;
  eta_nodes.linspace( 0.0, Param::eta_top, N_eta );
  hzeta_nodes.linspace( 0.0, Param::hzeta_right, N_hzeta );

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
  Mesh::invert_zeta find_zeta;
  Newton<double> newton_zeta( &find_zeta );
  //
  for ( unsigned i = 0; i < N_X; ++i )
  {
    unsigned kmin(0); double min(99);
    for ( unsigned k = 0; k < N_X; ++k )
    {
      if ( std::abs( Mesh::X( hzeta_nodes[k] ) - X_nodes[i] ) < min )
      {
        min = std::abs( Mesh::X( hzeta_nodes[k] ) - X_nodes[i] );
        kmin = k;
      }
    }
    find_zeta.X0 = X_nodes[ i ];
    Vector<double> guess( 1, 1.0 );
    guess[ 0 ] = hzeta_nodes[ kmin ];
    newton_zeta.iterate( guess );
    hzeta_nodes[ i ] = guess[ 0 ];
  } 

  // step sizes in the remapped domain : these should be constants
  const double dY( Y_nodes[ 1 ] - Y_nodes[ 0 ] );
  const double dX( X_nodes[ 1 ] - X_nodes[ 0 ] );

  /* ----- Solve the base flow ODE ----- */

  // Setup the base flow ODE problem
  Base_Flow::equation equation;
  Base_Flow::plate_BC plate_BC;
  Base_Flow::far_BC far_BC;
  ODE_BVP<double> base( &equation, eta_nodes, &plate_BC, &far_BC );
  // Set the initial guess
#ifdef BASE_2D
  for (std::size_t j=0; j < Param::M; ++j )
	{
		double eta = eta_nodes[ j ];				                      // eta value at node j
		base.solution()( j, f )  	= eta + exp( -eta );
    base.solution()( j, fd ) 	= 1.0 - exp( -eta ); 
		base.solution()( j, fdd ) = exp( -eta );
	}
#endif
#ifdef BASE_3D
  for (std::size_t j=0; j < Param::M; ++j )
	{
		double eta = eta_nodes[ j ];					                   // eta value at node j
		base.solution()( j, f )  	= eta + exp( -eta );
    base.solution()( j, fd ) 	= 1.0 - exp( -eta ); 
		base.solution()( j, fdd )  = exp( -eta );
    base.solution()( j, g )  	= 0.35 * (1.0 - exp( -eta ));
    base.solution()( j, gd ) 	= 1 - exp( -eta ) - exp( -1 / (eta * eta) ); 
		base.solution()( j, gdd )  = exp( -eta ) - 0.5 * tanh( eta ) + 0.5 * tanh( eta - 2.0 );
	}
#endif

  // Solve the system with KB = 0 and beta = 0 to provide a good initial guess
  Base_Flow::K = 0.0;
  Base_Flow::beta = 0.0;
  base.solve_bvp();

  // Solve the ODE system for the required parameters
  Base_Flow::K = Param::KB;
  Base_Flow::beta = Param::beta;
  base.solve_bvp();

  // Store the solution in a mesh
  OneD_node_mesh<double> Base_soln( eta_nodes, 6 );
#ifdef BASE_2D
  for (std::size_t j=0; j < Param::M; ++j )
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
#ifdef BASE_3D
  for (std::size_t j=0; j < Param::M; ++j )
	{
		Base_soln( j, UB )      =   base.solution()( j, fd );
    Base_soln( j, UBd )     =   base.solution()( j, fdd );
    Base_soln( j, PhiB )    =   base.solution()( j, f ) 
                              + ( 2.0 - Param::beta ) * base.solution()( j, g );
    Base_soln( j, ThetaB )  =   ( 1.0 - Param::beta ) * base.solution()( j, fdd )
                              - ( 2.0 - Param::beta ) * base.solution()( j, gdd );
    Base_soln( j, ThetaBd ) =   ( 1.0 - Param::beta ) * ( -(base.solution()( j, f ) + 
                                (2.0 - Param::beta) * base.solution()( j, g )) * 
                                base.solution()( j, fdd ) - Param::beta * ( 1.0 - 
                                base.solution()( j, fd ) * base.solution()( j, fd ) ) )
                              - ( 2.0 - Param::beta ) * ( -(base.solution()( j, f ) + 
                                (2.0 - Param::beta) * base.solution()( j, g )) * 
                                base.solution()( j, gdd ) - Param::beta * ( 1.0 - 
                                base.solution()( j, gd ) * base.solution()( j, gd ) ) - 
                                2.0 * (1.0 - Param::beta ) * (base.solution()( j, fd ) - 
                                base.solution()( j, gd )) * base.solution()( j, gd ) );
    Base_soln( j, PsiB )    =   ( 1.0 - Param::beta ) * base.solution()( j, fd )
                              - ( 2.0 - Param::beta ) * base.solution()( j, gd );
	}
#endif
  Base_soln.output( "./DATA/Base_soln.dat" );         // Output the solution to a file
  // Output the wall shear to the screen
#ifdef BASE_2D
  cout << "Base flow: 2D Falkner-Skan with transpiration" << endl; 
#endif
#ifdef BASE_3D
  cout << "Base flow: 3D alternative with transpiration" << endl;
#endif
  cout << "Base transpiration KB = " << Param::KB << endl;
  cout << "Hartree parameter beta = " << Param::beta << endl;
  cout << "U'(eta=0) =" << base.solution()( 0, fdd ) << endl;


  cout << "FINISHED" << endl;
}
