#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <sstream>

#include "Core"
#include "Sparse"

// Enumerations
enum{ f, fd, fdd, g, gd, gdd };                               // Base ODE
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };                 // Base ODE 
enum{ Ubar, Ubard, Phibar, Thetabar, Thetabard, Psibar };     // Far-field ODE           
enum{ Phi, Psi, U, Theta };                                   // PDE 

// Either BASE_2D or BASE_3D for 2D or 3D base flows
#define BASE_2D
// Either UNIFORM or NONUNIFORM for uniform of non-uniform mesh
#define NONUNIFORM
// Either DIRICHLET or NEUMANN boundary conditions on Phi and Psi at eta=eta_inf
#define DIRICHLET

namespace TSL
{
    namespace Param
    {
      double hzeta_right( 16.0 );       // Size of the domain in the zeta_hat direction
      double eta_top( 128.0 );          // Size of the domain in the eta direction
      const std::size_t N( 400 );       // Number of intervals in the zeta_hat direction
      const std::size_t M( 400 );       // Number of intervals in the eta direction
      const std::size_t Nvar( 4 );      // Number of variables
      double beta( 0.1 );               // Hartree parameter
      double KB( 0.0 );                 // Base flow transpiration ( +ve = blowing )
      double zeta0( 1.0 );              // Ridge/transpiration width
      double zeta0_2 = zeta0 * zeta0;   // Square of the ridge/transpiration width
      double A( 0.0 );                  // Mass flux parameter
      double K( 2.0 );                  // Transpiration parameter ( +ve = blowing )
      double gamma( 20.0 );             // Steepness factor

    } // End of namespace Param

    namespace Example
    {
      std::string output_path;          // Output path

      std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
      {
        // Return the column number for the kth variable at node (i,j)
        return Param::Nvar * ( i * ( Param::M + 1 ) + j ) + k;
      }

      double Phi_w( const double& hzeta )
      {
        // Return the transpiration function
        return -Param::K * 0.5 * ( 1. - tanh( Param::gamma * ( hzeta - 1. ) ) );
      }

      double Phi_w_hzeta( const double& hzeta )
      {
        // Return derivative of transpiration wrt hzeta
        double sech_squared = pow( cosh( Param::gamma * ( hzeta - 1 ) ) , -2. ); 
        return Param::K * 0.5 * Param::gamma * sech_squared;
      }                                       

    } // End of namespace Example


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

/*
    const double g0( 2 );// bigger => more points around zeta_hat=1
    const double B0( 4 ); // bigger => more points near zeta=0 and (less for large zeta)
    //
    double X( const double& zeta )
    {
        return B0 * zeta * 0.5 * ( 1+tanh( g0 * (1-zeta) ) ) 
              + (zeta + 2 * B0 ) * 0.5 * ( 1+tanh( g0 * (zeta-1) ) );
    }
    double Xd( const double& zeta )
    {
        return 0.5 * B0 * ( tanh( g0 * (1-zeta) ) + 1 ) 
             - 0.5 * B0 * g0 * zeta * std::pow( cosh( g0 * (1-zeta) ), -2 )
             + 0.5 * ( tanh( g0 * (zeta-1) ) + 1 ) 
             + 0.5 * g0 * (zeta + 2 * B0 ) * std::pow( cosh( g0 * (1-zeta) ), -2 );
    }
    double Xdd( const double& zeta )
    {
      return - B0 * g0 * std::pow( cosh( g0 * (1-zeta) ), -2 ) 
             + g0 * std::pow( cosh( g0 * (zeta-1) ),-2) 
             - B0 * g0 * g0 * zeta * std::pow(cosh(g0*(1-zeta)),-2)*tanh(g0*(1 - zeta)) 
             - g0 * g0 * (2*B0 + zeta)*std::pow(cosh(g0*(zeta-1)),-2)*tanh(g0*(zeta-1));
    }
*/

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
#ifdef BASE_3D
      class equation : public Equation<double>
      {
        public:
          double beta;                     // Hartree parameter      
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
                        -( 2.0 * ( 1.0 - beta ) * u[ fd ] 
                        - ( 2.0 - beta) * u[ gd ] ) * u[ gd ];
          }
      }; // End 3D alternative equation class

      class plate_BC : public Residual<double>
      {
        public:
          double KB;                        // Transpiration parameter

          plate_BC() : Residual<double> ( 4, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + KB;
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

    namespace Far_ODE
    {
      class far_equation : public Equation<double>
      {
        public:
          double beta;
          
          OneD_node_mesh<double> Base_soln;
          // The far-field equation is 6th order
          far_equation() : Equation<double> ( 6 ) {}
          // Define the equation
          void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
          { 
            double eta( coord( 0 ) );
            Vector<double> Base( Base_soln.get_interpolated_vars( eta ) );

            F[ Phibar ]     =   ( 2. - beta ) * u[ Ubar ] + u[ Psibar ];
            F[ Psibar ]     =   u[ Thetabar ];
            F[ Ubar ]       =   u[ Ubard ];
            F[ Ubard ]      =   2. * beta * Base[ UB ] * u[ Ubar ] 
                              - Base[ UBd ] * u[ Phibar ] - Base[ PhiB ] * u[ Ubard ]
                              + 2. * Base[ PsiB ] * u[ Ubar ];
            F[ Thetabar ]   =   u[ Thetabard ];
            F[ Thetabard ]  =   2. * ( 1. - beta ) * ( Base[ UB ] * u[ Ubard ] 
                              + u[ Ubar ] * Base[ UBd ] ) 
                              - Base[ PhiB ] * u[ Thetabard ] 
                              - Base[ ThetaBd ] * u[ Phibar ] 
                              + Base[ PsiB ] * u[ Thetabar ]
                              - u[ Psibar ] * Base[ ThetaB ]
                              - ( 2. - beta ) * ( Base[ UB ] * u[ Thetabar ] 
                              + Base[ ThetaB ] * u[ Ubar ] ) ;
          }
      }; // End of far-field ODE class

      class far_plate_BC : public Residual<double>
      {
        public:

          far_plate_BC() : Residual<double> ( 3, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ Ubar ];
            B[ 1 ] = z[ Phibar ];
            B[ 2 ] = z[ Psibar ];
          }
      }; // End of far-field plate BC class

      class far_far_BC : public Residual<double>
      {
        public:

          far_far_BC() : Residual<double> ( 3, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ Ubar ];
            B[ 1 ] = z[ Thetabar ];
            B[ 2 ] = z[ Psibar ] - 1.0;
          }
      }; // End of far-field far BC class


    }
   
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{ 
  cout << "*** ---------- Blowing with pressure gradient ---------- ***" << endl;
  cout << "  * We are solving using a " << Param::N + 1 << " x " << Param::M + 1 
       << " mesh with zeta_hat_inf = " << Param::hzeta_right << " and eta_inf = " 
       << Param::eta_top << "." << endl;

  /* ----- Make the output directory ----- */
  std::ostringstream ss;
  ss << "./DATA/K_" << Param::K << "_" << "beta_" << Param::beta << "_" << Param::N + 1 
     << "x" << Param::M + 1 << "_" << Param::hzeta_right << "_" << Param::eta_top << "/";
  Example::output_path = ss.str();               
  int status = mkdir( Example::output_path.c_str(), S_IRWXU );
  cout << "  * Output directory " + Example::output_path + 
          " has been made successfully." << endl;

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

  cout << "*** Solving the base flow ODE ***" << endl;

  // Setup the base flow ODE problem
  Base_Flow::equation equation;
  Base_Flow::plate_BC plate_BC;
  Base_Flow::far_BC far_BC;
  equation.beta = 0.1;
  plate_BC.KB = 0.0;
  ODE_BVP<double> base( &equation, eta_nodes, &plate_BC, &far_BC );

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
#ifdef BASE_3D
  for (std::size_t j=0; j < N_eta; ++j )
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
  if ( Param::beta >= 0.1 ) { arc_step = 0.01; }

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
#ifdef BASE_3D
  for (std::size_t j=0; j < N_eta; ++j )
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
  // Output the solution to a file
  Base_soln.output( Example::output_path + "Base_soln.dat" ); 
  // Output the wall shear to the screen
#ifdef BASE_2D
  cout << "  * Base flow: 2D Falkner-Skan with transpiration" << endl; 
  
#endif
#ifdef BASE_3D
  cout << "  * Base flow: 3D alternative with transpiration" << endl;
#endif
  cout << "  * This number should be zero for the 2D ODE solution and non-zero for the " 
       << " 3D solution: " << ( 1. - Param::beta ) * Base_soln.integral2(UB) 
                                                   - Base_soln.integral2(PsiB) << endl;
  cout << "  * Base transpiration KB = " << plate_BC.KB << endl;
  cout << "  * Hartree parameter beta = " << equation.beta << endl;
  cout << "  * UB'(eta=0) =" << base.solution()( 0, fdd ) << endl;
  cout << "  * We have solved the ODE problem, it is output to " + Example::output_path +
          "Base_soln.dat" << endl;

  /* ----- Solve the far-field ODE ----- */

  // Setup the far-field ODE problem
  Far_ODE::far_equation far_equation;
  Far_ODE::far_plate_BC far_plate_BC;
  Far_ODE::far_far_BC   far_far_BC; 
  far_equation.beta = Param::beta;
  far_equation.Base_soln = Base_soln;
   
  // Create boundary value problem
  ODE_BVP<double> far_ode( &far_equation, eta_nodes, &far_plate_BC, &far_far_BC );

  // Solve the system
  far_ode.solve_bvp();
  //cout << "Thetabar(0) = " << far_ode.solution()( 0, Thetabar ) << endl;

  /* ----- Solve for the perturbation quantities ----- */

  cout << "*** Solving the perturbation equations ***" << endl;
  cout << "  * Perturbation transpiration K = " << Param::K << endl;
  // Set the current guess states  
  TwoD_node_mesh<double> Q( X_nodes, Y_nodes, 4 );
  // We use the mesh below to write the data on the original zeta-eta domain
  TwoD_node_mesh<double> Q_output( hzeta_nodes, eta_nodes, 8 );

  // output a measure of the solution
  TrackerFile metric( Example::output_path + "A_file.dat" );  
  metric.push_ptr( &Param::zeta0, "zeta0" );
  metric.push_ptr( &Param::A, "A" );
  double U_eta( 0.0 );
  double eta_half( 0.0 );
  metric.push_ptr( &U_eta, "U_eta(0,0)");
  metric.push_ptr( &eta_half, "eta at which U=1/2 on zeta=0" );
  metric.header();
  

  // Vector for the RHS of the matrix problem
  Vector<double> B( 4 * N_eta * N_hzeta + 1, 0.0 );

  do                                                    // Iterate over values of zeta_0
  {

  /* Iterate to a solution */
  double max_residual( 0.0 );                           // Maximum residual
  std::size_t max_iterations( 20 );                     // Maximum number of iterations
  std::size_t iteration( 0 );                           // Initialise iteration counter

  do
  { 
    // N_eta x N_hzeta mesh with 4 unknowns at each node + 1 for mass flux parameter A
    SparseMatrix<double> A( 4 * N_eta * N_hzeta + 1, 4 * N_eta * N_hzeta + 1 ); 
    //Sparse_matrix<double> A( 4 * N_eta * N_hzeta + 1, 4 * N_eta * N_hzeta + 1 ); 
    cout << "  * Assembling sparse matrix problem" << endl; 

    Timer timer;
    timer.start();

    using namespace Example;
    std::size_t row( 0 );                               // Initialise row counter

    /* hzeta = 0 boundary ( left boundary ) */
      std::size_t i( 0 );
      
      for ( std::size_t j = 0; j < Param::M + 1 ; ++j )
      {
          double hzeta( hzeta_nodes[ 0 ] );
          double Xd( Mesh::Xd(hzeta) ); 
          double eta( eta_nodes[ j ] );
          double Yd( Mesh::Yd( eta ) );
          Vector<double> Base( Base_soln.get_interpolated_vars( eta ) );
          // PhiB' = (2-beta)*UB - PsiB
          double PhiBd( ( 2.0 - Param::beta ) * Base[ UB ] - Base[ PsiB ] );
          //double UBd( Base[ UBd ] );      // TODO do we need this    
          

          // Phi_hzeta = 0
          A( row, col( i, j, Phi ) )      = -3.*Xd/(2*dX); 
          A( row, col( i + 1, j, Phi ) )  =  4.*Xd/(2*dX);
          A( row, col( i + 2, j, Phi ) )  = -1.*Xd/(2*dX);

          B[ row ]                        = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
                                            -Q(i+2,j,Phi) )/(2*dX) );
          ++row;

          // Psi = 0
          A( row, col( i, j, Psi ) )      =   1;
          B[ row ]                        = - ( Q( i, j, Psi ) );
          ++row;

          // U_hzeta = 0
          A( row, col( i, j, U ) )        = -3.*Xd/(2*dX); 
          A( row, col( i + 1, j, U ) )    =  4.*Xd/(2*dX);
          A( row, col( i + 2, j, U ) )    = -1.*Xd/(2*dX);

          B[ row ]                        = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U)  
                                             -Q(i+2,j,U) )/(2*dX) );
          ++row;

          // Theta = 0
          A( row, col( i, j, Theta ) )    =   1;
          B[ row ]                        = - Q( i, j, Theta );
          ++row;

      } // end for loop over LHS eta nodes

    /* Interior points between the hzeta boundaries */
    for ( std::size_t i = 1; i < Param::N; ++i )
    {
      // hzeta location
      double hzeta( hzeta_nodes[ i ] );
      double Xd( Mesh::Xd( hzeta ) );
      double Xdd( Mesh::Xdd( hzeta ) );
      // Wall transpiration
      double Phi_w( Example::Phi_w( hzeta ) );
      double Phi_w_hzeta( Example::Phi_w_hzeta( hzeta ) );

      /* eta = 0 boundary ( bottom boundary ) */
      std::size_t j( 0 );						
      double eta( eta_nodes[ j ] );
      double Yd( Mesh::Yd(eta) );
        	
      // Phi = Phi_w
      A( row, col( i, j, Phi ) )        =  1.;
      B[ row ]                          = -Q( i, j, Phi ) + Phi_w;
      ++row;
      // Psi = 0
      A( row, col( i, j, Psi ) )        =  1.;
      B[ row ]                          = -Q( i, j, Psi );
      ++row;
      // U = 0
      A( row, col( i, j, U ) )          =  1.;
      B[ row ]                          = -Q( i, j, U );
      ++row;
      // Theta - Psi_eta = -( 1 / ( zeta0^2 ) ) * Phi_w_hzeta
      A( row, col( i, j, Theta ) )      =  1.;
      A( row, col( i, j, Psi ) )        =  3.*Yd / ( 2 * dY );
      A( row, col( i, j + 1, Psi ) )    = -4.*Yd / ( 2 * dY );
      A( row, col( i, j + 2, Psi ) )    =  1.*Yd / ( 2 * dY );
      B[ row ]                          = -Q(i,j,Theta) + Yd*( -3*Q(i,j,Psi) 
                                          + 4*Q(i,j+1,Psi) - Q(i,j+2,Psi) ) / (2*dY)
                                          - ( 1. / ( Param::zeta0_2 ) ) 
                                          * Phi_w_hzeta;
      ++row; 


      /* Main interior grid points */
      for ( std::size_t j = 1; j < Param::M; ++j )
      {
        // eta location
        double eta( eta_nodes[ j ] );
        double Yd( Mesh::Yd( eta ) );
        double Ydd( Mesh::Ydd( eta ) );
        // Base solution
        Vector<double> Base( Base_soln.get_interpolated_vars( eta ) );
        // PhiB' = (2-beta)*UB - PsiB
        double PhiBd( ( 2.0 - Param::beta )*Base[ UB ] - Base[ PsiB ] );
        // PhiB'' = (2-beta)*UB' - PsiB' = (2-beta)*UB' - ThetaB
        double PhiBdd( ( 2.0 - Param::beta )*Base[ UBd ] -  Base[ ThetaB ] );
        // PsiB' = Theta_B
        double PsiBd( Base[ ThetaB ] );
        // PsiB'' = ThetaB'
        double PsiBdd( Base[ ThetaBd ] ); 
        // UB'' = beta * [ UB^2 - 1] - PhiB * UB'
        double UBdd(  Param::beta * ( Base[ UB ] * Base[ UB ]  - 1. ) 
                    - Base[ PhiB ] * Base[ UBd ] );
        // ThetaB'' = 2(1-beta)*UB*UB' - PhiB*ThetaB' - PsiB*ThetaB - (2-beta)*UB*ThetaB
        double ThetaBdd( 2. * ( 1. - Param::beta ) * Base[ UB ] * Base[ UBd ]
                        - Base[ PhiB ] * Base[ ThetaBd ] - Base[ PsiB ] * Base[ ThetaB ]
                        - ( 2. - Param::beta ) * Base[ UB ] * Base[ ThetaB ] );
        
        // Laplacian coefficients for finite-differencing
        
        // X(i,j-1)
        double laplace_1 =  ( Yd*Yd/(dY*dY) - Ydd/ (2.*dY) ) ;
        // X(i-1,j)
        double laplace_3 = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) ) / ( Param::zeta0_2 );
        // X(i,j)
        double laplace_4 = -2.*( Yd*Yd / (dY*dY) 
                           + Xd*Xd/( Param::zeta0_2 * dX * dX ) );
        // X(i+1,j)
        double laplace_5 = ( Xdd/(2.*dX) + Xd*Xd/(dX*dX) ) / ( Param::zeta0_2 );
        // X(i,j+1)
        double laplace_7 = ( Yd*Yd/(dY*dY) + Ydd/ (2.*dY) );

        // Guessed/known components and various derivative values
        Vector<double> Guess( Q.get_nodes_vars( i, j ) );
        Vector<double> Guess_eta( ( Q.get_nodes_vars( i, j + 1 ) 
                                  - Q.get_nodes_vars( i, j - 1 ) ) * ( Yd /( 2 * dY )) );
        Vector<double> Guess_hzeta( ( Q.get_nodes_vars( i + 1, j ) 
                                    - Q.get_nodes_vars( i - 1, j ) ) 
                                    * ( Xd /( 2 * dX )) );
        Vector<double> Guess_laplace( Q.get_nodes_vars( i, j - 1 ) * laplace_1 
                                   +  Q.get_nodes_vars( i - 1, j ) * laplace_3
                                   +  Q.get_nodes_vars( i, j ) * laplace_4
                                   +  Q.get_nodes_vars( i + 1, j ) * laplace_5
                                   +  Q.get_nodes_vars( i, j + 1 ) * laplace_7 );
        
        //////////////////
        // Phi equation //
        //////////////////

        // Laplacian of Phi        
        A( row, col( i, j - 1, Phi ) )      = laplace_1;
        A( row, col( i - 1, j, Phi ) )      = laplace_3;
        A( row, col( i, j, Phi ) )          = laplace_4;
        A( row, col( i + 1, j, Phi ) )      = laplace_5;
        A( row, col( i, j + 1, Phi ) )      = laplace_7;
        // -(2-beta)*U_eta
        A( row, col( i, j + 1, U ) )        = -( 2. - Param::beta )*Yd/( 2 * dY );
        A( row, col( i, j - 1, U ) )        =  ( 2. - Param::beta )*Yd/( 2 * dY );
        // Theta_hzeta
        A( row, col( i + 1, j, Theta ) )    =  Xd / ( 2 * dX );
        A( row, col( i - 1, j, Theta ) )    = -Xd / ( 2 * dX );

        // Residual
        B[ row ]      = - Guess_laplace[ Phi ] + ( 2. - Param::beta ) * Guess_eta[ U ]
                        - Guess_hzeta[ Theta ]; 
        ++row;

        //////////////////
        // Psi equation //
        //////////////////

        // Laplacian of Psi
        A( row, col( i, j - 1, Psi ) )      = laplace_1;
        A( row, col( i - 1, j, Psi ) )      = laplace_3;
        A( row, col( i, j, Psi ) )          = laplace_4;
        A( row, col( i + 1, j, Psi ) )      = laplace_5;
        A( row, col( i, j + 1, Psi ) )      = laplace_7;

        // -(2-beta)*U_hzeta / (zeta0^2)
        A( row, col( i + 1, j, U ) )        = - ( 2. - Param::beta ) * Xd 
                                              / ( 2. * dX * Param::zeta0_2 );
        A( row, col( i - 1, j, U ) )        =   ( 2. - Param::beta ) * Xd
                                              / ( 2. * dX * Param::zeta0_2 );

        // -Theta_eta
        A( row, col( i, j + 1, Theta ) )    = - Yd / ( 2 * dY ); 
        A( row, col( i, j - 1, Theta ) )    =   Yd / ( 2 * dY ); 

        // Residual
        B[ row ]      = - Guess_laplace[ Psi ] + ( 2. - Param::beta ) 
                        * ( Guess_hzeta[ U ] )
                        / ( Param::zeta0_2 )
                        + Guess_eta[ Theta ];

        ++row;

        ////////////////
        // U equation //
        ////////////////

        // Laplacian of U
        A( row, col( i, j - 1, U ) )        = laplace_1;
        A( row, col( i - 1, j, U ) )        = laplace_3;
        A( row, col( i, j, U ) )            = laplace_4;
        A( row, col( i + 1, j, U ) )        = laplace_5;
        A( row, col( i, j + 1, U ) )        = laplace_7;

        // -2 * beta * ( UB + UG ) * U
        A( row, col( i, j, U ) )           += - 2.* Param::beta * ( Base[ UB ] 
                                              + Guess[ U ] );

        // ( hzeta * PsiB + PsiG ) * U_hzeta
        A( row, col( i + 1, j, U ) )       +=   ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                              * Xd / ( 2 * dX ); 
        A( row, col( i - 1, j, U ) )       += - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                              * Xd / ( 2 * dX );

        // [ PhiB + PhiG ] * U_eta 
        A( row, col( i, j + 1, U ) )       +=   ( Base[ PhiB ] + Guess[ Phi ] )
                                              * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )       += - ( Base[ PhiB ] + Guess[ Phi ] )
                                              * Yd / ( 2 * dY );
          
        // [ UG_hzeta ] * Psi
        A( row, col( i, j, Psi ) )          =   Guess_hzeta[ U ];

        // ( UB' + UG_eta ) * Phi
        A( row, col( i, j, Phi ) )          =    Base[ UBd ] + Guess_eta[ U ];
         
        // Residual
        B[ row ]        = - Guess_laplace[ U ] 
                          + Param::beta * ( 2. * Base[ UB ] + Guess[ U ] ) * Guess[ U ]
                          - ( hzeta * Base[ PsiB ] + Guess[ Psi ] ) 
                          * ( Guess_hzeta[ U ] ) 
                          - (Base[PhiB] + Guess[Phi]) * Guess_eta[ U ]
                          - Base[UBd] * Guess[Phi] ;
        ++row;

        ////////////////////
        // Theta equation //
        ////////////////////

        // Laplacian of Theta
        A( row, col( i, j - 1, Theta ) )     = laplace_1;
        A( row, col( i - 1, j, Theta ) )     = laplace_3;
        A( row, col( i, j, Theta ) )         = laplace_4;
        A( row, col( i + 1, j, Theta ) )     = laplace_5;
        A( row, col( i, j + 1, Theta ) )     = laplace_7;
        // -2 * (1-beta) * (UB+UG) * [hzeta] * U_eta
        A( row, col( i, j + 1, U ) )         = - 2. * ( 1. - Param::beta )
                                                 * ( Base[ UB ] + Guess[ U ] ) 
                                                 * ( hzeta ) 
                                                 * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )         =   2. * ( 1. - Param::beta )
                                                 * ( Base[ UB ] + Guess[ U ] ) 
                                                 * ( hzeta ) 
                                                 * Yd / ( 2 * dY );

        // -2 * (1-beta) * (UB + UG) * ( hzeta ) * U
        A( row, col( i, j, U ) )             = - 2. * ( 1. - Param::beta )
                                                 * ( Base[ UBd ] + Guess_eta[ U ] ) 
                                                 * ( hzeta );

        // (2 * (1-beta) * (eta + H) * UG_hzeta / (zeta0^2)) * U
        A( row, col( i, j, U ) )            +=  2. * ( 1. - Param::beta )
                                                * eta * Guess_hzeta[ U ]
                                                / ( Param::zeta0_2 );

        // 2 * (1-beta) * eta * (UB + UG) * U_hzeta / ( zeta0^2 )
        A( row, col( i + 1, j, U ) )         =  2. * ( 1. - Param::beta ) 
                                                * eta 
                                                * ( Base[ UB ] + Guess[ U ] )
                                                * Xd / ( 2 * dX * Param::zeta0_2 );
        A( row, col( i - 1, j, U ) )         = -2. * ( 1. - Param::beta ) 
                                                * eta 
                                                * ( Base[ UB ] + Guess[ U ] )
                                                * Xd / ( 2 * dX * Param::zeta0_2 );

        // ( PhiB + PhiG ) * Theta_eta
        A( row, col( i, j + 1, Theta ) )    +=  ( Base[ PhiB ] + Guess[ Phi ] ) * Yd 
                                                / ( 2 * dY );
        A( row, col( i, j - 1, Theta ) )    += -( Base[ PhiB ] + Guess[ Phi ] ) * Yd 
                                                / ( 2 * dY );

        // (hzeta * ThetaB' + ThetaG_eta ) * Phi
        A( row, col( i, j, Phi ) )           =   hzeta * Base[ ThetaBd ] 
                                               + Guess_eta[ Theta ];

        // (hzeta * PsiB + PsiG ) * Theta_hzeta
        A( row, col( i + 1, j, Theta ) )    +=  ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                               * Xd / ( 2 * dX );
        A( row, col( i - 1, j, Theta ) )    += -( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                               * Xd / ( 2 * dX );

        // [ThetaB + ThetaG_hzeta] * Psi
        A( row, col( i, j, Psi ) )           =  Base[ ThetaB ] + Guess_hzeta[ Theta ];

        // (2-beta) * ( UB + UG ) * Theta
        A( row, col( i, j, Theta ) )        +=   ( 2. - Param::beta ) * ( Base[ UB ] 
                                               + Guess[ U ] );

        // (2-beta) * ( hzeta * ThetaB + ThetaG ) * U
        A( row, col( i, j, U ) )            +=   ( 2. - Param::beta ) 
                                               * ( hzeta * Base[ Theta ] 
                                               + Guess[ Theta ] );
          
        // Residual
        B[ row ]      = - Guess_laplace[ Theta ]
                        + 2.*( 1. - Param::beta ) 
                        * ( hzeta * ( Base[ UB ] + Guess[ U ] ) 
                        * Guess_eta[ U ] + hzeta * Base[ UBd ] * Guess[ U ] 
                        - eta * ( Base[ UB ] + Guess[ U ] ) 
                        * ( Guess_hzeta[ U ] ) 
                        / ( Param::zeta0_2 ) )  
                        - ( Base[ PhiB ] + Guess[ Phi ] ) * Guess_eta[ Theta ]
                        - hzeta * Base[ ThetaBd ] * Guess[ Phi ]
                        - ( hzeta * Base[ PsiB ] + Guess[ Psi ] ) 
                        * ( Guess_hzeta[ Theta ] ) - Guess[ Psi ] * Base[ ThetaB ]
                        - ( 2. - Param::beta ) * ( ( Base[ UB ] + Guess[ U ] ) 
                        * Guess[ Theta ] + hzeta * Base[ ThetaB ] * Guess[ U ] );
        ++row;


      }

      /* eta = eta_inf boundary ( top boundary ) */
        j = Param::M ;
        eta = eta_nodes[ j ];
        Yd = Mesh::Yd( eta );

#ifdef DIRICHLET
        // Phi = A*(...)
        A( row, col( i, j, Phi ) )        =   1.0;
        A( row, 4 * N_hzeta * N_eta )     = - eta / ( eta * eta 
                                            + Param::zeta0_2 * hzeta 
                                            * hzeta );

        B[ row ]        = - Q( i, j, Phi )
                          + eta * Param::A / ( eta * eta 
                          + Param::zeta0_2 * hzeta * hzeta );

        ++row;

        // Psi = A*(...)
        A( row, col( i, j, Psi ) )        =   1.0;
        A( row, 4 * N_hzeta * N_eta )     = - hzeta / ( eta * eta 
                                            + Param::zeta0_2 * hzeta 
                                            * hzeta );

        B[ row ]        = - Q( i, j, Psi ) + hzeta * Param::A 
                          / ( eta * eta 
                          + Param::zeta0_2 * hzeta * hzeta );

        ++row;
#endif
#ifdef NEUMANN
        // Phi_eta = A*(...) - derivative condition
        A( row, col( i, j, Phi ) )        =   3.0 *Yd/ (2*dY);
        A( row, col( i, j - 1, Phi ) )    = - 4.0 *Yd/ (2*dY);
        A( row, col( i, j - 2, Phi ) )    =   1.0 *Yd/ (2*dY);
        A( row, 4 * N_hzeta * N_eta )     = - ( Param::zeta0_2 
                                            * hzeta * hzeta 
                                            - eta * eta ) 
                                            / pow( ( Param::zeta0_2 
                                            * hzeta * hzeta 
                                            + eta * eta ), 2);
    
        B[ row ]        = - ( 3 * Q( i, j, Phi ) - 4 * Q( i, j-1, Phi ) 
                          + Q( i, j-2, Phi ) ) * Yd / ( 2 * dY ) 
                          + Param::A * ( Param::zeta0_2 * hzeta * hzeta 
                          - eta * eta ) 
                          / pow( ( Param::zeta0_2 * hzeta * hzeta 
                          + eta * eta ), 2);
        ++row;

        // Psi_eta = A*(...)
        A( row, col( i, j, Psi ) )        =   3.0 *Yd/ (2*dY);
        A( row, col( i, j - 1, Psi ) )    = - 4.0 *Yd/ (2*dY);
        A( row, col( i, j - 2, Psi ) )    =   1.0 *Yd/ (2*dY);
        A( row, 4 * N_hzeta * N_eta )     =   2. * hzeta * eta
                                            / pow( ( Param::zeta0_2 
                                            * hzeta * hzeta 
                                            + eta * eta ) , 2 );

        B[ row ]        = - ( 3 * Q( i, j, Psi ) - 4 * Q( i, j-1, Psi ) 
                          + Q( i, j-2, Psi ) ) * Yd / ( 2 * dY )  
                          - Param::A * 2. * hzeta * eta
                          / pow( ( Param::zeta0_2 * hzeta * hzeta 
                          + eta * eta ) , 2 );
        ++row;
#endif 
          
        // U = 0
        A( row, col( i, j, U ) )            =   1;
        B[ row ]                            = - ( Q( i, j, U ) );
        ++row;
        
        // Theta = 0
        A( row, col( i, j, Theta ) )        =   1;
        B[ row ]                            = - ( Q(i,j,Theta) );
        ++row;

    } // End of for loop over interior nodes

    /* hzeta = hzeta_inf boundary ( right boundary ) */

    for ( std::size_t j = 0; j < Param::M + 1; ++j )
    {
      //offset for global problem
      std::size_t i( N_hzeta-1 );
      double hzeta( hzeta_nodes[ i ] );
      double Xd( Mesh::Xd(hzeta) );
      double eta( eta_nodes[j] );
        
        
      Vector<double> Base( base.solution().get_interpolated_vars( eta ) );

      // hzeta * Phi_hzeta + 2 * Phi = A*(...)
      A( row, col( i, j, Phi ) )          =   hzeta * 3. * Xd / ( 2 * dX ) + 2.;
      A( row, col( i - 1, j, Phi ) )      = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, Phi ) )      =   hzeta * 1. * Xd / ( 2 * dX );
      A( row, 4 * N_hzeta * N_eta )       = - 2 * pow( eta, 3 ) 
                                            / pow( Param::zeta0_2 
                                            * hzeta * hzeta 
                                            + eta * eta, 2 ); 

      B[ row ]        = - hzeta * ( 3 * Q( i, j, Phi) - 4 * Q( i - 1, j, Phi) 
                        + Q( i - 2, j, Phi) ) * Xd / ( 2 * dX ) 
                        - 2 * Q( i, j, Phi )
                        + 2 * Param::A * pow( eta, 3. ) 
                        / pow( Param::zeta0_2 * hzeta * hzeta 
                        + eta * eta, 2 );
      ++row;

      // hzeta * Psi_hzeta + Psi = A*(...) 
      A( row, col( i, j, Psi ) )          =   hzeta * 3. * Xd / ( 2 * dX ) + 1.;
      A( row, col( i - 1, j, Psi ) )      = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, Psi ) )      =   hzeta * 1. * Xd / ( 2 * dX );
      A( row, 4 * N_hzeta * N_eta )       = - 2. * hzeta * pow( eta, 2. )
                                            / pow( Param::zeta0_2 
                                            * hzeta * hzeta + eta * eta, 2 );

      B[ row ]        = - hzeta * ( 3 * Q( i, j, Psi ) - 4 * Q( i - 1, j, Psi ) 
                        + Q( i - 2, j, Psi) ) * Xd / ( 2 * dX ) 
                        - Q( i, j, Psi)
                        + 2. * Param::A * hzeta * pow( eta, 2. )
                        / pow( Param::zeta0_2 * hzeta * hzeta 
                        + eta * eta, 2 ) ;
      ++row;

      // hzeta * U_hzeta + 2 * U = 0
      A( row, col( i, j, U ) )            =   hzeta * 3. * Xd / ( 2 * dX ) + 2.;
      A( row, col( i - 1, j, U ) )        = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, U ) )        =   hzeta * 1. * Xd / ( 2 * dX );

      B[ row  ]       = - hzeta * ( 3 * Q( i, j, U ) - 4 * Q( i - 1, j, U ) 
                        + Q( i - 2, j, U) ) * Xd / ( 2 * dX ) - 2 * Q( i, j, U );
        
      ++row;

      // hzeta * Theta_hzeta + Theta = 0
      A( row, col( i, j, Theta )  )       =   hzeta * 3. * Xd / ( 2 * dX ) + 1.;
      A( row, col( i - 1, j, Theta ) )    = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, Theta ) )    =   hzeta * 1. * Xd / ( 2 * dX );

      B[ row ]        = - hzeta * ( 3 * Q( i, j, Theta ) - 4 * Q( i - 1, j, Theta ) 
                        + Q( i - 2, j, Theta ) ) * Xd / ( 2 * dX ) - Q( i, j, Theta ) ;

      ++row;

    }

    double hzeta( hzeta_nodes[ Param::N ] );    // hzeta = hzeta_inf
    //double zeta0( Param::zeta0 );
    /* A coefficient condition */
    // zeta0^2 * hzeta_max theta( zeta=zeta_max, eta=0 ) = A*otheta(0)
    Vector<double> Base( base.solution().get_interpolated_vars( 0.0 ) );
    A( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta ) = - far_ode.solution()( 0, Thetabar );
    A( 4 * N_eta * N_hzeta, 4 * N_eta * (N_hzeta - 1) + 3 ) =   Param::zeta0_2 * hzeta;
    // RHS
    B[ row ] = - Q( N_hzeta-1, 0, Theta ) * Param::zeta0_2 * hzeta 
               + Param::A * far_ode.solution()( 0, Thetabar ) ;   
    
    max_residual = B.norm_inf();
    cout << "***                                              Maximum residual = " 
         << B.norm_inf() << endl;

    //Timer timer;
    //timer.start();
    Vector<double> x;
    x = A.solve( B );      
    B = x;   
    timer.print();
    timer.stop();

    // Update the known values using the correction which we just found
    
    for ( std::size_t i = 0; i < Param::N + 1; ++i )
    {
      for ( std::size_t j = 0; j < Param::M + 1; ++j )
      {
        Q( i, j, Phi )    += B[ col( i, j, Phi ) ];
        Q( i, j, Psi )    += B[ col( i, j, Psi ) ];
        Q( i, j, U )      += B[ col( i, j, U ) ];
        Q( i, j, Theta )  += B[ col( i, j, Theta ) ];
      }
    }
    Param::A += B[ 4 * ( Param::N + 1 ) * ( Param::M + 1 ) ];
  
    cout << "***    Iteration = " << iteration 
           << "    Maximum correction = " << B.norm_inf() << endl;

    ++iteration;
  }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) ); // End iteration
  
  if ( iteration >= max_iterations )
  {
    cout << "STOPPED AFTER TOO MANY ITERATIONS" << endl;
  }

  /* push the data back into the unmapped domain */
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      double hzeta=hzeta_nodes[i];
      for ( std::size_t j = 0; j < N_eta; ++j )
      {
        double eta=eta_nodes[j];
        // first 4 values output are the without the underlying 2D base flow
        Q_output( i, j, 0 ) = Q( i, j, Phi);
        Q_output( i, j, 1 ) = Q( i, j, Psi);
        Q_output( i, j, 2 ) = Q( i, j, U);
        Q_output( i, j, 3 ) = Q( i, j, Theta);
        // second 4 values are the "full" solution, but still with the zeta0 scaling
        Q_output( i, j, 4 ) =   Q( i, j, Phi)   
                            + Base_soln.get_interpolated_vars( eta )[PhiB];
        Q_output( i, j, 5 ) = Q( i, j, Psi)   
                            + hzeta * Base_soln.get_interpolated_vars( eta )[PsiB];
        Q_output( i, j, 6 ) =   Q( i, j, U)
                            + Base_soln.get_interpolated_vars( eta )[UB];
        Q_output( i, j, 7 ) = Q( i, j, Theta) 
                            + hzeta * Base_soln.get_interpolated_vars( eta )[ThetaB];
      }
    }
    // Convert to string //TODO we need a utility function to do this
    std::stringstream ss;
    ss << Param::zeta0;
    std::string zeta0_str = ss.str(); 

    Q_output.dump_gnu( Example::output_path + "Qout_" + zeta0_str + ".dat" );

    Vector<double> Base( Base_soln.get_interpolated_vars( 0.0 ) ); 
    U_eta = -( 3 * Q_output(0,0,U+4) - 4 * Q_output(0,1,U+4) 
            + Q_output(0,2,U+4) ) * Mesh::Yd(0.0)/(2*dY);
  
    // Find value of eta on zeta=0 at which U=1/2
    std::size_t lower = 0;
    std::size_t upper = 1;
    for (std::size_t j=0; j < Param::M; ++j)
    {
    if ( Q_output(0,j,U+4) < 0.5 && Q_output(0,j+1,U+4) > 0.5 ) { lower = j; upper=j+1; } 
    }
    // linearly interpolate
    eta_half =  ( 0.5 - Q_output(0,lower,U+4) ) * ( eta_nodes[upper] - eta_nodes[lower] ) 
              / ( Q_output(0,upper,U+4) - Q_output(0,lower,U+4)  ) + eta_nodes[lower];

    metric.update();

    // Get the wall shear values as a function of hzeta
    OneD_node_mesh<double> wall_shear( hzeta_nodes, 1 );
    for ( std::size_t i=0; i < Param::N + 1; ++i )
    {
      wall_shear( i, 0 ) = -( 3 * Q_output(i,0,U+4) - 4 * Q_output(i,1,U+4) 
            + Q_output(i,2,U+4) ) * Mesh::Yd(0.0)/(2*dY);    
    }

    wall_shear.output( Example::output_path + "Wall_shear_zeta0_" + zeta0_str + ".dat" );


    cout << "  * zeta0 = " << Param::zeta0 << ", A = " << Param::A << endl;
    Param::zeta0 += 1.0; // 0.5 is best for blowing
    Param::zeta0_2 = Param::zeta0 * Param::zeta0;

  }while( Param::zeta0 < 40.5 );


  cout << "FINISHED" << endl;
}
