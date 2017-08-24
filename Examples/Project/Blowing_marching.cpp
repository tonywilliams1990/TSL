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
// Either NO_SPEED_UP or SPEED_UP for normal or reusing the factorised matrix
#define SPEED_UP

namespace TSL
{
    namespace Param
    {
      double hzeta_right( 30.0 );        // Size of the domain in the zeta_hat direction
      double eta_top( 30.0 );           // Size of the domain in the eta direction
      const std::size_t N( 200 );       // Number of intervals in the zeta_hat direction
      const std::size_t M( 200 );       // Number of intervals in the eta direction
      const std::size_t Nvar( 4 );      // Number of variables
      double beta( 0.0 );               // Hartree parameter
      double KB( 0.0 );                 // Base flow transpiration ( +ve = blowing )
      double zeta0( 1.0 );              // Ridge/transpiration width
      double zeta0_2 = zeta0 * zeta0;   // Square of the ridge/transpiration width
      double A( 0.0 );                  // Mass flux parameter
      double K( 1.0 );                  // Transpiration parameter ( +ve = blowing )
      double gamma( 20.0 );             // Steepness factor
      double x_step( 0.01 );            // Size of the downstream grid spacing
      const std::size_t N_x( 2000 );    // Number of intervals in the downstream x direction
      //=> x_max = x_step * N_x
      double x_max( x_step * N_x );     // Maximum downstream location
      double x( 0.0 );                  // Current downstream x location
      double x_d( 5.0 );                // Downstream location of the injection
      //TODO change the output description at the start + include more info

    } // End of namespace Param

    namespace Example
    {
      std::string output_path;          // Output path

      std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
      {
        // Return the column number for the kth variable at node (i,j)
        return Param::Nvar * ( i * ( Param::M + 1 ) + j ) + k;
      }

      double Phi_w( const double& hzeta, const double& x )
      {
        /*double x_pow = pow( x, 2. * ( 1. - Param::beta ) / ( 2. - Param::beta ) );
        double e_x = exp( - ( x - Param::x_d ) * ( x - Param::x_d ) );
        double e_z = exp( - ( 2. - Param::beta ) * x_pow * Param::zeta0_2 * hzeta * hzeta );
        return - Param::K * sqrt( ( 2. - Param::beta ) * x_pow ) * e_x * e_z;*/
        /*return - Param::K * sqrt( 2 * x )
              * exp( - 2 * x * Param::zeta0_2 * hzeta * hzeta )
              * exp( - ( x - Param::x_d ) * ( x - Param::x_d ) );*/
        return - Param::K * ( 1.0 - exp( - x * x ) )
                          * exp( - Param::zeta0_2 * hzeta * hzeta );
      }

      double Phi_w_hzeta( const double& hzeta, const double& x )
      {
        /*double x_pow = pow( x, 2. * ( 1. - Param::beta ) / ( 2. - Param::beta ) );
        return - 2 * ( 2. - Param::beta ) * x_pow * Param::zeta0_2 * hzeta * Example::Phi_w( hzeta, x );*/
        /*return - 4 * x * Param::zeta0_2 * hzeta * Example::Phi_w( hzeta, x );*/
        return - 2 * Param::zeta0_2 * hzeta * Example::Phi_w( hzeta, x );
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

      const double a1( 10.0 );
      const double a2( 4.0 );

      // X = a1 + zeta - a1 * exp( -zeta / a2 )
      double X( const double& zeta )
      {
        return a1 + zeta - a1 * std::exp( - zeta / a2 );
      }
      double Xd( const double& zeta )
      {
        return 1 + ( a1 / a2 ) * std::exp( - zeta / a2 );
      }
      double Xdd( const double& zeta )
      {
        return - ( a1 / ( a2 * a2 ) ) * std::exp( - zeta / a2 );
      }

      const double b1( 10.0 );
      const double b2( 4.0 );

      // Y = b1 + zeta - b1 * exp( -zeta / b2 )
      double Y( const double& eta )
      {
        return b1 + eta - b1 * std::exp( - eta / b2 );
      }
      double Yd( const double& eta )
      {
        return 1 + ( b1 / b2 ) * std::exp( - eta / b2 );
      }
      double Ydd( const double& eta )
      {
        return - ( b1 / ( b2 * b2 ) ) * std::exp( - eta / b2 );
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
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{
  cout << "*** ---------- Blowing downstream marching ---------- ***" << endl;
  cout << "  * We are solving using a " << Param::N + 1 << " x " << Param::M + 1
       << " mesh with zeta_hat_inf = " << Param::hzeta_right << " and eta_inf = "
       << Param::eta_top << "." << endl;

  /* ----- Make the output directory ----- */
  std::ostringstream ss;
  ss << "./DATA/Marching_K_" << Param::K << "_beta_" << Param::beta << "_"
     << Param::N_x << "x" << Param::N + 1 << "x" << Param::M + 1 << "_"
     << Param::x_max << "_" << Param::hzeta_right << "_" << Param::eta_top << "/";
  Example::output_path = ss.str();
  int status = mkdir( Example::output_path.c_str(), S_IRWXU );
  if ( status == 0 ) {
  cout << "  * Output directory " + Example::output_path +
          " has been made successfully." << endl;
  }


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

  /* ----- Solve for the perturbation quantities ----- */

  cout << "*** Solving the perturbation equations ***" << endl;
  cout << "  * Perturbation transpiration K = " << Param::K << endl;
  // Current guess states (n+1,G)
  TwoD_node_mesh<double> Q( X_nodes, Y_nodes, 4 );
  // Solution at previous downstream step (n)
  TwoD_node_mesh<double> Q_old( X_nodes, Y_nodes, 4 );

  // We use the mesh below to write the data on the original zeta-eta domain
  TwoD_node_mesh<double> Q_output( hzeta_nodes, eta_nodes, 8 );

  // output a measure of the solution
  TrackerFile metric( Example::output_path + "A_file.dat" );
  metric.push_ptr( &Param::zeta0, "zeta0" );
  metric.push_ptr( &Param::A, "A" );
  metric.push_ptr( &Param::x, "x" );
  double U_eta( 0.0 );
  double eta_half( 0.0 );
  double integral_U2( 0.0 );
  metric.push_ptr( &U_eta, "U_eta(0,0)");
  metric.push_ptr( &eta_half, "eta at which U=1/2 on zeta=0" );
  metric.push_ptr( &integral_U2, "integral U^2 over the cross-section" );
  metric.header();

  // Vector for the RHS of the matrix problem
  Vector<double> B( 4 * N_eta * N_hzeta, 0.0 );

  do                                                    // Iterate over values of zeta_0
  {
  /* Iterate to a solution */
  double max_residual( 0.0 );                           // Maximum residual
#ifdef NO_SPEED_UP
  std::size_t max_iterations( 20 );                     // Maximum number of iterations
#endif
  std::size_t iteration( 0 );                           // Initialise iteration counter
#ifdef SPEED_UP
  std::size_t max_iterations( 200 );                    // Maximum number of iterations
  Eigen::SparseMatrix<double, Eigen::ColMajor, long long> A_Eigen( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta );
  Eigen::SparseLU< Eigen::SparseMatrix<double, Eigen::ColMajor, long long> > solver;
#endif
  do
  {
    SparseMatrix<double> A( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta );
    cout << "  * Assembling sparse matrix problem" << endl;

    Timer timer;
    timer.start();

    using namespace Example;
    std::size_t row( 0 );                               // Initialise row counter

    /* hzeta = 0 boundary ( left boundary => symmetry conditions ) */
      std::size_t i( 0 );

      for ( std::size_t j = 0; j < Param::M + 1 ; ++j )
      {
          double hzeta( hzeta_nodes[ 0 ] );
          double Xd( Mesh::Xd(hzeta) );
          double eta( eta_nodes[ j ] );
          Vector<double> Base( Base_soln.get_interpolated_vars( eta ) );
          /*
          // Phi_hzeta = 0
          A( row, col( i, j, Phi ) )      = -3.*Xd/(2*dX);
          A( row, col( i + 1, j, Phi ) )  =  4.*Xd/(2*dX);
          A( row, col( i + 2, j, Phi ) )  = -1.*Xd/(2*dX);

          B[ row ]                        = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi)
                                            -Q(i+2,j,Phi) )/(2*dX) )
                                            -( Xd*( -3*Q_old(i,j,Phi) + 4*Q_old(i+1,j,Phi)
                                            -Q_old(i+2,j,Phi) )/(2*dX) );
          ++row;

          // Psi = 0
          A( row, col( i, j, Psi ) )      =   1;
          B[ row ]                        = - Q( i, j, Psi )
                                            - Q_old( i, j, Psi);
          ++row;

          // U_hzeta = 0
          A( row, col( i, j, U ) )        = -3.*Xd/(2*dX);
          A( row, col( i + 1, j, U ) )    =  4.*Xd/(2*dX);
          A( row, col( i + 2, j, U ) )    = -1.*Xd/(2*dX);

          B[ row ]                        = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U)
                                            -Q(i+2,j,U) ) / (2*dX) )
                                            -( Xd*( -3*Q_old(i,j,U) + 4*Q_old(i+1,j,U)
                                            -Q_old(i+2,j,U) ) / (2*dX) );
          ++row;

          // Theta = 0
          A( row, col( i, j, Theta ) )    =   1;
          B[ row ]                        = - Q( i, j, Theta )
                                            - Q_old( i, j, Theta );
          ++row;*/

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
      /*double Phi_w( Example::Phi_w( hzeta, Param::x + Param::x_step ) );
      double Phi_w_hzeta( Example::Phi_w_hzeta( hzeta, Param::x + Param::x_step ) );
      double Phi_w_old( Example::Phi_w( hzeta, Param::x ) );
      double Phi_w_old_hzeta( Example::Phi_w_hzeta( hzeta, Param::x ) );*/
      double Phi_w( Example::Phi_w( hzeta, Param::x ) );
      double Phi_w_hzeta( Example::Phi_w_hzeta( hzeta, Param::x ) );

      /* eta = 0 boundary ( bottom boundary ) */
      std::size_t j( 0 );
      double eta( eta_nodes[ j ] );
      double Yd( Mesh::Yd(eta) );

      /*
      // Phi = Phi_w
      A( row, col( i, j, Phi ) )        =   1;
      B[ row ]                          = - Q( i, j, Phi ) - Q_old( i, j, Phi )
                                          + Phi_w + Phi_w_old;
      ++row;
      // Psi = 0
      A( row, col( i, j, Psi ) )        =   1;
      B[ row ]                          = - Q( i, j, Psi ) - Q_old( i, j, Psi );
      ++row;
      // U = 0
      A( row, col( i, j, U ) )          =   1;
      B[ row ]                          = - Q( i, j, U ) - Q_old( i, j, U );
      ++row;
      // Theta - Psi_eta = -( 1 / ( zeta0^2 ) ) * Phi_w_hzeta
      A( row, col( i, j, Theta ) )      =  1;
      A( row, col( i, j, Psi ) )        =  3 * Yd / ( 2 * dY );
      A( row, col( i, j + 1, Psi ) )    = -4 * Yd / ( 2 * dY );
      A( row, col( i, j + 2, Psi ) )    =  1 * Yd / ( 2 * dY );
      B[ row ]        = - Q( i, j, Theta ) + Yd * ( -3 * Q( i, j, Psi )
                        + 4 * Q( i, j + 1, Psi ) - Q( i, j + 2, Psi ) )
                        / ( 2 * dY )
                        - Q_old( i, j, Theta ) + Yd * ( -3 * Q_old( i, j, Psi )
                        + 4 * Q_old( i, j + 1, Psi ) - Q_old( i, j + 2, Psi ) )
                        / ( 2 * dY )
                        - ( 1. / ( Param::zeta0_2 ) ) * Phi_w_hzeta
                        - ( 1. / ( Param::zeta0_2 ) ) * Phi_w_old_hzeta;
      ++row;*/

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
        double eta( eta_nodes[ j ] );
        //double x( Param::x );
        //double dx( Param::x_step );
        double Yd( Mesh::Yd( eta ) );
        double Ydd( Mesh::Ydd( eta ) );
        Vector<double> Base( Base_soln.get_interpolated_vars( eta ) );

        // Laplacian coefficients for finite-differencing
        // X(i,j-1)
        double laplace_1 =  Yd*Yd/(dY*dY) - Ydd/(2.*dY);
        // X(i-1,j)
        double laplace_3 = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) ) / ( Param::zeta0_2 );
        // X(i,j)
        double laplace_4 = -2.*( Yd*Yd / (dY*dY)
                           + Xd*Xd/( Param::zeta0_2 * dX * dX ) );
        // X(i+1,j)
        double laplace_5 = ( Xdd/(2.*dX) + Xd*Xd/(dX*dX) ) / ( Param::zeta0_2 );
        // X(i,j+1)
        double laplace_7 =  Yd*Yd/(dY*dY) + Ydd/(2.*dY);

        // Guessed/known components and various derivative values
        Vector<double> Guess( Q.get_nodes_vars( i, j ) );
        Vector<double> Guess_eta( ( Q.get_nodes_vars( i, j + 1 )
                                  - Q.get_nodes_vars( i, j - 1 ) )
                                  * ( Yd /( 2 * dY )) );
        Vector<double> Guess_hzeta( ( Q.get_nodes_vars( i + 1, j )
                                    - Q.get_nodes_vars( i - 1, j ) )
                                    * ( Xd /( 2 * dX )) );
        Vector<double> Guess_laplace( Q.get_nodes_vars( i, j - 1 ) * laplace_1
                                +  Q.get_nodes_vars( i - 1, j ) * laplace_3
                                +  Q.get_nodes_vars( i, j ) * laplace_4
                                +  Q.get_nodes_vars( i + 1, j ) * laplace_5
                                +  Q.get_nodes_vars( i, j + 1 ) * laplace_7 );

        // Components at previous x step and various derivative values
        Vector<double> Old( Q_old.get_nodes_vars( i, j ) );
        Vector<double> Old_eta( ( Q_old.get_nodes_vars( i, j + 1 )
                                - Q_old.get_nodes_vars( i, j - 1 ) )
                                * ( Yd /( 2 * dY )) );
        Vector<double> Old_hzeta( ( Q_old.get_nodes_vars( i + 1, j )
                                  - Q_old.get_nodes_vars( i - 1, j ) )
                                  * ( Xd /( 2 * dX )) );
        Vector<double> Old_laplace( Q_old.get_nodes_vars( i, j - 1 ) * laplace_1
                              +  Q_old.get_nodes_vars( i - 1, j ) * laplace_3
                              +  Q_old.get_nodes_vars( i, j ) * laplace_4
                              +  Q_old.get_nodes_vars( i + 1, j ) * laplace_5
                              +  Q_old.get_nodes_vars( i, j + 1 ) * laplace_7 );

        //////////////////
        // Phi equation //
        //////////////////

        // Laplacian of Phi
        A( row, col( i, j - 1, Phi ) )      = laplace_1;
        A( row, col( i - 1, j, Phi ) )      = laplace_3;
        A( row, col( i, j, Phi ) )          = laplace_4;
        A( row, col( i + 1, j, Phi ) )      = laplace_5;
        A( row, col( i, j + 1, Phi ) )      = laplace_7;
        // -(2-beta) * U_eta
        A( row, col( i, j + 1, U ) )        = -( 2. - Param::beta )*Yd/( 2 * dY );
        A( row, col( i, j - 1, U ) )        =  ( 2. - Param::beta )*Yd/( 2 * dY );
        // Theta_hzeta
        A( row, col( i + 1, j, Theta ) )    =  Xd / ( 2 * dX );
        A( row, col( i - 1, j, Theta ) )    = -Xd / ( 2 * dX );
        // -(4x/dx) * U_eta
        A( row, col( i, j + 1, U ) )       += -( 4 * Param::x / Param::x_step )
                                              * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )       +=  ( 4 * Param::x / Param::x_step )
                                              * Yd / ( 2 * dY );

        // Residual
        B[ row ]      = - Guess_laplace[ Phi ]
                        + ( 2. - Param::beta ) * Guess_eta[ U ]
                        - Guess_hzeta[ Theta ]
                        - Old_laplace[ Phi ]
                        + ( 2. - Param::beta ) * Old_eta[ U ]
                        - Old_hzeta[ Theta ]
                        + ( 4 * Param::x / Param::x_step )
                        * ( Guess_eta[ U ] - Old_eta[ U ] );
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
        // -(2-beta) * U_hzeta / zeta0^2
        A( row, col( i + 1, j, U ) )        = - ( 2. - Param::beta ) * Xd
                                              / ( 2. * dX * Param::zeta0_2 );
        A( row, col( i - 1, j, U ) )        =   ( 2. - Param::beta ) * Xd
                                              / ( 2. * dX * Param::zeta0_2 );
        // -Theta_eta
        A( row, col( i, j + 1, Theta ) )    = - Yd / ( 2 * dY );
        A( row, col( i, j - 1, Theta ) )    =   Yd / ( 2 * dY );
        // -(4x/dx) * U_hzeta / zeta0^2
        A( row, col( i + 1, j, U ) )       += - ( 4 * Param::x / Param::x_step )
                                              * Xd / ( 2 * dX * Param::zeta0_2 );
        A( row, col( i - 1, j, U ) )       +=   ( 4 * Param::x / Param::x_step )
                                              * Xd / ( 2 * dX * Param::zeta0_2 );

        // Residual
        B[ row ]      = - Guess_laplace[ Psi ] + ( 2. - Param::beta )
                        * ( Guess_hzeta[ U ] ) / ( Param::zeta0_2 )
                        + Guess_eta[ Theta ]
                        - Old_laplace[ Psi ] + ( 2. - Param::beta )
                        * ( Old_hzeta[ U ] ) / ( Param::zeta0_2 )
                        + Old_eta[ Theta ]
                        + ( 4 * Param::x / Param::x_step ) * ( Guess_hzeta[ U ]
                        - Old_hzeta[ U ] ) / Param::zeta0_2;
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
        // -beta * ( Uold + UG + 2*UB ) * U
        A( row, col( i, j, U ) )           += - Param::beta * ( Old[ U ]
                                              + Guess[ U ] + 2 * Base[ UB ] );
        // -(4x/dx) * ( UB + UG ) * U
        A( row, col( i, j, U ) )           += - ( 4 * Param::x / Param::x_step )
                                              * ( Base[ UB ] + Guess[ U ] );
        // 0.5 * ( Psiold + PsiG + 2 * hzeta * PsiB ) * U_hzeta
        A( row, col( i + 1, j, U ) )       +=   0.5 * ( Old[ Psi ] + Guess[ Psi ]
                                              + 2 * hzeta * Base[ PsiB ] )
                                              * Xd / ( 2 * dX );
        A( row, col( i - 1, j, U ) )       += - 0.5 * ( Old[ Psi ] + Guess[ Psi ]
                                              + 2 * hzeta * Base[ PsiB ] )
                                              * Xd / ( 2 * dX );
        // 0.5 * ( Phiold + PhiG + 2 * PhiB ) * U_eta
        A( row, col( i, j + 1, U ) )       +=   0.5 * ( Old[ Phi ] + Guess[ Phi ]
                                              + 2 * Base[ PhiB ] )
                                              * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )       += - 0.5 * ( Old[ Phi ] + Guess[ Phi ]
                                              + 2 * Base[ PhiB ] )
                                              * Yd / ( 2 * dY );
        // 0.5 * ( Uold_hzeta + UG_hzeta ) * Psi
        A( row, col( i, j, Psi ) )          =   0.5 * ( Old_hzeta[ U ]
                                              + Guess_hzeta[ U ] );
        // 0.5 * ( Uold_eta + UG_eta + 2 * UB' ) * Phi
        A( row, col( i, j, Phi ) )          =   0.5 * ( Old_eta[ U ]
                                              + Guess_eta[ U ] + 2 * Base[ UBd ] );

        // Residual
        B[ row ]        = - Old_laplace[ U ] - Guess_laplace[ U ]
                          + Param::beta * ( 2. * Base[ UB ] + 0.5 * ( Old[ U ]
                          + Guess[ U ] ) ) * ( Old[ U ] + Guess[ U ] )
                          - 0.5 * ( Old[ Psi ] + Guess[ Psi ]
                          + 2 * hzeta * Base[ PsiB ] )
                          * ( Old_hzeta[ U ] + Guess_hzeta[ U ] )
                          - 0.5 * ( Old[ Phi ] + Guess[ Phi ]
                          + 2 * Base[ PhiB ] )
                          * ( Old_eta[ U ] + Guess_eta[ U ] )
                          - Base[ UBd ] * ( Old[ Phi ] + Guess[ Phi ] )
                          + ( 4 * Param::x / Param::x_step ) * ( Base[ UB ]
                          + 0.5 * ( Old[ U ] + Guess[ U ] ) )
                          * ( Guess[ U ] - Old[ U ] ) ;
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
        // -(1-beta) * hzeta * (2*UB + Uold + UG) * U_eta
        A( row, col( i, j + 1, U ) )         = - ( 1. - Param::beta ) * hzeta
                                                 * ( 2 * Base[ UB ] + Old[ U ]
                                                 +  Guess[ U ] )
                                                 * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )         =   ( 1. - Param::beta ) * hzeta
                                                 * ( 2 * Base[ UB ] + Old[ U ]
                                                 +  Guess[ U ] )
                                                 * Yd / ( 2 * dY );
        // -(1-beta) * hzeta * (2*UB' + Uold_eta + UG_eta) * U
        A( row, col( i, j, U ) )             = - ( 1. - Param::beta )
                                                 * hzeta * ( 2 * Base[ UBd ]
                                                 + Old_eta[ U ]
                                                 + Guess_eta[ U ] );
        // ((1-beta) * eta * ( Uold_hzeta + UG_hzeta ) / zeta0^2) * U
        A( row, col( i, j, U ) )            +=   ( 1. - Param::beta ) * eta
                                                 * ( Old_hzeta[ U ]
                                                 + Guess_hzeta[ U ] )
                                                 / Param::zeta0_2;
        // ((1-beta) * eta * (2*UB + Uold + UG) / zeta0^2 ) * U_hzeta
        A( row, col( i + 1, j, U ) )         =   ( 1. - Param::beta ) * eta
                                                 * ( 2 * Base[ UB ] + Old[ U ]
                                                 + Guess[ U ] ) * Xd
                                                 / ( 2 * dX * Param::zeta0_2 );
        A( row, col( i - 1, j, U ) )         = - ( 1. - Param::beta ) * eta
                                                 * ( 2 * Base[ UB ] + Old[ U ]
                                                 + Guess[ U ] ) * Xd
                                                 / ( 2 * dX * Param::zeta0_2 );
        // 0.5 * ( 2 * PhiB + Phiold + PhiG ) * Theta_eta
        A( row, col( i, j + 1, Theta ) )    +=   0.5 * ( 2 * Base[ PhiB ]
                                                 + Old[ Phi ] + Guess[ Phi ] )
                                                 * Yd / ( 2 * dY );
        A( row, col( i, j - 1, Theta ) )    += - 0.5 * ( 2 * Base[ PhiB ]
                                                 + Old[ Phi ] + Guess[ Phi ] )
                                                 * Yd / ( 2 * dY );
        // 0.5 * ( 2 * hzeta * ThetaB' + Thetaold_eta + ThetaG_eta ) * Phi
        A( row, col( i, j, Phi ) )           =   0.5 * ( 2 * hzeta
                                                 * Base[ ThetaBd ]
                                                 + Old_eta[ Theta ]
                                                 + Guess_eta[ Theta ] );
        // 0.5 * ( 2 * hzeta * PsiB + Psiold + PsiG ) * Theta_hzeta
        A( row, col( i + 1, j, Theta ) )    +=   0.5 * ( 2 * hzeta * Base[ PsiB ]
                                                 + Old[ Psi ] + Guess[ Psi ] )
                                                 * Xd / ( 2 * dX );
        A( row, col( i - 1, j, Theta ) )    += - 0.5 * ( 2 * hzeta * Base[ PsiB ]
                                                 + Old[ Psi ] + Guess[ Psi ] )
                                                 * Xd / ( 2 * dX );
        // 0.5 * ( 2 * ThetaB + Thetaold_hzeta + ThetaG_hzeta) * Psi
        A( row, col( i, j, Psi ) )           =   0.5 * ( 2 * Base[ ThetaB ]
                                                 + Old_hzeta[ Theta ]
                                                 + Guess_hzeta[ Theta ] );
        // 0.5 * (2-beta) * ( 2*UB + Uold + UG ) * Theta
        A( row, col( i, j, Theta ) )        +=   0.5 * ( 2. - Param::beta )
                                                 * ( 2 * Base[ UB ] + Old[ U ]
                                                 + Guess[ U ] );
        // 0.5 * (2-beta) * ( 2 * hzeta * ThetaB + Thetaold + ThetaG ) * U
        A( row, col( i, j, U ) )            +=   0.5 * ( 2. - Param::beta )
                                                 * ( 2 * hzeta * Base[ Theta ]
                                                 + Old[ Theta ]+ Guess[ Theta ] );

        // -(2x/dx) * ( 2 * UBd + Uold_eta + UG_eta ) * Psi
        A( row, col( i, j, Psi ) )          += - ( 2 * Param::x / Param::x_step )
                                                 * ( 2 * Base[ UBd ]
                                                 + Old_eta[ U ]
                                                 + Guess_eta[ U ] );
        // -(2x/dx) * ( PsiG - Psiold ) * U_eta
        A( row, col( i, j + 1, U ) )        += - ( 2 * Param::x / Param::x_step )
                                                 * ( Guess[ Psi ] - Old[ Psi ] )
                                                 * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )        +=   ( 2 * Param::x / Param::x_step )
                                                 * ( Guess[ Psi ] - Old[ Psi ] )
                                                 * Yd / ( 2 * dY );
        // ( (2x/dx) * ( Uold_hzeta + UG_hzeta ) / zeta0^2 ) * Phi
        A( row, col( i, j, Phi ) )          +=   ( 2 * Param::x / Param::x_step )
                                                 * ( Old_hzeta[ U ]
                                                 + Guess_hzeta[ U ] )
                                                 / Param::zeta0_2;
        // ( (2x/dx) * ( PhiG - Phiold ) / zeta0^2 ) * U_hzeta
        A( row, col( i + 1, j, U ) )        +=   ( 2 * Param::x / Param::x_step )
                                                 * ( Guess[ Phi ]
                                                 - Old[ Phi ] ) * Xd
                                                 / ( 2 * dX * Param::zeta0_2 );
        A( row, col( i - 1, j, U ) )        += - ( 2 * Param::x / Param::x_step )
                                                 * ( Guess[ Phi ]
                                                 - Old[ Phi ] ) * Xd
                                                 / ( 2 * dX * Param::zeta0_2 );
        // -(2x/dx) * ( 2 * UB + Uold + UG ) * Theta
        A( row, col( i, j, Theta ) )        += - ( 2 * Param::x / Param::x_step )
                                                 * ( 2 * Base[ UB ] + Old[ U ]
                                                 + Guess[ U ] );
        // -(2x/dx) * ( ThetaG - Thetaold ) * U
        A( row, col( i, j, U ) )            += - ( 2 * Param::x / Param::x_step )
                                                 * ( Guess[ Theta ]
                                                 - Old[ Theta ] );
        // (2x/dx) * ( 2 * hzeta * ThetaB + Thetaold + ThetaG ) * U
        A( row, col( i, j, U ) )            +=   ( 2 * Param::x / Param::x_step )
                                                 * ( 2 * hzeta * Base[ ThetaB ]
                                                 + Old[ Theta ]
                                                 + Guess[ Theta ] );
        // (2x/dx) * ( UG - Uold ) * Theta
        A( row, col( i, j, Theta ) )        +=   ( 2 * Param::x / Param::x_step )
                                                 * ( Guess[ U ] - Old[ U ] );

        // Residual
        B[ row ]      = - Old_laplace[ Theta ] - Guess_laplace[ Theta ]
                        + ( 2 * Param::x / Param::x_step ) * ( 2 * Base[ UBd ]
                        + Old_eta[ U ] + Guess_eta[ U ] ) * ( Guess[ Psi ] - Old[ Psi ] )
                        - ( ( 2 * Param::x / Param::x_step ) * ( Old_hzeta[ U ]
                        + Guess_hzeta[ U ] ) * ( Guess[ Phi ] - Old[ Phi ] )
                        / Param::zeta0_2 )
                        + ( 2 * Param::x / Param::x_step ) * ( 2 * Base[ UB ]
                        + Old[ U ] + Guess[ U ] ) * ( Guess[ Theta ] - Old[ Theta ] )
                        - ( 2 * Param::x / Param::x_step ) * ( 2 * hzeta
                        * Base[ ThetaB ] + Old[ Theta ] + Guess[ Theta ] )
                        * ( Guess[ U ] - Old[ U ] )
                        + ( 1. - Param::beta ) * hzeta * ( 2 * Base[ UB ]
                        + Old[ U ] + Guess[ U ] ) * ( Old_eta[ U ] + Guess_eta[ U ] )
                        + 2 * ( 1. - Param::beta ) * hzeta * Base[ UBd ]
                        * ( Old[ U ] + Guess[ U ] )
                        - (( 1. - Param::beta ) * eta * ( 2 * Base[ UB ]
                        + Old[ U ] + Guess[ U ] ) * ( Old_hzeta[ U ]
                        + Guess_hzeta[ U ] ) / Param::zeta0_2 )
                        - 0.5 * ( 2 * Base[ PhiB ] + Old[ Phi ] + Guess[ Phi ] )
                        * ( Old_eta[ Theta ] + Guess_eta[ Theta ] )
                        - hzeta * Base[ ThetaBd ] * ( Old[ Phi ] + Guess[ Phi ] )
                        - 0.5 * ( 2 * Base[ ThetaB ] + Old_hzeta[ Theta ]
                        + Guess_hzeta[ Theta ] ) * ( Old[ Psi ] + Guess[ Psi ] )
                        - hzeta * Base[ PsiB ] * ( Old_hzeta[ Theta ]
                        + Guess_hzeta[ Theta ] )
                        - ( 2. - Param::beta ) * ( 0.5 * ( 2 * Base[ UB ]
                        + Old[ U ] + Guess[ U ] ) * ( Old[ Theta ]
                        + Guess[ Theta ] ) + hzeta * Base[ ThetaB ]
                        * ( Old[ U ] + Guess[ U ] ) );
        ++row;

      }

      /* eta = eta_inf boundary ( top boundary ) */
        j = Param::M ;
        eta = eta_nodes[ j ];
        Yd = Mesh::Yd( eta );

        /*
        // Phi_eta*( eta^2 + zeta_0^2*hzeta^2) + [ 2*eta - (eta^2 + zeta_0^2*hzeta^2)/eta ]*Phi = 0
        A( row, col( i, j, Phi ) )        =   3.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 1, Phi ) )    = - 4.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 2, Phi ) )    =   1.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j, Phi ) )       +=   2.0 * eta - ((eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / eta );

        B[ row ]                          = - ( 3 * Q( i, j, Phi ) - 4 * Q( i, j-1, Phi )
                                            + Q( i, j-2, Phi ) ) * Yd * (eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / ( 2 * dY )
                                            + (((eta * eta + Param::zeta0_2 * hzeta * hzeta) / eta )
                                            - 2.0 * eta) * Q( i, j, Phi )
                                            - ( 3 * Q_old( i, j, Phi ) - 4 * Q_old( i, j-1, Phi )
                                            + Q_old( i, j-2, Phi ) ) * Yd * (eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / ( 2 * dY )
                                            + (((eta * eta + Param::zeta0_2 * hzeta * hzeta) / eta )
                                            - 2.0 * eta) * Q_old( i, j, Phi );
        ++row;

        // Psi_eta*( eta^2 + zeta_0^2*hzeta^2) + 2 * eta * Psi = 0
        A( row, col( i, j, Psi ) )        =   3.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 1, Psi ) )    = - 4.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 2, Psi ) )    =   1.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j, Psi ) )       +=   2 * eta;

        B[ row ]                          =  - (( 3 * Q( i, j, Psi ) - 4 * Q( i, j-1, Psi )
                                             + Q( i, j-2, Psi ) ) * Yd * (eta * eta
                                             + Param::zeta0_2 * hzeta * hzeta) / ( 2 * dY ))
                                             - 2 * eta * Q( i, j, Psi )
                                             - (( 3 * Q_old( i, j, Psi ) - 4 * Q_old( i, j-1, Psi )
                                             + Q_old( i, j-2, Psi ) ) * Yd * (eta * eta
                                             + Param::zeta0_2 * hzeta * hzeta) / ( 2 * dY ))
                                             - 2 * eta * Q_old( i, j, Psi );
        ++row;

        // U = 0
        A( row, col( i, j, U ) )            =   1;
        B[ row ]                            = - Q( i, j, U ) - Q_old( i, j, U );
        ++row;

        // Theta = 0
        A( row, col( i, j, Theta ) )        =   1;
        B[ row ]                            = - Q( i, j, Theta )
                                              - Q_old( i, j, Theta );
        ++row;*/

        // Phi_eta*( eta^2 + zeta_0^2*hzeta^2) + [ 2*eta - (eta^2 + zeta_0^2*hzeta^2)/eta ]*Phi = 0
        A( row, col( i, j, Phi ) )        =   3.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 1, Phi ) )    = - 4.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 2, Phi ) )    =   1.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j, Phi ) )       +=   2.0 * eta - ((eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / eta );

        B[ row ]                          = - ( 3 * Q( i, j, Phi ) - 4 * Q( i, j-1, Phi )
                                            + Q( i, j-2, Phi ) ) * Yd * (eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / ( 2 * dY )
                                            + (((eta * eta + Param::zeta0_2 * hzeta * hzeta) / eta )
                                            - 2.0 * eta) * Q( i, j, Phi );
        ++row;

        // Psi_eta*( eta^2 + zeta_0^2*hzeta^2) + 2 * eta * Psi = 0
        A( row, col( i, j, Psi ) )        =   3.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 1, Psi ) )    = - 4.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j - 2, Psi ) )    =   1.0 * Yd * ( eta * eta
                                            + Param::zeta0_2 * hzeta * hzeta) / (2*dY);
        A( row, col( i, j, Psi ) )       +=   2 * eta;

        B[ row ]                          =  - (( 3 * Q( i, j, Psi ) - 4 * Q( i, j-1, Psi )
                                             + Q( i, j-2, Psi ) ) * Yd * (eta * eta
                                             + Param::zeta0_2 * hzeta * hzeta) / ( 2 * dY ))
                                             - 2 * eta * Q( i, j, Psi );
        ++row;

        // U = 0
        A( row, col( i, j, U ) )            =   1;
        B[ row ]                            = - ( Q( i, j, U ) );
        ++row;

        // Theta = 0
        A( row, col( i, j, Theta ) )        =   1;
        B[ row ]                            = - ( Q( i, j, Theta ) );
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
      /*
      // (eta^2 + zeta_0^2 * hzeta^2) * Phi_hzeta + 2 * zeta_0^2 * hzeta * Phi = 0
      A( row, col( i, j, Phi ) )          =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 3. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 1, j, Phi ) )      = - (eta*eta + Param::zeta0_2*hzeta*hzeta) * 4. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 2, j, Phi ) )      =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 1. * Xd
                                            / ( 2 * dX );
      A( row, col( i, j, Phi ) )         +=   2 * Param::zeta0_2 * hzeta;

      B[ row ]        = - (eta*eta + Param::zeta0_2*hzeta*hzeta) * ( 3 * Q( i, j, Phi)
                        - 4 * Q( i - 1, j, Phi) + Q( i - 2, j, Phi) ) * Xd / ( 2 * dX )
                        - 2 * Param::zeta0_2 * hzeta * Q( i, j, Phi )
                        - (eta*eta + Param::zeta0_2*hzeta*hzeta) * ( 3 * Q_old( i, j, Phi)
                        - 4 * Q_old( i - 1, j, Phi) + Q_old( i - 2, j, Phi) ) * Xd / ( 2 * dX )
                        - 2 * Param::zeta0_2 * hzeta * Q_old( i, j, Phi );
      ++row;

      // (eta^2 + zeta_0^2 * hzeta^2)*Psi_hzeta + (2*zeta_0^2*hzeta-(eta^2 + zeta_0^2*hzeta^2)/hzeta)*Psi = 0
      A( row, col( i, j, Psi ) )          =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 3. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 1, j, Psi ) )      = - (eta*eta + Param::zeta0_2*hzeta*hzeta) * 4. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 2, j, Psi ) )      =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 1. * Xd
                                            / ( 2 * dX );
      A( row, col( i, j, Psi ) )         +=   2 * Param::zeta0_2 * hzeta
                                            - ((eta*eta + Param::zeta0_2*hzeta*hzeta) / hzeta);


      B[ row ]        = - ((eta*eta + Param::zeta0_2*hzeta*hzeta) * ( 3 * Q( i, j, Psi )
                        - 4 * Q( i - 1, j, Psi ) + Q( i - 2, j, Psi) ) * Xd / ( 2 * dX ))
                        - 2 * Param::zeta0_2 * hzeta  * Q( i, j, Psi)
                        + ((eta*eta + Param::zeta0_2*hzeta*hzeta) / hzeta)  * Q( i, j, Psi)
                        - ((eta*eta + Param::zeta0_2*hzeta*hzeta) * ( 3 * Q_old( i, j, Psi )
                        - 4 * Q_old( i - 1, j, Psi ) + Q_old( i - 2, j, Psi) ) * Xd / ( 2 * dX ))
                        - 2 * Param::zeta0_2 * hzeta  * Q_old( i, j, Psi)
                        + ((eta*eta + Param::zeta0_2*hzeta*hzeta) / hzeta)  * Q_old( i, j, Psi);
      ++row;

      // hzeta * U_hzeta + 2 * U = 0
      A( row, col( i, j, U ) )            =   hzeta * 3. * Xd / ( 2 * dX ) + 2.;
      A( row, col( i - 1, j, U ) )        = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, U ) )        =   hzeta * 1. * Xd / ( 2 * dX );

      B[ row  ]       = - hzeta * ( 3 * Q( i, j, U ) - 4 * Q( i - 1, j, U )
                        + Q( i - 2, j, U) ) * Xd / ( 2 * dX ) - 2 * Q( i, j, U )
                        - hzeta * ( 3 * Q_old( i, j, U ) - 4 * Q_old( i - 1, j, U )
                        + Q_old( i - 2, j, U) ) * Xd / ( 2 * dX )
                        - 2 * Q_old( i, j, U );
      ++row;

      // hzeta * Theta_hzeta + Theta = 0
      A( row, col( i, j, Theta )  )       =   hzeta * 3. * Xd / ( 2 * dX ) + 1.;
      A( row, col( i - 1, j, Theta ) )    = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, Theta ) )    =   hzeta * 1. * Xd / ( 2 * dX );

      B[ row ]        = - hzeta * ( 3 * Q( i, j, Theta ) - 4 * Q( i - 1, j, Theta )
                        + Q( i - 2, j, Theta ) ) * Xd / ( 2 * dX ) - Q( i, j, Theta )
                        - hzeta * ( 3 * Q_old( i, j, Theta ) - 4 * Q_old( i - 1, j, Theta )
                        + Q_old( i - 2, j, Theta ) ) * Xd / ( 2 * dX ) - Q_old( i, j, Theta );
      ++row;*/

      // (eta^2 + zeta_0^2 * hzeta^2) * Phi_hzeta + 2 * zeta_0^2 * hzeta * Phi = 0
      A( row, col( i, j, Phi ) )          =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 3. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 1, j, Phi ) )      = - (eta*eta + Param::zeta0_2*hzeta*hzeta) * 4. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 2, j, Phi ) )      =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 1. * Xd
                                            / ( 2 * dX );
      A( row, col( i, j, Phi ) )         +=   2 * Param::zeta0_2 * hzeta;

      B[ row ]        = - (eta*eta + Param::zeta0_2*hzeta*hzeta) * ( 3 * Q( i, j, Phi)
                        - 4 * Q( i - 1, j, Phi) + Q( i - 2, j, Phi) ) * Xd / ( 2 * dX )
                        - 2 * Param::zeta0_2 * hzeta * Q( i, j, Phi );
      ++row;

      // (eta^2 + zeta_0^2 * hzeta^2)*Psi_hzeta + (2*zeta_0^2*hzeta-(eta^2 + zeta_0^2*hzeta^2)/hzeta)*Psi = 0
      A( row, col( i, j, Psi ) )          =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 3. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 1, j, Psi ) )      = - (eta*eta + Param::zeta0_2*hzeta*hzeta) * 4. * Xd
                                            / ( 2 * dX );
      A( row, col( i - 2, j, Psi ) )      =   (eta*eta + Param::zeta0_2*hzeta*hzeta) * 1. * Xd
                                            / ( 2 * dX );
      A( row, col( i, j, Psi ) )         +=   2 * Param::zeta0_2 * hzeta
                                            - ((eta*eta + Param::zeta0_2*hzeta*hzeta) / hzeta);


      B[ row ]        = - ((eta*eta + Param::zeta0_2*hzeta*hzeta) * ( 3 * Q( i, j, Psi )
                        - 4 * Q( i - 1, j, Psi ) + Q( i - 2, j, Psi) ) * Xd / ( 2 * dX ))
                        - 2 * Param::zeta0_2 * hzeta  * Q( i, j, Psi)
                        + ((eta*eta + Param::zeta0_2*hzeta*hzeta) / hzeta)  * Q( i, j, Psi);
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

    max_residual = B.norm_inf();
    cout << "***                                              Maximum residual = "
         << B.norm_inf() << endl;

    Vector<double> x;
#ifdef SPEED_UP
    // Convert things to Eigen to see if we can get some speed benefits
    // Only decompose the matrix on the first iteration as it can be reused after
    // This means more iterations but they should be quicker
    if ( iteration == 0 )
    {
      A_Eigen = A.convert_to_Eigen();
      solver.compute( A_Eigen );
    }

    Eigen::Matrix<double, -1, 1> B_Eigen( 4 * N_eta * N_hzeta );
    B_Eigen = B.convert_to_Eigen_Matrix();

    Eigen::Matrix<double, -1, 1> x_Eigen( 4 * N_eta * N_hzeta );
    x_Eigen = solver.solve( B_Eigen );
    x = x.convert_to_Vector( x_Eigen );
#endif
#ifdef NO_SPEED_UP
    x = A.solve( B );
#endif
    B = x; //TODO do we need to do this - could just use x to update Q(i,j,Var) and max correction
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

    cout << "***    Iteration = " << iteration
         << "    Maximum correction = " << B.norm_inf() << endl;

    ++iteration;
  }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) ); // End iteration

  if ( iteration >= max_iterations )
  {
    cout << "STOPPED AFTER TOO MANY ITERATIONS" << endl;
    //TODO need to make it so it actually stops!!! break?
  }

  /* push the data back into the unmapped domain */
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      double hzeta=hzeta_nodes[i];
      for ( std::size_t j = 0; j < N_eta; ++j )
      {
        double eta=eta_nodes[j];
        // first 4 values output are the without the underlying base flow
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
    // Output data
    Q_output.dump_gnu( Example::output_path + "Qout_"
                     + Utility::stringify( Param::zeta0, 3 ) + "_x_"
                     + Utility::stringify( Param::x, 3 ) + ".dat" );

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

    Param::A = Q_output( 0, Param::M, Phi ) * Param::eta_top;

    // Integral of U^2 over the cross-section
    integral_U2 = Q_output.square_integral2D( U );

    metric.update();

    // Get the wall shear values as a function of hzeta
    OneD_node_mesh<double> wall_shear( hzeta_nodes, 1 );
    for ( std::size_t i=0; i < Param::N + 1; ++i )
    {
      wall_shear( i, 0 ) = -( 3 * Q_output(i,0,U+4) - 4 * Q_output(i,1,U+4)
            + Q_output(i,2,U+4) ) * Mesh::Yd(0.0)/(2*dY);
    }

    wall_shear.output( Example::output_path + "Wall_shear_zeta0_"
                     + Utility::stringify( Param::zeta0, 3 ) + "_x_"
                     + Utility::stringify( Param::x, 3 ) + ".dat" );

    // Put the current solution into Q_old for next iteration
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      for ( std::size_t j = 0; j < N_eta; ++j )
      {
        Q_old( i, j, Phi )    = Q( i, j, Phi);
        Q_old( i, j, Psi )    = Q( i, j, Psi);
        Q_old( i, j, U )      = Q( i, j, U);
        Q_old( i, j, Theta )  = Q( i, j, Theta);
      }
    }

    // Increment x
    cout << "  * x = " << Param::x << ", A = " << Param::A << endl;
    Param::x += Param::x_step;

  }while( Param::x < Param::x_step * Param::N_x );

  cout << "FINISHED" << endl;
}
