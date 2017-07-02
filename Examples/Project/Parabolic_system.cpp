#include <cassert>
#include "math.h"

#include "Core"

// Enumerations
enum{ f, fd, fdd, g, gd, gdd };                               // Base ODE
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };                 // Base ODE
enum{ Phi, U, Ud, Psi, Theta, Thetad};                        // PDE enumeration

// Either BASE_2D or BASE_3D for 2D or 3D base flows
#define BASE_2D

namespace TSL
{
  namespace Example
  {
    double hzeta_inf( 16 );			      // Size of the domain in the hzeta direction
    double eta_top( 128 );				      // Size of the domain in the eta direction
    unsigned N( 2000 );					      // Number of intervals in the hzeta direction
    unsigned M( 1000 );					      // Number of intervals in the eta direction
    std::string output_path("./DATA/Parabolic_System/");        // Data output path
    double hzeta( hzeta_inf );        // Value of hzeta at the current step
    double K( 2.0 );                  // Blowing intensity
    double beta( 0.1 );               // Hartree parameter
    double gamma( 20.0 );             // Steepness factor

    double Phi_w( const double& hzeta )
    {
        // Return the transpiration function
        return -Example::K * 0.5 * ( 1. - tanh( Example::gamma * ( hzeta - 1. ) ) );
    }

    /* ------------------------- Base-flow ODE (2D/3D) -------------------------*/

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

    /* -------------------- Governing PDE -------------------------*/

    // Parabolic system
    class Parabolic_equations : public Equation_2matrix<double>
    {
    public:
#ifdef BASE_2D
      /// The problem is 3rd order and real
      Parabolic_equations() : Equation_2matrix<double> ( 3 ) {}

      /// Define the system of equations
      void residual_fn( const Vector<double>& z, Vector<double>& f ) const
      {
        // The 3rd order system for ( Phi, U, U' )
        f[ 0 ] =  z[ U ];
        f[ 1 ] =  z[ Ud ];
        f[ 2 ] =  Example::beta * ( z[ U ] * z[ U ] - 1. ) - z[ Phi ] * z[ Ud ];
      }

      void matrix0( const Vector<double>& z, Matrix<double>& m ) const
      {
        // identity matrix
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = 1.0;
        m( 2, 2 ) = 1.0;
      }

      void get_jacobian_of_matrix0_mult_vector( const Vector<double> &state,
                                                const Vector<double> &vec,
                                                Matrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }

      /// Define the unsteady terms by providing the mass matrix
      void matrix1( const Vector<double>& z, Matrix<double>& m ) const
      {
        m( 0, 1 ) = ( 1. - Example::beta ) * Example::hzeta;
        m( 2, 1 ) = ( 1. - Example::beta ) * Example::hzeta * z[ U ];
      }
#endif
#ifdef BASE_3D
      /// The problem is 6th order and real
      Parabolic_equations() : Equation_2matrix<double> ( 6 ) {}

      /// Define the system of equations
      void residual_fn( const Vector<double>& z, Vector<double>& f ) const
      {
        // The 6th order system for ( Phi, U, U', Psi, Theta, Theta' )
        f[ 0 ] =  ( 2. - Example::beta ) * z[ U ] - z[ Psi ];
        f[ 1 ] =  z[ Ud ];
        f[ 2 ] =  Example::beta * ( z[ U ] * z[ U ] - 1. ) - z[ Phi ] * z[ Ud ];
        f[ 3 ] =  z[ Theta ];
        f[ 4 ] =  z[ Thetad ];
        f[ 5 ] =  2. * ( 1. - Example::beta) * z[ U ] * z[ Ud ] - z[ Phi ] * z[ Thetad ]
                - z[ Psi ] * z[ Theta ] - ( 2. - Example::beta ) * z[ U ] * z[ Theta ];
      }

      void matrix0( const Vector<double>& z, Matrix<double>& m ) const
      {
        // identity matrix
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = 1.0;
        m( 2, 2 ) = 1.0;
        m( 3, 3 ) = 1.0;
        m( 4, 4 ) = 1.0;
        m( 5, 5 ) = 1.0;
      }

      void get_jacobian_of_matrix0_mult_vector( const Vector<double> &state,
                                                const Vector<double> &vec,
                                                Matrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }

      /// Define the unsteady terms by providing the mass matrix
      void matrix1( const Vector<double>& z, Matrix<double>& m ) const
      {
        m( 0, 3 ) = Example::hzeta;
        m( 2, 1 ) = Example::hzeta * z[ Psi ];
        m( 5, 4 ) = Example::hzeta * z[ Psi ];
      }
#endif
    };

    // Boundary conditions at eta = 0
    class Parabolic_bottom_BC : public Residual_with_coords<double>
    {
    public:
#ifdef BASE_2D
      // 2 BCs for 3 unknowns and one independent coordinate (time)
      Parabolic_bottom_BC() : Residual_with_coords<double> ( 2, 3, 1 ) {}

      void residual_fn( const Vector<double>& z, Vector<double>& b ) const
      {
        b[ 0 ] = z[ Phi ] - Example::Phi_w( Example::hzeta );
        b[ 1 ] = z[ U ];
      }
#endif
#ifdef BASE_3D
      // 3 BCs for 6 unknowns and one independent coordinate (time)
      Parabolic_bottom_BC() : Residual_with_coords<double> ( 3, 6, 1 ) {}

      void residual_fn( const Vector<double>& z, Vector<double>& b ) const
      {
        b[ 0 ] = z[ Phi ] - Example::Phi_w( Example::hzeta );
        b[ 1 ] = z[ U ];
        b[ 2 ] = z[ Psi ];
      }
#endif
    };

    // Boudnary conditions at eta = eta_inf
    class Parabolic_top_BC : public Residual_with_coords<double>
    {
    public:
#ifdef BASE_2D
      // 2 BCs for 3 unknowns and one independent coordinate (time)
      Parabolic_top_BC() : Residual_with_coords<double> ( 1, 3, 1 ) {}

      void residual_fn( const Vector<double>& z, Vector<double>& b ) const
      {
        b[ 0 ] = z[ U ] - 1.0;
      }
#endif
#ifdef BASE_3D
      // 3 BCs for 6 unknowns and one independent coordinate (time)
      Parabolic_top_BC() : Residual_with_coords<double> ( 3, 6, 1 ) {}

      void residual_fn( const Vector<double>& z, Vector<double>& b ) const
      {
        b[ 0 ] = z[ U ] - 1.0;
        b[ 1 ] = z[ Psi ] - ( 1.0 - Example::beta );
        b[ 2 ] = z[ Theta ];
      }
#endif
    };

  } // end Example namespace
} // end TSL namespace

using namespace TSL;
using namespace std;

int main()
{
  cout << "*** ---------- Parabolic system ---------- ***" << endl;
  //TODO cout information about output directory
  double hzeta_inf = Example::hzeta_inf;
  double eta_top = Example::eta_top;
  // Number of points to solve for
  unsigned N_zeta = Example::N + 1;
  unsigned N_eta  = Example::M + 1;
  // Nodal positions
  Vector<double> eta_nodes, hzeta_nodes;
  eta_nodes.linspace( 0.0, eta_top, N_eta );
  hzeta_nodes.linspace( 0.0, hzeta_inf, N_zeta );//TODO is this right/do we need it?

  /* --------- Solve the base-flow ODE ---------- */

  cout << "*** Solving the base flow ODE ***" << endl;

  // Setup the base flow ODE problem
  Example::Base_Flow::equation base_equation;
  Example::Base_Flow::plate_BC plate_BC;
  Example::Base_Flow::far_BC far_BC;
  base_equation.beta = 0.1;
  plate_BC.KB = 0.0;
  ODE_BVP<double> base( &base_equation, eta_nodes, &plate_BC, &far_BC );

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

  double arc_step( -0.01 );
  double max_arc_step( 0.1 );

  // Solve the system with beta = 0.1 then arc-length continue until beta = Example::beta
  if ( Example::beta >= 0.1 ) { arc_step = 0.01; }

  base.init_arc( &base_equation.beta, arc_step, max_arc_step );
  do
  {
    arc_step = base.arclength_solve( arc_step );
  }while( base_equation.beta < Example::beta );
  base_equation.beta = Example::beta;
  base.solve_bvp();

  // Store the solution in a mesh
  OneD_node_mesh<double> Base_soln( eta_nodes, 6 );
#ifdef BASE_2D
  for (std::size_t j=0; j < N_eta; ++j )
	{
		Base_soln( j, UB )      =   base.solution()( j, fd );
    Base_soln( j, UBd )     =   base.solution()( j, fdd );
    Base_soln( j, PhiB )    =   base.solution()( j, f );
    Base_soln( j, ThetaB )  =   ( 1.0 - Example::beta ) * base.solution()( j, fdd );
    Base_soln( j, ThetaBd ) =   ( 1.0 - Example::beta ) * ( - base.solution()( j, f ) *
                                base.solution()( j, fdd ) - Example::beta * ( 1.0 -
                                base.solution()( j, fd ) * base.solution()( j, fd ) ) );
    Base_soln( j, PsiB )    =   ( 1.0 - Example::beta ) * base.solution()( j, fd );
	}
#endif
#ifdef BASE_3D
  for (std::size_t j=0; j < N_eta; ++j )
	{
		Base_soln( j, UB )      =   base.solution()( j, fd );
    Base_soln( j, UBd )     =   base.solution()( j, fdd );
    Base_soln( j, PhiB )    =   base.solution()( j, f )
                              + ( 2.0 - Example::beta ) * base.solution()( j, g );
    Base_soln( j, ThetaB )  =   ( 1.0 - Example::beta ) * base.solution()( j, fdd )
                              - ( 2.0 - Example::beta ) * base.solution()( j, gdd );
    Base_soln( j, ThetaBd ) =   ( 1.0 - Example::beta ) * ( -(base.solution()( j, f ) +
                                ( 2.0 - Example::beta) * base.solution()( j, g )) *
                                base.solution()( j, fdd ) - Example::beta * ( 1.0 -
                                base.solution()( j, fd ) * base.solution()( j, fd ) ) )
                              - ( 2.0 - Example::beta ) * ( -(base.solution()( j, f ) +
                                ( 2.0 - Example::beta) * base.solution()( j, g )) *
                                base.solution()( j, gdd ) - Example::beta * ( 1.0 -
                                base.solution()( j, gd ) * base.solution()( j, gd ) ) -
                                2.0 * (1.0 - Example::beta ) * (base.solution()( j, fd ) -
                                base.solution()( j, gd )) * base.solution()( j, gd ) );
    Base_soln( j, PsiB )    =   ( 1.0 - Example::beta ) * base.solution()( j, fd )
                              - ( 2.0 - Example::beta ) * base.solution()( j, gd );
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
       << " 3D solution: " << ( 1. - Example::beta ) * Base_soln.integral2(UB)
                                                   - Base_soln.integral2(PsiB) << endl;
  cout << "  * Hartree parameter beta = " << base_equation.beta << endl;
  cout << "  * UB'(eta=0) =" << base.solution()( 0, fdd ) << endl;
  cout << "  * We have solved the ODE problem, it is output to " + Example::output_path +
          "Base_soln.dat" << endl;

  /*----------------- Solve the parabolic system of PDEs ---------------------*/

  cout << "  * Solving the parabolic PDE system using a " << N_zeta
       << " x " << N_eta << " mesh." << endl;

  // Define the system and BCs
  Example::Parabolic_equations problem;
  Example::Parabolic_bottom_BC BC_bottom;
  Example::Parabolic_top_BC BC_top;

  // Define the domain
  //double bottom( 0.0 );
  //double top( Example::eta_top );  

  // hzeta step
  double dhzeta = - ( Example::hzeta_inf ) / (Example::N);
  cout << "  * dhzeta = " << dhzeta << endl;

  // Construct our IBVP
  PDE_IBVP<double> Parabolic( &problem, eta_nodes, &BC_bottom, &BC_top );
  Parabolic.max_itns() = 100;
  Parabolic.coord( ) = hzeta_inf;
  // Initialise the solution using the ODE similarity solution
#ifdef BASE_2D
  for ( unsigned j = 0; j < N_eta; ++j )
  {
    Parabolic.solution()( j, Phi )    = Base_soln( j, PhiB );
    Parabolic.solution()( j, U )      = Base_soln( j, UB );
    Parabolic.solution()( j, Ud )     = Base_soln( j, UBd );
  }

  // Output mesh
  TwoD_node_mesh<double> Q( hzeta_nodes, eta_nodes, 3 );
#endif
#ifdef BASE_3D
  // Initialise the solution using the ODE similarity solution
  for ( unsigned j = 0; j < N_eta; ++j )
  {
    Parabolic.solution()( j, Phi )    = Base_soln( j, PhiB );
    Parabolic.solution()( j, U )      = Base_soln( j, UB );
    Parabolic.solution()( j, Ud )     = Base_soln( j, UBd );
    Parabolic.solution()( j, Psi )    = Base_soln( j, PsiB );
    Parabolic.solution()( j, Theta )  = Base_soln( j, ThetaB );
    Parabolic.solution()( j, Thetad ) = Base_soln( j, ThetaBd );
  }

  // Output mesh
  TwoD_node_mesh<double> Q( hzeta_nodes, eta_nodes, 6 );
#endif

  // Wall shear values and delta transp values ( at eta = eta_inf )
  Vector<double> hzeta_vals;
  Vector<double> wall_shear_vals;
  Vector<double> delta_transp_vals;
  hzeta_vals.push_back( Example::hzeta );
  wall_shear_vals.push_back( Parabolic.solution()( 0, Ud ) );
  delta_transp_vals.push_back( Parabolic.solution()( N_eta - 1, Phi )
                        - Base_soln( N_eta - 1, PhiB ) );

#ifdef BASE_2D
  // Put the solution at hzeta_inf into the output mesh
  for ( unsigned j = 0; j <= Example::M; ++j )
  {
    Q( Example::N, j, Phi )   = Parabolic.solution()( j, Phi );
    Q( Example::N, j, U )     = Parabolic.solution()( j, U );
    Q( Example::N, j, Ud )    = Parabolic.solution()( j, Ud );
  }
#endif
#ifdef BASE_3D
  // Put the solution at hzeta_inf into the output mesh
  for ( unsigned j = 0; j <= Example::M; ++j )
  {
    Q( Example::N, j, Phi )    = Parabolic.solution()( j, Phi );
    Q( Example::N, j, U )      = Parabolic.solution()( j, U );
    Q( Example::N, j, Ud )     = Parabolic.solution()( j, Ud );
    Q( Example::N, j, Psi )    = Parabolic.solution()( j, Psi );
    Q( Example::N, j, Theta )  = Parabolic.solution()( j, Theta );
    Q( Example::N, j, Thetad ) = Parabolic.solution()( j, Thetad );
  }
#endif

  //Example::hzeta += dhzeta;
  Example::hzeta = hzeta_nodes[ N_zeta - 2 ];
  cout << "hzeta = " << Parabolic.coord( ) << endl;

  /* --- Step in hzeta --- */
  std::size_t i = 1;
  do
  {
   // Take a hzeta step
    try
    {
      Parabolic.step2( dhzeta );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }

    // Put the new solution into the output mesh
    for ( unsigned j = 0; j <= Example::M; ++j )
    {
#ifdef BASE_2D
      Q( Example::N - i, j, Phi )   = Parabolic.solution()( j, Phi );
      Q( Example::N - i, j, U )     = Parabolic.solution()( j, U );
      Q( Example::N - i, j, Ud )    = Parabolic.solution()( j, Ud );
#endif
#ifdef BASE_3D
      Q( Example::N - i, j, Phi )   = Parabolic.solution()( j, Phi );
      Q( Example::N - i, j, U )     = Parabolic.solution()( j, U );
      Q( Example::N - i, j, Ud )    = Parabolic.solution()( j, Ud );
      Q( Example::N - i, j, Psi )   = Parabolic.solution()( j, Psi );
      Q( Example::N - i, j, Theta ) = Parabolic.solution()( j, Theta );
      Q( Example::N - i, j, Theta ) = Parabolic.solution()( j, Theta );
#endif
    }
    cout << "hzeta = " << Parabolic.coord( ) << endl;
    hzeta_vals.push_back( Example::hzeta );
    wall_shear_vals.push_back( Parabolic.solution()( 0, Ud ) );
    delta_transp_vals.push_back( Parabolic.solution()( N_eta - 1, Phi )
                          - Base_soln( N_eta - 1, PhiB ) );
    //cout << "delta_transp = " << Parabolic.solution()( N_eta - 1, Phi )
    //                    - Base_soln( N_eta - 1, PhiB ) << endl;
    // Update the value of hzeta
    Example::hzeta += dhzeta;
    //Example::hzeta = hzeta_nodes[ N_zeta - 2 - i ];
    ++i;
  }while( Example::hzeta >= hzeta_nodes[ 0 ] );

  cout << "U'(hzeta, eta=0) = " << Parabolic.solution()( 0, Ud ) << endl;

#ifdef BASE_2D
  // Output the data mesh to a file
  Q.dump_gnu( Example::output_path + "Qout_K_"
            + Utility::stringify( abs(Example::K), 3 ) + "_beta_"
            + Utility::stringify( Example::beta , 3 )+ "_2D.dat" );
  // Wall shear
  OneD_node_mesh<double> Wall_shear( hzeta_vals, 1 );
  Wall_shear.set_vars_from_vector( wall_shear_vals );
  Wall_shear.output( Example::output_path + "Wall_shear_K_"
                   + Utility::stringify( abs( Example::K ), 3 ) + "_beta_"
                   + Utility::stringify( Example::beta , 3 )+ "_2D.dat" );
  // delta_transp
  OneD_node_mesh<double> delta_transp( hzeta_vals, 1 );
  delta_transp.set_vars_from_vector( delta_transp_vals );
  double M_pi;
  M_pi = delta_transp.integral2() * ( 2 / 3.14159265359 );
  cout << "M / pi = " << M_pi << endl;

#endif
#ifdef BASE_3D
  // Output the data mesh to a file
  Q.dump_gnu( Example::output_path + "Qout_K_"
            + Utility::stringify( abs(Example::K), 3 ) + "_beta_"
            + Utility::stringify( Example::beta , 3 )+ "_3D.dat" );

  OneD_node_mesh<double> Wall_shear( hzeta_vals, 1 );
  Wall_shear.set_vars_from_vector( wall_shear_vals );
  Wall_shear.output( Example::output_path + "Wall_shear_K_"
                   + Utility::stringify( abs( Example::K ), 3 ) + "_beta_"
                   + Utility::stringify( Example::beta , 3 )+ "_3D.dat" );
#endif

  cout << "FINISHED SOLUTION" << endl;
}
