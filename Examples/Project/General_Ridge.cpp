#include <cassert>
#include <cmath>

#include "Core"

// Enumerations
enum{ f, fd, fdd, g, gd, gdd };                               // Base ODE
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };                 // Base ODE 
enum{ Ubar, Ubard, Phibar, Thetabar, Thetabard, Psibar };     // Far-field ODE           
enum{ Phi, Psi, U, Theta };                                   // PDE 

// Either BASE_2D or BASE_3D for 2D or 3D base flows
#define BASE_2D
// Either UNIFORM or NONUNIFORM for uniform of non-uniform mesh
#define NONUNIFORM

namespace TSL
{
    namespace Param
    {
      double hzeta_right( 16.0 );       // Size of the domain in the zeta_hat direction
      double eta_top( 128.0 );          // Size of the domain in the eta direction
      const std::size_t N( 50 );       // Number of intervals in the zeta_hat direction
      const std::size_t M( 50 );       // Number of intervals in the eta direction
      const std::size_t Nvar( 4 );      // Number of variables
      double beta( 0.0 );               // Hartree parameter
      double KB( 0.0 );                 // Base flow transpiration ( +ve = blowing )
      double zeta0( 1.0 );              // Ridge/transpiration width
      double A( 0.0 );                  // Mass flux parameter
      double K( 2.5 );                  // Transpiration parameter ( +ve = blowing )
      double gamma( 20.0 );             // Steepness factor

    } // End of namespace Param

    namespace Example
    {
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

      double H( const double& hzeta )
      {
        //TODO return ridge profile as function of hzeta
        return 0.0;
      }

      double Hd( const double& hzeta )
      {
        //TODO return derivative of ridge profile wrt hzeta
        return 0.0;
      }

      double Hdd( const double& hzeta )
      {
        //TODO return 2nd derivative of ridge profile wrt hzeta
        return 0.0;
      }

      Vector<double> laplace_vals( const double& Hd, const double& Hdd, const double& Xd,
                                   const double& Yd, const double& Xdd, const double& Ydd,
                                   const double& dX, const double& dY, const double& zeta0)
      {
        // A function to return the values for the coefficients of the finite-difference
        // Laplace operator
        Vector<double> a(9,0.0);

        // X(i-1,j-1)
        a[0] = -2.*Hd*Yd*Xd / (4.*zeta0*zeta0*dY*dX);
        // X(i,j-1)
        a[1] = ( 1. + (Hd*Hd)/(zeta0*zeta0) ) * ( Yd*Yd/(dY*dY) - Ydd/(2.*dY) ) 
              + Hdd*Yd/(2.*zeta0*zeta0*dY);
        // X(i+1,j-1)
        a[2] =  2.*Hd*Yd*Xd / (4.*zeta0*zeta0*dY*dX);
        // X(i-1,j)
        a[3] = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) )/(zeta0*zeta0);
        // X(i,j)
        a[4] = -2.*( Yd*Yd*( 1. + Hd*Hd/(zeta0*zeta0) )/(dY*dY) 
                     + Xd*Xd/(zeta0*zeta0*dX*dX) );
        // X(i+1,j)
        a[5] = ( Xdd/(2.*dX) + Xd*Xd/(dX*dX) ) / (zeta0*zeta0);
        // X(i-1,j+1)
        a[6] = 2.*Hd*Yd*Xd / (4.*zeta0*zeta0*dY*dX);
        // X(i,j+1)
        a[7] = ( 1. + (Hd*Hd)/(zeta0*zeta0) ) * ( Yd*Yd/(dY*dY) + Ydd/(2.*dY) ) 
              - Hdd*Yd/(2.*zeta0*zeta0*dY);
        // X(i+1,j+1)
        a[8] = -2.*Hd*Yd*Xd / (4.*zeta0*zeta0*dY*dX);

        return a;
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
            //TODO need to check that this stuff is working
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
  if ( Param::beta > 0.1 ) { arc_step = 0.01; }

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
  cout << "Base transpiration KB = " << plate_BC.KB << endl;
  cout << "Hartree parameter beta = " << equation.beta << endl;
  cout << "UB'(eta=0) =" << base.solution()( 0, fdd ) << endl;

  cout << "We have solved the ODE problem, it is output to ./DATA/Base_soln.dat" << endl;

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
  
  
  cout << "Thetabar(0) = " << far_ode.solution()( 0, Thetabar ) << endl;

  /* ----- Solve for the perturbation quantities ----- */

  cout << "*** Solving the perturbation equations ***" << endl;

  //TODO need to create TwoD_Node_Mesh class
  // Set the current guess states  
  TwoD_node_mesh<double> Q( X_nodes, Y_nodes, 4 );
  // We use the mesh below to write the data on the original zeta-eta domain
  TwoD_node_mesh<double> Q_output( hzeta_nodes, eta_nodes, 8 );

  // Track some of the variables
  TrackerFile metric( "./DATA/A_file.dat" );
  metric.push_ptr( &Param::A, "A" );
  metric.header();
  

  // Vector for the RHS of the matrix problem
  Vector<double> B( 4 * N_eta * N_hzeta + 1, 0.0 );

  /* Iterate to a solution */
  double max_residual( 0.1 );                           // Maximum residual
  std::size_t max_iterations( 5 );                     // Maximum number of iterations
  std::size_t iteration( 0 );                           // Initialise iteration counter

  do
  { 
    // N_eta x N_hzeta mesh with 4 unknowns at each node + 1 for mass flux parameter A
    SparseMatrix<double> A( 4 * N_eta * N_hzeta + 1, 4 * N_eta * N_hzeta + 1 ); 
    cout << "Assembling sparse matrix problem" << endl; 

    using namespace Example;
    std::size_t row( 0 );                               // Initialise row counter

    /* hzeta = 0 boundary ( left boundary )*/
    std::size_t i( 0 );

    for ( std::size_t j = 0; j < Param::M + 1; ++j )
    {
      double hzeta( hzeta_nodes[ 0 ] );
      double Xd( Mesh::Xd( hzeta ) );
      double eta( eta_nodes[ j ] );
      double Yd( Mesh::Yd( eta ) );
      double Hd( Example::Hd( hzeta ) );
      Vector<double> Base( Base_soln.get_interpolated_vars( eta ) );
      // PhiB' = (2-beta)*UB - PsiB
      double PhiBd( ( 2.0 - Param::beta )*Base[ UB ] - Base[ PsiB ] );
      double UBd( Base[ UBd ] );

      // Phi_hzeta - H'( PhiB' + Phi_eta ) = 0
      if( j == 0 ) // eta = 0 ( bottom left corner )
      {
        A( row, col( i, j, Phi ) )      = -3.*Xd/(2*dX) + 3.*Hd*Yd/(2*dY); 
        A( row, col( i + 1, j, Phi ) )  =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, Phi ) )  = -1.*Xd/(2*dX);
        A( row, col( i, j + 1, Phi ) )  = -4.*Hd*Yd/(2*dY);
        A( row, col( i, j + 2, Phi ) )  =  1.*Hd*Yd/(2*dY);
        B[ row ]                        = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
                                             -Q(i+2,j,Phi) )/(2*dX) ) 
                                          + Hd*Yd*( -3*Q(i,j,Phi) + 4*Q(i,j+1,Phi)
                                             -Q(i,j+2,Phi) )/(2*dY)
                                          + Hd*PhiBd;
        ++row;
      }
      else if( j == Param::M ) // eta = eta_inf ( top left corner )
      {
        A( row, col( i, j, Phi ) )      = -3.*Xd/(2*dX) - 3.*Hd*Yd/(2*dY); 
        A( row, col( i + 1, j, Phi ) )  =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, Phi ) )  = -1.*Xd/(2*dX);
        A( row, col( i, j - 1, Phi ) )  =  4.*Hd*Yd/(2*dY);
        A( row, col( i, j - 2, Phi ) )  = -1.*Hd*Yd/(2*dY);
        B[ row ]                        = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
                                             -Q(i+2,j,Phi) )/(2*dX) )
                                          + Hd*Yd*( 3*Q(i,j,Phi) - 4*Q(i,j-1,Phi)
                                             +Q(i,j-2,Phi) )/(2*dY)
                                          + Hd*PhiBd;
        ++row;
      }
      else // Rest of the non-corner nodes
      {
        A( row, col( i, j, Phi ) )      = -3.*Xd/(2*dX);
        A( row, col( i + 1, j, Phi ) )  =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, Phi ) )  = -1.*Xd/(2*dX);
        A( row, col( i, j + 1, Phi ) )  = -Hd*Yd/(2*dY);
        A( row, col( i, j - 1, Phi ) )  =  Hd*Yd/(2*dY);
        B[ row ]                        = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
                                             -Q(i+2,j,Phi) )/(2*dX) )
                                          + Hd*Yd*( Q(i,j+1,Phi) - Q(i,j-1,Phi) )/(2*dY)
                                          + Hd*PhiBd;
        ++row;
      }

      // Psi = 0
      A( row, col( i, j, Psi ) )        =  1.;
      B[ row ]                          = -Q(i,j,Psi);
      ++row;

      // U_hzeta - H'( UB' + U_eta ) = 0
      if( j == 0 ) // eta = 0 ( bottom left corner )
      {
        A( row, col( i, j, U ) )        = -3.*Xd/(2*dX) + 3.*Hd*Yd/(2*dY); 
        A( row, col( i + 1, j, U ) )    =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, U ) )    = -1.*Xd/(2*dX);
        A( row, col( i, j + 1, U ) )    = -4.*Hd*Yd/(2*dY);
        A( row, col( i, j + 2, U ) )    =  1.*Hd*Yd/(2*dY);
        B[ row ]                        = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U) 
                                             -Q(i+2,j,U) )/(2*dX) ) 
                                          + Hd*Yd*( -3*Q(i,j,U) + 4*Q(i,j+1,U)
                                             -Q(i,j+2,U) )/(2*dY)
                                          + Hd*UBd;
        ++row;
      }
      else if( j == Param::M ) // eta = eta_inf ( top left corner )
      {
        A( row, col( i, j, U ) )        = -3.*Xd/(2*dX) - 3.*Hd*Yd/(2*dY); 
        A( row, col( i + 1, j, U ) )    =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, U ) )    = -1.*Xd/(2*dX);
        A( row, col( i, j - 1, U ) )    =  4.*Hd*Yd/(2*dY);
        A( row, col( i, j - 2, U ) )    = -1.*Hd*Yd/(2*dY);
        B[ row ]                        = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U) 
                                             -Q(i+2,j,U) )/(2*dX) )
                                          + Hd*Yd*( 3*Q(i,j,U) - 4*Q(i,j-1,U)
                                             +Q(i,j-2,U) )/(2*dY)
                                          + Hd*UBd;
        ++row;
      }
      else // Rest of the non-corner nodes
      {
        A( row, col( i, j, U ) )        = -3.*Xd/(2*dX);
        A( row, col( i + 1, j, U ) )    =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, U ) )    = -1.*Xd/(2*dX);
        A( row, col( i, j + 1, U ) )    = -Hd*Yd/(2*dY);
        A( row, col( i, j - 1, U ) )    =  Hd*Yd/(2*dY);
        B[ row ]                        = -( Xd *( -3 * Q( i, j, U ) + 4 * Q( i+1, j, U ) 
                                          -Q( i + 2, j, U ) ) / ( 2 * dX ) )
                                          + Hd*Yd*( Q(i,j+1,U) - Q(i,j-1,U) )/(2*dY)
                                          + Hd*UBd;
        ++row;
      }

      // Theta = 0      
      A( row, col( i, j, Theta ) )      =  1.;
      B[ row ]                          = -Q(i,j,Theta);
      ++row;


    } // End of loop for LHS eta nodes

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
      // Ridge profile
      double H( Example::H( hzeta ) );
      double Hd( Example::Hd( hzeta ) );
      double Hdd( Example::Hdd( hzeta ) );
      // Ridge/transpiration width
      double zeta0( Param::zeta0 );

      /* eta = 0 boundary ( bottom boundary ) */
      std::size_t j( 0 );
      double eta( eta_nodes[ j ] );
      double Yd( Mesh::Yd( eta ) );
            
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
                                          - ( 1. / ( zeta0 * zeta0 ) ) * Phi_w_hzeta;
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
        Vector<double> coeff = laplace_vals( Hd, Hdd, Xd, Yd, Xdd, Ydd, dX, dY, zeta0);
        // Guessed/known components and various derivative values
        Vector<double> Guess( Q.get_nodes_vars( i, j ) );
        Vector<double> Guess_eta( ( Q.get_nodes_vars( i, j + 1 ) 
                                  - Q.get_nodes_vars( i, j - 1 ) ) * ( Yd /( 2 * dY )) );
        Vector<double> Guess_hzeta( ( Q.get_nodes_vars( i + 1, j ) 
                                    - Q.get_nodes_vars( i - 1, j ) ) 
                                    * ( Xd /( 2 * dX )) );
        Vector<double> Guess_laplace( Q.get_nodes_vars( i - 1, j - 1 ) * coeff[0]
                                   +  Q.get_nodes_vars( i, j - 1 ) * coeff[1] 
                                   +  Q.get_nodes_vars( i + 1, j - 1 ) * coeff[2]
                                   +  Q.get_nodes_vars( i - 1, j ) * coeff[3]
                                   +  Q.get_nodes_vars( i, j ) * coeff[4]
                                   +  Q.get_nodes_vars( i + 1, j ) * coeff[5]
                                   +  Q.get_nodes_vars( i - 1, j + 1 ) * coeff[6]
                                   +  Q.get_nodes_vars( i, j + 1 ) * coeff[7]
                                   +  Q.get_nodes_vars( i + 1, j + 1 ) * coeff[8] );
        
        //////////////////
        // Phi equation //
        //////////////////

        // Laplacian of Phi        
        A( row, col( i - 1, j - 1, Phi ) )  = coeff[0];
        A( row, col( i, j - 1, Phi ) )      = coeff[1];
        A( row, col( i + 1, j - 1, Phi ) )  = coeff[2];
        A( row, col( i - 1, j, Phi ) )      = coeff[3];
        A( row, col( i, j, Phi ) )          = coeff[4];
        A( row, col( i + 1, j, Phi ) )      = coeff[5];
        A( row, col( i - 1, j + 1, Phi ) )  = coeff[6];
        A( row, col( i, j + 1, Phi ) )      = coeff[7];
        A( row, col( i + 1, j + 1, Phi ) )  = coeff[8];
        // -(2-beta)*U_eta
        A( row, col( i, j + 1, U ) )        = -( 2. - Param::beta )*Yd/( 2 * dY );
        A( row, col( i, j - 1, U ) )        =  ( 2. - Param::beta )*Yd/( 2 * dY );
        // Theta_hzeta
        A( row, col( i + 1, j, Theta ) )    =  Xd / ( 2 * dX );
        A( row, col( i - 1, j, Theta ) )    = -Xd / ( 2 * dX );
        // -H' * Theta_eta
        A( row, col( i, j + 1, Theta ) )    = -Hd*Yd / ( 2 * dY );
        A( row, col( i, j - 1, Theta ) )    =  Hd*Yd / ( 2 * dY );
        // Residual
        B[ row ]      = - Guess_laplace[ Phi ] + ( 2. - Param::beta ) * Guess_eta[ U ]
                        - Guess_hzeta[ Theta ]
                        + Hd * ( hzeta * Base[ ThetaBd ] + Guess_eta[ Theta ] )
                        + ( Hdd * PhiBd - Hd * Hd * PhiBdd )/( zeta0 * zeta0 ); 
        ++row;

        //////////////////
        // Psi equation //
        //////////////////

        // Laplacian of Psi
        A( row, col( i - 1, j - 1, Psi ) )  = coeff[0];
        A( row, col( i, j - 1, Psi ) )      = coeff[1];
        A( row, col( i + 1, j - 1, Psi ) )  = coeff[2];
        A( row, col( i - 1, j, Psi ) )      = coeff[3];
        A( row, col( i, j, Psi ) )          = coeff[4];
        A( row, col( i + 1, j, Psi ) )      = coeff[5];
        A( row, col( i - 1, j + 1, Psi ) )  = coeff[6];
        A( row, col( i, j + 1, Psi ) )      = coeff[7];
        A( row, col( i + 1, j + 1, Psi ) )  = coeff[8];

        // -(2-beta)*U_hzeta / (zeta0^2)
        A( row, col( i + 1, j, U ) )        = - ( 2. - Param::beta ) * Xd 
                                              / ( 2. * dX * zeta0 * zeta0 );
        A( row, col( i - 1, j, U ) )        =   ( 2. - Param::beta ) * Xd
                                              / ( 2. * dX * zeta0 * zeta0 );

        // (2-beta)*H'*U_eta / (zeta0^2)
        A( row, col( i, j + 1, U ) )        =   ( 2. - Param::beta ) * Hd * Yd
                                              / ( 2. * dY * zeta0 * zeta0 );
        A( row, col( i, j - 1, U ) )        = - ( 2. - Param::beta ) * Hd * Yd
                                              / ( 2. * dY * zeta0 * zeta0 );

        // -Theta_eta
        A( row, col( i, j + 1, Theta ) )    = - Yd / ( 2 * dY ); 
        A( row, col( i, j - 1, Theta ) )    =   Yd / ( 2 * dY ); 

        // Residual
        B[ row ]      = - Guess_laplace[ Psi ] + ( 2. - Param::beta ) 
                        * ( Guess_hzeta[ U ] - Hd * ( Base[ UBd ] + Guess_eta[ U ]  ) )
                        + Guess_eta[ Theta ] - Hd * Hd * hzeta * PsiBdd 
                        / ( zeta0 * zeta0 )
                        + PsiBd * ( 2. * Hd + hzeta * Hdd ) / ( zeta0 * zeta0 );

        ++row;

        ////////////////
        // U equation //
        ////////////////

        // Laplacian of U
        A( row, col( i - 1, j - 1, U ) )    = coeff[0];
        A( row, col( i, j - 1, U ) )        = coeff[1];
        A( row, col( i + 1, j - 1, U ) )    = coeff[2];
        A( row, col( i - 1, j, U ) )        = coeff[3];
        A( row, col( i, j, U ) )            = coeff[4];
        A( row, col( i + 1, j, U ) )        = coeff[5];
        A( row, col( i - 1, j + 1, U ) )    = coeff[6];
        A( row, col( i, j + 1, U ) )        = coeff[7];
        A( row, col( i + 1, j + 1, U ) )    = coeff[8];

        // -2 * beta * ( UB + UG ) * U
        A( row, col( i, j, U ) )           += - 2.* Param::beta * ( Base[ UB ] 
                                              + Guess[ U ] );

        // ( hzeta * PsiB + PsiG ) * U_hzeta
        A( row, col( i + 1, j, U ) )       +=   ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                              * Xd / ( 2 * dX ); 
        A( row, col( i - 1, j, U ) )       += - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                              * Xd / ( 2 * dX );

        // [ PhiB + PhiG - ( hzeta * PsiB + PsiG ) * H' ] * U_eta 
        A( row, col( i, j + 1, U ) )       +=   ( Base[ PhiB ] + Guess[ Phi ] - 
                                                ( hzeta * PsiB + Guess[ Psi ] ) * Hd )
                                              * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )       += - ( Base[ PhiB ] + Guess[ Phi ] - 
                                                ( hzeta * PsiB + Guess[ Psi ] ) * Hd )
                                              * Yd / ( 2 * dY );
        
        // [ UG_hzeta - H' * ( UB' + UG_eta ) ] * Psi
        A( row, col( i, j, Psi ) )          =   ( Guess_hzeta[ U ] - Hd * ( Base[ UBd ]
                                                + Guess_eta[ U ] ) );

        // ( UB' + UG_eta ) * Phi
        A( row, col( i, j, Phi ) )          =   ( Base[ UBd ] + Guess_eta[ U ] );

        // Residual
        B[ row ]      = - Guess_laplace[ U ] 
                        + Param::beta * ( 2. * Base[ UB ] + Guess[ U ] ) * Guess[ U ]
                        - ( hzeta * Base[ Psi ] + Guess[ Psi ] ) * ( Guess_hzeta[ U ]
                        - Hd * ( Base[ UBd ] + Guess_eta[ U ] ) )
                        - Base[ PhiB ] * Guess_eta[ U ]
                        - ( Base[ UBd ] + Guess_eta[ U ] ) * Guess[ Phi ]
                        + ( Hdd * Base[ UBd ] - Hd * Hd * UBdd );
 
        ++row;

        ////////////////////
        // Theta equation //
        ////////////////////

        // Laplacian of Theta
        A( row, col( i - 1, j - 1, Theta ) ) = coeff[0];
        A( row, col( i, j - 1, Theta ) )     = coeff[1];
        A( row, col( i + 1, j - 1, Theta ) ) = coeff[2];
        A( row, col( i - 1, j, Theta ) )     = coeff[3];
        A( row, col( i, j, Theta ) )         = coeff[4];
        A( row, col( i + 1, j, Theta ) )     = coeff[5];
        A( row, col( i - 1, j + 1, Theta ) ) = coeff[6];
        A( row, col( i, j + 1, Theta ) )     = coeff[7];
        A( row, col( i + 1, j + 1, Theta ) ) = coeff[8];

        // -2 * (1-beta) * (UB+UG) * [hzeta + (eta+H)*H'/(zeta0^2)] * U_eta
        A( row, col( i, j + 1, U ) )         = - 2. * ( 1. - Param::beta )
                                                 * ( Base[ UB ] + Guess[ U ] ) 
                                                 * ( hzeta + ( eta + H ) * Hd 
                                                 / ( zeta0 * zeta0 ) ) * Yd / ( 2 * dY );
        A( row, col( i, j - 1, U ) )         =   2. * ( 1. - Param::beta )
                                                 * ( Base[ UB ] + Guess[ U ] ) 
                                                 * ( hzeta + ( eta + H ) * Hd 
                                                 / ( zeta0 * zeta0 ) ) * Yd / ( 2 * dY );

        // -2 * (1-beta) * [ (hzeta - H')*(UB' + UG_eta) + UG_hzeta  ] * U
        A( row, col( i, j, U ) )             = -2. * ( 1. - Param::beta )
                                                * ( ( hzeta - Hd ) * ( Base[ UBd ] 
                                                + Guess_eta[ U ] ) + Guess_hzeta[ U ] );

        // 2 * (1-beta) * (eta + H) * (UB + UG) * U_hzeta / ( zeta0^2 )
        A( row, col( i + 1, j, U ) )         =  2. * ( 1. - Param::beta ) * ( eta + H ) 
                                                * ( Base[ UB ] + Guess[ U ] )
                                                * Xd / ( 2 * dX * zeta0 * zeta0 );
        A( row, col( i - 1, j, U ) )         = -2. * ( 1. - Param::beta ) * ( eta + H ) 
                                                * ( Base[ UB ] + Guess[ U ] )
                                                * Xd / ( 2 * dX * zeta0 * zeta0 );

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

        // - H' * ( hzeta * PsiB + PsiG ) * Theta_eta
        A( row, col( i, j + 1, Theta ) )    += -( hzeta * Base[ PsiB ] + Guess[ Psi ] ) 
                                               * Hd * Yd / ( 2 * dY );
        A( row, col( i, j - 1, Theta ) )    +=  ( hzeta * Base[ PsiB ] + Guess[ Psi ] ) 
                                               * Hd * Yd / ( 2 * dY );

        // [ThetaB + ThetaG_hzeta - H' * ( hzeta * ThetaB' + ThetaG_eta ) ] * Psi
        A( row, col( i, j, Psi ) )           =  Base[ ThetaB ] + Guess_hzeta[ Theta ]
                                               - Hd * ( hzeta * Base[ ThetaBd ]
                                               + Guess_eta[ Theta ] );

        // (2-beta) * ( UB + UG ) * Theta
        A( row, col( i, j, Theta ) )        +=   ( 2. - Param::beta ) * ( Base[ UB ] 
                                               + Guess[ U ] );

        // (2-beta) * ( hzeta * ThetaB + ThetaG ) * U
        A( row, col( i, j, U ) )            +=   ( 2. - Param::beta ) 
                                               * ( hzeta * Base[ Theta ] 
                                               + Guess[ Theta ] );

        // Residual
        B[ row ]      = - Guess_laplace[ Theta ]
                        + ( Base[ ThetaBd ] * ( 2. * Hd + hzeta * Hdd ) 
                        - hzeta * ThetaBdd * Hd * Hd  ) / ( zeta0 * zeta0 )
                        + 2.*( 1. - Param::beta ) * ( hzeta * ( Base[ UB ] + Guess[ U ] ) 
                        * Guess_eta[ U ] + hzeta * Base[ UBd ] * Guess[ U ] 
                        - ( eta + H ) * ( Base[ UB ] + Guess[ U ] ) 
                        * ( Guess_hzeta[ U ] - Hd * ( Base[ UBd ] + Guess_eta[ U ] ) ) 
                        / ( zeta0 * zeta0 ) )  
                        - ( Base[ PhiB ] + Guess[ Phi ] ) * Guess_eta[ Theta ]
                        - hzeta * Base[ ThetaBd ] * Guess[ Phi ]
                        - ( hzeta * Base[ PsiB ] + Guess[ Psi ] ) 
                        * ( Guess_hzeta[ Theta ] - Hd * ( hzeta * Base[ ThetaBd ] 
                        + Guess_eta[ U ] ) ) - Guess[ Psi ] * Base[ ThetaB]
                        - ( 2. - Param::beta ) * ( ( Base[ UB ] + Guess[ U ] ) 
                        * Guess[ Theta ] + hzeta * Base[ ThetaB ] * Guess[ U ] );
        ++row;


      }

      /* eta = eta_inf boundary ( top boundary ) */
      j = Param::M;
      eta = eta_nodes[ j ];
      Yd = Mesh::Yd( eta );

      // Phi_eta = A*(...)
      A( row, col( i, j, Phi ) )        =   3.0 *Yd/ (2*dY);
      A( row, col( i, j - 1, Phi ) )    = - 4.0 *Yd/ (2*dY);
      A( row, col( i, j - 2, Phi ) )    =   1.0 *Yd/ (2*dY);
      A( row, 4 * N_hzeta * N_eta )     = - ( zeta0 * zeta0 * hzeta * hzeta 
                                          - ( eta + H ) * ( eta + H ) ) 
                                          / pow( ( zeta0 * zeta0 * hzeta * hzeta 
                                          + ( eta + H ) * ( eta + H ) ), 2);
  
      B[ row ]        = - ( 3 * Q( i, j, Phi ) - 4 * Q( i, j-1, Phi ) 
                        + Q( i, j-2, Phi ) ) * Yd / ( 2 * dY ) 
                        + Param::A * ( zeta0 * zeta0 * hzeta * hzeta 
                        - ( eta + H ) * ( eta + H ) ) 
                        / pow( ( zeta0 * zeta0 * hzeta * hzeta 
                        + ( eta + H ) * ( eta + H ) ), 2);
      ++row;

      // Psi_eta = A*(...)
      A( row, col( i, j, Psi ) )        =   3.0 *Yd/ (2*dY);
      A( row, col( i, j - 1, Psi ) )    = - 4.0 *Yd/ (2*dY);
      A( row, col( i, j - 2, Psi ) )    =   1.0 *Yd/ (2*dY);
      A( row, 4 * N_hzeta * N_eta )     =   2. * hzeta * ( eta + H )
                                          / pow( ( zeta0 * zeta0 * hzeta * hzeta 
                                          + ( eta + H ) * ( eta + H ) ) , 2 );

      B[ row ]        = - ( 3 * Q( i, j, Psi ) - 4 * Q( i, j-1, Psi ) 
                        + Q( i, j-2, Psi ) ) * Yd / ( 2 * dY )  
                        - 2. * hzeta * ( eta + H )
                        / pow( ( zeta0 * zeta0 * hzeta * hzeta 
                        + ( eta + H ) * ( eta + H ) ) , 2 );
      ++row;

      // U = 0
      A( row, col( i, j, U ) )          =  1;
      B[ row ]                          = - Q( i, j, U );
      ++row;

      // Theta = 0
      A( row, col( i, j, Theta ) )      =  1;
      B[ row ]                          = - Q( i, j, Theta );
      ++row;

    } // End of for loop over interior nodes

    /* hzeta = hzeta_inf boundary ( right boundary ) */
    for ( std::size_t j = 0; j < Param::M + 1; ++j )
    {
      std::size_t i( Param::N );
      double hzeta( hzeta_nodes[ i ] );
      double Xd( Mesh::Xd( hzeta ) );
      double eta( eta_nodes[ j ] );
      double zeta0( Param::zeta0 ); 

      // hzeta * Phi_hzeta + 2 * Phi = A*(...)
      A( row, col( i, j, Phi ) )          =   hzeta * 3. * Xd / ( 2 * dX ) + 2.;
      A( row, col( i - 1, j, Phi ) )      = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, Phi ) )      =   hzeta * 1. * Xd / ( 2 * dX );
      A( row, 4 * N_hzeta * N_eta )       = - 2 * pow( eta, 3 ) 
                                          / pow( zeta0 * zeta0 * hzeta * hzeta 
                                          + eta * eta, 2 ); 

      B[ row ]        = - hzeta * ( 3 * Q( i, j, Phi) - 4 * Q( i - 1, j, Phi) 
                        + Q( i - 2, j, Phi) ) * Xd / ( 2 * dX ) 
                        - 2 * Q( i, j, Phi )
                        + 2 * Param::A * pow( eta, 3. ) 
                        / pow( zeta0 * zeta0 * hzeta * hzeta + eta * eta, 2 );
      ++row; 

      // hzeta * Psi_hzeta + Psi = A*(...) 
      A( row, col( i, j, Psi ) )          =   hzeta * 3. * Xd / ( 2 * dX ) + 1.;
      A( row, col( i - 1, j, Psi ) )      = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, Psi ) )      =   hzeta * 1. * Xd / ( 2 * dX );
      A( row, 4 * N_hzeta * N_eta )       = - 2. * hzeta * pow( eta, 2. )
                                          / pow( zeta0 * zeta0 * hzeta * hzeta 
                                          + eta * eta, 2 );

      B[ row ]        = - hzeta * ( 3 * Q( i, j, Psi ) - 4 * Q( i - 1, j, Psi ) 
                        + Q( i - 2, j, Psi) ) * Xd / ( 2 * dX ) 
                        - Q( i, j, Psi)
                        + 2. * Param::A * hzeta * pow( eta, 2. )
                        / pow( zeta0 * zeta0 * hzeta * hzeta + eta * eta, 2 ) ;
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
    double zeta0( Param::zeta0 );
    /* A coefficient condition */
    A( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta ) = - far_ode.solution()( 0, Thetabar );
    A( 4 * N_eta * N_hzeta, 4 * N_eta * (N_hzeta - 1) + 3 ) = zeta0 * zeta0 * hzeta;
    B[ row ] = - Q( N_hzeta-1, 0, Theta ) * zeta0 * zeta0 * hzeta 
               + Param::A * far_ode.solution()( 0, Thetabar ) ;
    
    max_residual = B.norm_inf();
    cout << "*** \t \t \t \t Maximum residual = " << max_residual << endl;

    Timer timer;
    timer.start();
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
  
    cout << "*** Iteration = " << iteration << " \t \t Maximum correction = " 
         << B.norm_inf() << endl;

    ++iteration;
  }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) ); // End iteration
  
  if ( iteration >= max_iterations )
  {
    cout << "STOPPED AFTER TOO MANY ITERATIONS" << endl;
  }

  // TODO push the data back into the unmapped domain
  
  // Update the metric file
  metric.update();


  cout << "FINISHED" << endl;
}
