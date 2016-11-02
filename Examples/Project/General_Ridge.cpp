#include <cassert>
#include <cmath>

#include "Core"

// Enumerations
enum{ f, fd, fdd, g, gd, gdd };                         // Base ODE enumeration
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };           // Base flow enumeration

enum{ Phi, Psi, U, Theta };                             // PDE enumeration

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
      const std::size_t N( 50 );       // Number of intervals in the zeta_hat direction
      const std::size_t M( 50 );       // Number of intervals in the eta direction
      const std::size_t Nvar( 4 );      // Number of variables
      double beta( 0.0 );               // Hartree parameter
      double KB( 0.0 );                 // Base flow transpiration ( +ve = blowing )
      double zeta0( 1.0 );              // Ridge/transpiration width
      double A( 0.0 );                  // Mass flux parameter
      double K( 0.0 );                  // Transpiration parameter ( +ve = blowing )
      double gamma( 20.0 );             // Steepness factor
      //TODO far-field ODE for A constraint

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
        return Param::K * 0.5 * ( 1. - tanh( Param::gamma * ( hzeta - 1. ) ) );
      }

      double Phi_w_hzeta( const double& hzeta )
      {
        // Return derivative of transpiration wrt hzeta
        double sech_squared = 1. - pow( tanh( Param::gamma * ( hzeta - 1. ) ), 2. ); 
        return - Param::K * 0.5 * Param::gamma * sech_squared;
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

      const double b1( 0.1 );
      const double b2( 0.3 );   // Y = (eta_hat + b1)^b2

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
  cout << "U'(eta=0) =" << base.solution()( 0, fdd ) << endl;

  cout << "We have solved the ODE problem, it is output to ./DATA/Base_soln.dat" << endl;

  /* ----- Solve for the perturbation quantities ----- */

  cout << "*** Solving the perturbation equations ***" << endl;

  //TODO need to create TwoD_Node_Mesh class
  // Set the current guess states  
  TwoD_node_mesh<double> Q( X_nodes, Y_nodes, 4 );
  // We use the mesh below to write the data on the original zeta-eta domain
  TwoD_node_mesh<double> Q_output( hzeta_nodes, eta_nodes, 8 );

  // Vector for the RHS of the matrix problem
  Vector<double> B( 4 * N_eta * N_hzeta + 1, 0.0 );

  /* Iterate to a solution */
  double max_residual( 0.1 );                           // Maximum residual
  std::size_t max_iterations( 1 );                     // Maximum number of iterations
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

      // Phi_zeta - H'( PhiB' + Phi_eta ) = 0
      if( j == 0 ) // eta = 0 ( bottom left corner )
      {
        A( row, col( i, j, Phi ) )      = -3.*Xd/(2*dX) + 3.*Hd*Yd/(2*dY); 
        A( row, col( i + 1, j, Phi ) )  =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, Phi ) )  = -1.*Xd/(2*dX);
        A( row, col( i, j + 1, Phi ) )  = -4.*Hd*Yd/(2*dY);
        A( row, col( i, j + 2, Phi ) )  =  1.*Hd*Yd/(2*dY);
        B[ col( i, j, Phi ) ]           = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
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
        B[ col( i, j, Phi ) ]           = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
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
        B[ col( i, j, Phi ) ]           = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
                                             -Q(i+2,j,Phi) )/(2*dX) )
                                          + Hd*Yd*( Q(i,j+1,Phi) - Q(i,j-1,Phi) )/(2*dY)
                                          + Hd*PhiBd;
        ++row;
      }

      // Psi = 0
      A( row, col( i, j, Psi ) )        =  1.;
      B[ col( i, j, Psi ) ]             = -Q(i,j,Psi);
      ++row;

      // U_zeta - H'( UB' + U_eta ) = 0
      if( j == 0 ) // eta = 0 ( bottom left corner )
      {
        A( row, col( i, j, U ) )        = -3.*Xd/(2*dX) + 3.*Hd*Yd/(2*dY); 
        A( row, col( i + 1, j, U ) )    =  4.*Xd/(2*dX);
        A( row, col( i + 2, j, U ) )    = -1.*Xd/(2*dX);
        A( row, col( i, j + 1, U ) )    = -4.*Hd*Yd/(2*dY);
        A( row, col( i, j + 2, U ) )    =  1.*Hd*Yd/(2*dY);
        B[ col( i, j, U ) ]             = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U) 
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
        B[ col( i, j, U ) ]             = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U) 
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
        B[ col( i, j, U ) ]             = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U) 
                                             -Q(i+2,j,U) )/(2*dX) )
                                          + Hd*Yd*( Q(i,j+1,U) - Q(i,j-1,U) )/(2*dY)
                                          + Hd*UBd;
        ++row;
      }

      // Theta = 0      
      A( row, col( i, j, Theta ) )      =  1.;
      B[ col( i, j, Theta ) ]           = -Q(i,j,Theta);
      ++row;


    } // End of loop for LHS eta nodes

    /* Interior points between the hzeta boundaries */
    for ( std::size_t i = 1; i < Param::N; ++i )
    {
      // hzeta location
      double hzeta( hzeta_nodes[ 0 ] );
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
      // Base solution
      Vector<double> Base( Base_soln.get_interpolated_vars( eta ) );//TODO do we need this?
            
      // Phi = Phi_w
      A( row, col( i, j, Phi ) )        =  1.;
      B[ col( i, j, Phi ) ]             = -Q(i,j,Phi) + Phi_w;
      ++row;
      // Psi = 0
      A( row, col( i, j, Psi ) )        =  1.;
      B[ col( i, j, Psi ) ]             = -Q(i,j,Psi);
      ++row;
      // U = 0
      A( row, col( i, j, U ) )          =  1.;
      B[ col( i, j, U ) ]               = -Q(i,j,U);
      ++row;
      // Theta - Psi_eta = -( 1 / ( zeta0^2 ) ) * Phi_w_hzeta
      A( row, col( i, j, Theta ) )      =  1.;
      A( row, col( i, j, Psi ) )        =  3.*Yd/(2*dY);
      A( row, col( i, j + 1, Psi ) )    = -4.*Yd/(2*dY);
      A( row, col( i, j + 2, Psi ) )    =  1.*Yd/(2*dY);
      B[ col( i, j, Theta ) ]           = -Q(i,j,Theta) + Yd*( -3*Q(i,j,Psi) 
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
        
        // Laplacian coefficients for finite-differencing
        Vector<double> coeff = laplace_vals( Hd, Hdd, Xd, Yd, Xdd, Ydd, dX, dY, zeta0);
        // Guessed/known components and various derivative values
        Vector<double> Guess( Q.get_nodes_vars( i, j ) );
        Vector<double> Guess_eta( ( Q.get_nodes_vars( i, j + 1 ) 
                                  - Q.get_nodes_vars( i, j - 1 ) ) * (Yd/(2.*dY))  );
        Vector<double> Guess_hzeta( ( Q.get_nodes_vars( i + 1, j ) 
                                    - Q.get_nodes_vars( i - 1, j ) ) * (Xd/(2.*dX)) );
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
        // -(2-beta)U_eta
        A( row, col( i, j + 1, U ) )        = -( 2. - Param::beta )*Yd/(2.*dY);
        A( row, col( i, j - 1, U ) )        =  ( 2. - Param::beta )*Yd/(2.*dY);
        // Theta_hzeta
        A( row, col( i + 1, j, Theta ) )    =  Xd / (2.*dX);
        A( row, col( i - 1, j, Theta ) )    = -Xd / (2.*dX);
        // -H' * Theta_eta
        A( row, col( i, j + 1, Theta ) )    = -Hd*Yd / (2.*dY);
        A( row, col( i, j - 1, Theta ) )    =  Hd*Yd / (2.*dY);
        // Residual
        B[ col( i, j, Phi ) ]     = -Guess_laplace[ Phi ] 
                                  + ( 2. - Param::beta ) * Guess_eta[ U ]
                                  - Guess_hzeta[ Theta ]
                                  + Hd * ( hzeta * Base[ ThetaBd ] + Guess_eta[ Theta ] )
                                  + ( Hdd * PhiBd - Hd * Hd * PhiBdd )/( zeta0 * zeta0 ); 
        ++row;

        //////////////////
        // Psi equation //
        //////////////////


        ++row;

        ////////////////
        // U equation //
        ////////////////


        ++row;

        ////////////////////
        // Theta equation //
        ////////////////////


        ++row;


      }

      /* eta = eta_inf boundary ( top boundary ) */
      j = Param::M;
      eta = eta_nodes[ j ];
      Yd = Mesh::Yd( eta );
      

    }

    /* hzeta = hzeta_inf boundary ( right boundary ) */
    i = Param::N;

    /* A coefficient condition */

    ++iteration;
  }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) ); // End iteration
  
  if ( iteration >= max_iterations )
  {
    cout << "STOPPED AFTER TOO MANY ITERATIONS" << endl;
  }
  

  cout << "FINISHED" << endl;
}
