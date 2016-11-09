#include <cassert>
#include <cmath>

#include "Core"


// PDE enumeration, this means that Phi=0, Psi=1, U=2, Theta=3
enum{ Phi, Psi, U, Theta };

// ODE enumeration -- similarly for the ODE
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB, ou, oud, ophi, otheta, othetad, opsi }; 

// either NONUNIFORM or UNIFORM meshes
#define NONUNIFORM // => X=zeta, Xd=1, Xdd=0, Y=eta, Yd=1, Ydd=0
namespace TSL
{
  namespace Example
  { 

    double K( 2.5 );	              // Wall transpiration factor ( +ve = blowing )
    double gamma( 20.0 );           // Phi_w steepness factor
    double A( 0.0 );                // Displacement coefficient
    double beta( 0.0 );             // Hartree parameter
    double KB( 0.0 );               // Base flow transpiration ( +ve = blowing )
    double hzeta_right( 16 );       // The size of the domain in the hzeta direction
    double eta_top( 128 );          // The size of the domain in the eta direction
    std::size_t N( 100 );           // Number of intervals in the hzeta direction
    std::size_t M( 100 );           // Number of intervals in the eta direction
    std::size_t Nvar( 4 );          // Number of variables
    std::string output_path("./DATA/RIDGE_TEST/"); // Data output path
    
    std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
    {
        // Return the column number for the kth variable at node (i,j)
        return Example::Nvar * ( i * ( Example::M + 1 ) + j ) + k;
    }

    double Phi_w( const double& hzeta )
    {
	      // Return a 'top-hat' function for the wall transpiration
	      return -Example::K * 0.5 * (1.0 - tanh( Example::gamma * ( hzeta - 1.0 ) ) ) ;
    }

    double Phi_w_hzeta( const double& hzeta )
    {
        // Return the derivative of the wall transpiration
        return Example::K * 0.5 * Example::gamma 
                          * std::pow( cosh( Example::gamma * ( hzeta - 1. ) ), -2 );
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
      a[1] =   ( 1. + (Hd*Hd)/(zeta0*zeta0) ) * ( Yd*Yd/(dY*dY) - Ydd/(2.*dY) ) 
             + Hdd*Yd/(2.*zeta0*zeta0*dY);
      // X(i+1,j-1)
      a[2] =  2.*Hd*Yd*Xd / (4.*zeta0*zeta0*dY*dX);
      // X(i-1,j)
      a[3] = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) )/(zeta0*zeta0);
      // X(i,j)
      a[4] = - 2.*( Yd*Yd*( 1. + Hd*Hd/(zeta0*zeta0) )/(dY*dY)
             + Xd*Xd/(zeta0*zeta0*dX*dX) );
      // X(i+1,j)
      a[5] = ( Xdd/(2.*dX) + Xd*Xd/(dX*dX) ) / (zeta0*zeta0);
      // X(i-1,j+1)
      a[6] = 2.*Hd*Yd*Xd / (4.*zeta0*zeta0*dY*dX);
      // X(i,j+1)
      a[7] =   ( 1. + (Hd*Hd)/(zeta0*zeta0) ) * ( Yd*Yd/(dY*dY) + Ydd/(2.*dY) ) 
             - Hdd*Yd/(2.*zeta0*zeta0*dY);
      // X(i+1,j+1)
      a[8] = -2.*Hd*Yd*Xd / (4.*zeta0*zeta0*dY*dX);

      return a;
    }


#ifdef UNIFORM
    // USE A UNIFORM MESH    
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
    //
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
    // USE A NON-UNIFORM MESH 
    /*const double a1( 0.1 );
    const double a2( 0.5 );   // X = (zeta + a1)^a2 

    double X( const double& zeta )
    {
        return pow(zeta + a1, a2) - pow(a1,a2);
    }
    double Xd( const double& zeta )
    {
        return a2 * pow(zeta + a1, a2 - 1);
    }
    double Xdd( const double& zeta )
    {
      return a2 * (a2 - 1) * pow(zeta + a1, a2 - 2);
    }*/

    const double gamma0( 2 );// bigger => more points around hzeta=1
    const double B0( 4 ); // bigger => more points near zeta=0 and (less for large zeta)
    //
    double X( const double& zeta )
    {
        return B0*zeta*0.5*( 1+tanh(gamma0*(1-zeta)) ) + (zeta+2*B0)*0.5*( 1+tanh(gamma0*(zeta-1)) );
    }
    double Xd( const double& zeta )
    {
        return 0.5*B0*( tanh(gamma0*(1-zeta))+1 ) - 0.5*B0*gamma0*zeta*std::pow( cosh(gamma0*(1-zeta)),-2 )
        + 0.5*( tanh(gamma0*(zeta-1))+1 ) + 0.5*gamma0*(zeta+2*B0)*std::pow( cosh(gamma0*(1-zeta)),-2 );
    }
    double Xdd( const double& zeta )
    {
      return -B0*gamma0*std::pow(cosh(gamma0*(1-zeta)),-2) 
        + gamma0*std::pow(cosh(gamma0*(zeta-1)),-2) 
        - B0*gamma0*gamma0*zeta*std::pow(cosh(gamma0*(1-zeta)),-2)*tanh(gamma0*(1 - zeta)) 
        - gamma0*gamma0*(2*B0 + zeta)*std::pow(cosh(gamma0*(zeta-1)),-2)*tanh(gamma0*(zeta-1));
    }

    //
    const double b1( 0.3 );
    const double b2( 0.3 );   // Y = (eta + b1)^b2

    double Y( const double& eta )
    {
      return std::pow(eta + b1, b2) - std::pow(b1,b2);
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

    class Base_equation : public Equation<double>
    {
      public:
        double beta;

        Base_equation() : Equation<double>( 12 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &g ) const
        {
          // CORNER BL EQNS IN LARGE ZETA LIMIT
          g[ UB ]       =  z[ UBd ];
          g[ UBd ]      =  beta * ( z[ UB ] * z[ UB ] - 1. ) - z[ PhiB ] * z[ UBd ];
          g[ PhiB ]     =  (2. - beta ) * z[ UB ] - z[ PsiB ];
          g[ ThetaB ]   =  z[ ThetaBd ];
          g[ ThetaBd ]  =  2. * ( 1. - beta ) * z[ UB ] * z[ UBd ] 
                         - z[ PhiB ] * z[ ThetaBd ] 
                         - z[ PsiB ] * z[ ThetaB ] 
                         - 2. * ( 1. - beta ) * z[ UB ] * z[ ThetaB ];
          g[ PsiB ]     =  z[ ThetaB ];
          // CORRECTION TO THE FAR-FIELD FOR THE POTENTIAL HALF LINE SINK OUTER FLOW
          g[ ou ]       =  z[ oud ];
          g[ oud ]      =  2.* beta *z[ UB ] * z[ ou ] - z[ UBd ] * z[ ophi ] 
                         - z[ PhiB ] * z[ oud ] + 2. * z[ PsiB ] * z[ ou ];
          g[ ophi ]     =  ( 2. - beta ) * z[ ou ] + z[ opsi ];
          g[ otheta ]   =  z[ othetad ];
          g[ othetad ]  =  2. * ( 1. - beta ) * ( z[ UB ] * z[ oud ] + z[ ou ] * z[ UBd ] )
                         - z[ PhiB ] * z[ othetad ] - z[ ThetaBd ] * z[ ophi ] 
                         + z[ PsiB ] * z[ otheta ] - z[ opsi ] * z[ ThetaB ] 
                         - ( 2. - beta ) * ( z[ UB ] * z[ otheta ] 
                         + z[ ou ] * z[ ThetaB ]  );
          g[ opsi ]     =  z[ otheta ];
        }
      
      
    };

    class Base_plate_BC : public Residual<double>
    {
      public:
        double KB;

        Base_plate_BC() : Residual<double> ( 6, 12 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
          B[ 0 ] = z[ UB ];
          B[ 1 ] = z[ PsiB ];
          B[ 2 ] = z[ PhiB ] + KB;
          //
          B[ 3 ] = z[ ou ];
          B[ 4 ] = z[ ophi ];
          B[ 5 ] = z[ opsi ];
        }
    };

    class Base_free_BC : public Residual<double>
    {
      public:
        double beta;

        Base_free_BC() : Residual<double> ( 6, 12 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
          B[ 0 ] = z[ UB ] - 1.;
          B[ 1 ] = z[ PsiB ] - ( 1. - beta );
          B[ 2 ] = z[ ThetaB ];
          //
          B[ 3 ] = z[ ou ];
          B[ 4 ] = z[ otheta ];
          B[ 5 ] = z[ opsi ] - 1.; 
        }
    };


  } // end Example namespace
} // end CppNoddy namespace


using namespace TSL;
using namespace std;


int main()
{
  cout << "*** ---------- Blowing Code ---------- ***" << endl;

  double zeta0( 4.0 );                                      // zeta0 is the gap size
  
  // Original coordinates are zeta & eta, then we change to X=X(zeta) & Y=Y(eta)
  double hzeta_right = Example::hzeta_right; 
  double eta_top = Example::eta_top;                      
  //
  // Define the remapped (non-uniform mesh) domain
  double left   = Example::X( 0.0 );
  double right  = Example::X( hzeta_right );
  double bottom = Example::Y(0.0);
  double top    = Example::Y( eta_top );
  //
	// Number of points 
  std::size_t N_hzeta = Example::N + 1;
  std::size_t N_X( N_hzeta );
  std::size_t N_eta = Example::M + 1;
  std::size_t N_Y( N_eta );     
  // Nodal positions in the remapped domain (spanned by X,Y)
  Vector<double> X_nodes, Y_nodes;
  X_nodes.linspace( left, right, N_X );
  Y_nodes.linspace( bottom, top, N_Y );

  // Store the original coordinates for writing data on the original zeta-eta domain
  Vector<double> eta_nodes, hzeta_nodes;
  eta_nodes.linspace( 0.0, eta_top, N_eta );
  hzeta_nodes.linspace( 0.0, hzeta_right, N_hzeta );
  //
  // To find eta=eta(Y) and zeta=zeta(X) we will use Newton iteration 
  Example::invert_eta find_eta;
  Newton<double> newton_eta( &find_eta, 50 );
  //
  for ( std::size_t j = 0; j < N_Y; ++j )
  {
    std::size_t kmin(0); double min(99);
    for ( std::size_t k = 0; k < N_Y; ++k )
    {
      if ( std::abs( Example::Y( eta_nodes[k] ) - Y_nodes[j] ) < min )
      {
        min = std::abs( Example::Y( eta_nodes[k] ) - Y_nodes[j] );
        kmin = k;
      }
    }
    find_eta.Y0 = Y_nodes[j];
    Vector<double> guess(1,1.0);
    guess[0] = eta_nodes[kmin];
    newton_eta.iterate(guess);
    eta_nodes[j] = guess[0];
  }  
  //
  Example::invert_zeta find_zeta;
  Newton<double> newton_zeta( &find_zeta, 50 );
  //
  for ( std::size_t i = 0; i < N_X; ++i )
  {
    std::size_t kmin(0); double min(99);
    for ( std::size_t k = 0; k < N_X; ++k )
    {
      if ( std::abs( Example::X( hzeta_nodes[k] ) - X_nodes[i] ) < min )
      {
        min = std::abs( Example::X( hzeta_nodes[k] ) - X_nodes[i] );
        kmin = k;
      }
    }
    find_zeta.X0 = X_nodes[i];
    Vector<double> guess(1,1.0);
    guess[0] = hzeta_nodes[kmin];
    newton_zeta.iterate(guess);
    hzeta_nodes[i] = guess[0];
  } 
  //
  // Step sizes in the remapped domain : these should be constants
  const double dY( Y_nodes[ 1 ] - Y_nodes[ 0 ] );
  const double dX( X_nodes[ 1 ] - X_nodes[ 0 ] );
  //
  // Construct and solve the ODE for the far field in the BL  
  Example::Base_equation equation;
  Example::Base_free_BC freebc;
  Example::Base_plate_BC platebc;
  equation.beta = Example::beta;
  freebc.beta   = Example::beta;
  platebc.KB    = Example::KB;
  //

  std::size_t ODE_N_eta( N_eta );                  // Number of nodes in the ODE solution
  Vector<double> ODE_eta_nodes( eta_nodes );    // ODE nodes
  cout << "# WE ARE SOLVING USING A " << N_hzeta << " x " << N_eta << " mesh, and " 
       << ODE_N_eta << " points in the ODE solver.\n";  

  ODE_BVP<double> base( &equation, ODE_eta_nodes, &platebc, &freebc );
  base.max_iterations() = 50;
  // Initial guesstimate
  for ( std::size_t j = 0; j < ODE_N_eta; ++j )
  {
    double eta( ODE_eta_nodes[ j ] );
    base.solution()( j, UB )      =   1 - exp( -eta );
    base.solution()( j, UBd )     =   exp(-eta);
    base.solution()( j, PhiB )    =   eta - 1.5 * tanh(eta);
    base.solution()( j, PsiB )    =   1 - exp( -eta );                              
    base.solution()( j, ThetaB )  = - 5 * eta * exp( -2*eta );      
    base.solution()( j, ThetaBd ) = - 5 * exp( -eta ) + 10 * eta * exp( -eta );
  }
  // solve -- probably converges to the 2D solution -- user should check!
  base.solve_bvp();
  cout << "# We have solved the ODE problem, it is output to " + Example::output_path + "base.dat.\n";
  // initial guess  
  base.solution().output( Example::output_path + "base.dat" );

  //cout << "# THE NUMBER BELOW SHOULD BE CLOSE TO ZERO for the 2D ODE solution\n";
  //cout << base.solution().integral2(UB) - base.solution().integral2(PsiB) << endl;
  cout << "# We now have a solution for the far-field on-plate flow " <<  endl;  
  cout << "# UB'(0) = " << base.solution()( 0, UBd ) << endl;
  cout << "# Thetabar(0) = " << base.solution()( 0, otheta ) << endl;
  // set the current guess states  
  TwoD_node_mesh<double> Q( X_nodes, Y_nodes, 4 );
  // we use the mesh below to write the data on the original zeta-eta domain
  TwoD_node_mesh<double> Q_output( hzeta_nodes, eta_nodes, 8 );
  
  // output a measure of the solution
  TrackerFile metric( Example::output_path + "A_file.dat" );  
  metric.push_ptr( &zeta0, "zeta0" );
  metric.push_ptr( &Example::A, "A" );
  double U_eta( 0.0 );
  double eta_half( 0.0 );
  metric.push_ptr( &U_eta, "U_eta(0,0)");
  metric.push_ptr( &eta_half, "eta at which U=1/2 on zeta=0" );
  metric.header();


  /* Vector for the RHS of the matrix problem  */
  Vector<double> B( 4 * N_eta * N_hzeta + 1, 0.0 );
  
  do                                                   // Iterate over values of zeta_0
  {
    /* Iterate to a solution */
    double max_residual( 0.0 );                        // Maximum residual
    std::size_t max_iterations( 20 );                     // Maximum number of iterations  
    std::size_t iteration( 0 );                           // Initialise iteration counter
    do 
    {
      /* N_eta x N_hzeta mesh, with 4 unknowns at each node + 1 the coefficient "A". */
      SparseMatrix<double> A( 4 * N_eta * N_hzeta + 1, 4 * N_eta * N_hzeta + 1 );
      cout << "Assembling global sparse matrix problem.\n";
      
      using namespace Example;
      std::size_t row( 0 );                               // Initialise row counter

      /* hzeta = 0 boundary ( left boundary ) */
      std::size_t i( 0 );
      
      for ( std::size_t j = 0; j < Example::M + 1 ; ++j )
      {
          double hzeta( hzeta_nodes[ 0 ] );
          double Xd( Example::Xd(hzeta) ); 
          double eta( eta_nodes[ j ] );
          double Yd( Example::Yd( eta ) );
          double Hd( Example::Hd( hzeta ) );
          Vector<double> Base( base.solution().get_interpolated_vars( eta ) );
          // PhiB' = (2-beta)*UB - PsiB
          double PhiBd( ( 2.0 - Example::beta ) * Base[ UB ] - Base[ PsiB ] );
          double UBd( Base[ UBd ] );      // TODO do we need this    
          

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
          else if( j == Example::M ) // eta = eta_inf ( top left corner )
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
                                              + Hd*Yd*( Q(i,j+1,Phi) - Q(i,j-1,Phi) )
                                              /(2*dY) + Hd*PhiBd;
            ++row;
          }

          // Psi = 0
          A( row, col( i, j, Psi ) )      =   1;
          B[ row ]                        = - ( Q( i, j, Psi ) );
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
          else if( j == Example::M ) // eta = eta_inf ( top left corner )
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
            B[ row ]                        = -( Xd *( -3 * Q( i, j, U ) 
                                              + 4 * Q( i+1, j, U ) - Q( i + 2, j, U ) ) 
                                              / ( 2 * dX ) )
                                              + Hd*Yd*( Q(i,j+1,U) - Q(i,j-1,U) )/(2*dY)
                                              + Hd*UBd;
            ++row;
          }

          // Theta = 0
          A( row, col( i, j, Theta ) )    =   1;
          B[ row ]                        = - Q( i, j, Theta );
          ++row;

      } // end for loop over LHS eta nodes

      /* Interior points between the hzeta boundaries */ 
   
      for ( std::size_t i = 1; i < Example::N ; ++i )
      {
        

        // hzeta location
        double hzeta( hzeta_nodes[ i ] );
        double Xd( Example::Xd(hzeta) );
        double Xdd( Example::Xdd(hzeta) );
        // Wall transpiration
        double Phi_w( Example::Phi_w( hzeta ) );
        double Phi_w_hzeta( Example::Phi_w_hzeta( hzeta ) );
        // Ridge profile
        double H( Example::H( hzeta ) );
        double Hd( Example::Hd( hzeta ) );
        double Hdd( Example::Hdd( hzeta ) );

        /* eta = 0 boundary ( bottom boundary ) */
        std::size_t j( 0 );						
        double eta( eta_nodes[ j ] );
        double Yd( Example::Yd(eta) );
        	
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
      
        // Main interior grid points 

        for ( std::size_t j = 1; j < Example::M ; ++j )
        {
          // eta location
          double eta( eta_nodes[ j ] );
          double Yd( Example::Yd( eta ) );
          double Ydd( Example::Ydd( eta ) );
          // Base solution
          Vector<double> Base( base.solution().get_interpolated_vars( eta ) );
          // PhiB' = (2-beta)*UB - PsiB
          double PhiBd( ( 2.0 - Example::beta ) * Base[ UB ] - Base[ PsiB ] );
          // PhiB'' = (2-beta)*UB' - PsiB' = (2-beta)*UB' - ThetaB
          double PhiBdd( ( 2.0 - Example::beta ) * Base[ UBd ] -  Base[ ThetaB ] );
          // PsiB' = Theta_B
          double PsiBd( Base[ ThetaB ] );
          // PsiB'' = ThetaB'
          double PsiBdd( Base[ ThetaBd ] ); 
          // UB'' = beta * [ UB^2 - 1] - PhiB * UB'
          double UBdd(  Example::beta * ( Base[ UB ] * Base[ UB ]  - 1. ) 
                      - Base[ PhiB ] * Base[ UBd ] );
          // ThetaB'' = 2(1-beta)*UB*UB' - PhiB*ThetaB' - PsiB*ThetaB - (2-beta)*UB*ThetaB
          double ThetaBdd( 2. * ( 1. - Example::beta ) * Base[ UB ] * Base[ UBd ]
                          - Base[ PhiB ] * Base[ ThetaBd ] - Base[ PsiB ] * Base[ ThetaB ]
                          - ( 2. - Example::beta ) * Base[ UB ] * Base[ ThetaB ] );

          // Laplacian coefficients for finite-differencing
          Vector<double> coeff = laplace_vals( Hd, Hdd, Xd, Yd, Xdd, Ydd, dX, dY, zeta0);

          // Guessed/known components and various derivative values      
          Vector<double> Guess( Q.get_nodes_vars(i,j) );
          Vector<double> Guess_Y( ( Q.get_nodes_vars(i,j+1) 
                                  - Q.get_nodes_vars(i,j-1) ) / (2*dY) );
          Vector<double> Guess_X( ( Q.get_nodes_vars(i+1,j) 
                                  - Q.get_nodes_vars(i-1,j) ) / (2*dX) );
          
          Vector<double> Guess_eta( ( Q.get_nodes_vars( i, j + 1 ) 
                                    - Q.get_nodes_vars( i, j - 1 ) ) * ( Yd /( 2*dY )) );
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
          A( row, col( i, j + 1, U ) )        = -( 2. - Example::beta )*Yd/( 2 * dY );
          A( row, col( i, j - 1, U ) )        =  ( 2. - Example::beta )*Yd/( 2 * dY );
          // Theta_hzeta
          A( row, col( i + 1, j, Theta ) )    =  Xd / ( 2 * dX );
          A( row, col( i - 1, j, Theta ) )    = -Xd / ( 2 * dX );
          // -H' * Theta_eta
          A( row, col( i, j + 1, Theta ) )    = -Hd*Yd / ( 2 * dY );
          A( row, col( i, j - 1, Theta ) )    =  Hd*Yd / ( 2 * dY );
          // Residual
          B[ row ]      = - Guess_laplace[ Phi ] + ( 2. - Example::beta ) * Guess_eta[ U ]
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
          A( row, col( i + 1, j, U ) )        = - ( 2. - Example::beta ) * Xd 
                                                / ( 2. * dX * zeta0 * zeta0 );
          A( row, col( i - 1, j, U ) )        =   ( 2. - Example::beta ) * Xd
                                                / ( 2. * dX * zeta0 * zeta0 );

          // (2-beta)*H'*U_eta / (zeta0^2)
          A( row, col( i, j + 1, U ) )        =   ( 2. - Example::beta ) * Hd * Yd
                                                / ( 2. * dY * zeta0 * zeta0 );
          A( row, col( i, j - 1, U ) )        = - ( 2. - Example::beta ) * Hd * Yd
                                                / ( 2. * dY * zeta0 * zeta0 );

          // -Theta_eta
          A( row, col( i, j + 1, Theta ) )    = - Yd / ( 2 * dY ); 
          A( row, col( i, j - 1, Theta ) )    =   Yd / ( 2 * dY ); 

          // Residual
          B[ row ]      = - Guess_laplace[ Psi ] + ( 2. - Example::beta ) 
                          * ( Guess_hzeta[ U ] - Hd * ( Base[ UBd ] + Guess_eta[ U ]  ) )
                          / ( zeta0 * zeta0 )
                          + Guess_eta[ Theta ] 
                          - Hd * Hd * hzeta * PsiBdd / ( zeta0 * zeta0 )  
                          + PsiBd * ( 2. * Hd + hzeta * Hdd ) / ( zeta0 * zeta0 );

          ++row;

          ////////////////    
          // U equation //
          ////////////////

          // u_i,j+1
          A( row, col( i, j + 1, U ) )      =   Yd*Yd/(dY*dY) + Ydd/(2*dY) + (Guess[Phi] 
                                              + Base[PhiB])* Yd / (2*dY);
          // u_i,j
          A( row, col( i, j, U ) )          = - 2*Yd*Yd/(dY*dY) - 2*Xd*Xd
                                              / ( zeta0*zeta0*dX*dX );
          // u_i,j-1
          A( row, col( i, j - 1, U ) )      =   Yd*Yd/(dY*dY) - Ydd/(2*dY) - (Guess[Phi] 
                                              + Base[PhiB])* Yd / (2*dY);
          // u_i+1,j
          A( row, col( i + 1, j, U ) )      =   Xd*Xd/(zeta0*zeta0*dX*dX) + Xdd
                                              / ( 2 * zeta0 * zeta0 * dX ) 
                                              + (Guess[Psi] + Base[PsiB]*hzeta) * Xd 
                                              / ( 2 * dX );
          // u_i-1,j
          A( row, col( i - 1, j, U ) )      =   Xd*Xd/(zeta0*zeta0*dX*dX) - Xdd
                                              / ( 2*zeta0*zeta0*dX ) 
                                              - (Guess[Psi] + Base[PsiB]*hzeta) * Xd 
                                              / ( 2*dX );
          // phi_i,j
          A( row, col( i, j, Phi ) )        =   ( Yd * Guess_Y[U] + Base[UBd] );
          // psi_i,j
          A( row, col( i, j, Psi ) )        =   ( Xd * Guess_X[U] + 0);
          // RHS
          B[ row ]                          = - Guess_laplace[ U ] 
                                              - Yd * (Base[PhiB] + Guess[Phi]) 
                                              * Guess_Y[ U ] - Xd 
                                              * (hzeta*Base[PsiB] + Guess[Psi]) 
                                              * Guess_X[ U ] - Base[UBd] * Guess[Phi];
          ++row;

          ////////////////////
          // Theta equation //
          ////////////////////

          // theta_i,j+1
          A( row, col( i, j + 1, Theta ) )  =   Yd*Yd/(dY*dY) + Ydd/(2*dY) 
                                              + (Guess[Phi]+Base[PhiB]) * Yd / (2*dY);
          // theta_i,j  
          A( row, col( i, j, Theta ) )      = - 2 * Yd*Yd / (dY*dY) - 2*Xd*Xd 
                                              / (zeta0*zeta0*dX*dX) 
                                              +  2 * (Guess[U]+Base[UB]);
          // theta_i,j-1
          A( row, col( i, j - 1, Theta ) )  =   Yd*Yd/(dY*dY) - Ydd/(2*dY) 
                                              - (Guess[Phi]+Base[PhiB]) * Yd / (2*dY);
          // theta_i+1,j
          A( row, col( i + 1, j, Theta ) )  =   Xd*Xd/(zeta0*zeta0*dX*dX) + Xdd
                                              / ( 2*zeta0*zeta0*dX ) 
                                              + (Guess[Psi]+hzeta*Base[PsiB]) * Xd 
                                              / (2*dX);
          // theta_i-1,j
          A( row, col( i - 1, j, Theta ) )  =   Xd*Xd/(zeta0*zeta0*dX*dX) - Xdd
                                              / ( 2*zeta0*zeta0*dX ) 
                                              - (Guess[Psi]+hzeta*Base[PsiB]) * Xd 
                                              / (2*dX);
          // u_i,j+1
          A( row, col( i, j + 1, U ) )      = - hzeta*(Guess[U]+Base[UB])*Yd / dY;
          // u_i,j
          A( row, col( i, j, U ) )          = - 2*hzeta*(Yd*Guess_Y[U]+Base[UBd]) 
                                              + 2*eta*Xd*Guess_X[ U ]/(zeta0*zeta0) 
                                              + 2*(Guess[Theta]+hzeta*Base[ThetaB]) ;
          // u_i,j-1
          A( row, col( i, j - 1, U ) )      =   hzeta*(Guess[U]+Base[UB])*Yd / dY;
          // u_i+1,j
          A( row, col( i + 1, j, U ) )      =   eta*(Guess[ U ]+Base[UB])*Xd
                                              / ( zeta0*zeta0*dX );
          // u_i-1,j
          A( row, col( i - 1, j, U ) )      = - eta*(Guess[ U ]+Base[UB])*Xd
                                              / ( zeta0*zeta0*dX );
          // phi_i,j
          A( row, col( i, j, Phi ) )        =   ( Yd*Guess_Y[ Theta ] 
                                              + hzeta*Base[ThetaBd] );
          // psi_i,j
          A( row, col( i, j, Psi ) )        =   ( Xd*Guess_X[ Theta ] + Base[ThetaB] );
          // RHS
          B[ row ]    = - Guess_laplace[ Theta ] - Yd * (Base[PhiB] + Guess[Phi]) 
                        * Guess_Y[Theta] - Xd * (hzeta * Base[PsiB] + Guess[Psi]) 
                        * Guess_X[Theta] + 2 * hzeta * Yd * (Base[UB] + Guess[U]) 
                        * Guess_Y[U] - 2 * (eta / (zeta0 * zeta0)) * Xd 
                        * (Base[UB] + Guess[U]) * Guess_X[U] - 2 * (Guess[Theta] 
                        + hzeta * (Base[ThetaB] - Base[UBd])) * Guess[U]
					              - hzeta * Base[ThetaBd] * Guess[Phi] - Base[ThetaB] * Guess[Psi]
					              - 2 * Base[UB] * Guess[Theta]; 
          ++row;
        }

        /* eta = eta_inf boundary ( top boundary ) */
        j = Example::M ;
        eta = eta_nodes[ j ];
        Yd = Example::Yd(eta);

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
                          + Example::A * ( zeta0 * zeta0 * hzeta * hzeta 
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
                          - Example::A * 2. * hzeta * ( eta + H )
                          / pow( ( zeta0 * zeta0 * hzeta * hzeta 
                          + ( eta + H ) * ( eta + H ) ) , 2 );
        ++row;
          
        // U = 0
        A( row, col( i, j, U ) )            =   1;
        B[ row ]                            = - ( Q( i, j, U ) );
        ++row;
        
        // Theta = 0
        A( row, col( i, j, Theta ) )        =   1;
        B[ row ]                            = - ( Q(i,j,Theta) );
        ++row;

      } // end of for loop over interior nodes
        
      // hzeta = hzeta_inf ( right boundary )
      for ( std::size_t j = 0; j < Example::M + 1; ++j )
      {
        //offset for global problem
        std::size_t i( N_hzeta-1 );
        double hzeta( hzeta_nodes[ i ] );
        double Xd( Example::Xd(hzeta) );
        double eta( eta_nodes[j] );
        
        double zeta( zeta0*hzeta );
        
        Vector<double> Base( base.solution().get_interpolated_vars( eta ) );

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
                          + 2 * Example::A * pow( eta, 3. ) 
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
                          + 2. * Example::A * hzeta * pow( eta, 2. )
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
        
      // WE ALSO NEED A FINAL CLOSING CONDITION FOR THE DISPLACEMENT COEFFICIENT "A"
      
      // zeta0^2 * hzeta_max theta( zeta=zeta_max, eta=0 ) = A*otheta(0)
      Vector<double> Base( base.solution().get_interpolated_vars( 0.0 ) );
      A( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta ) = -Base[otheta];
      A( 4 * N_eta * N_hzeta, 4 * N_eta * (N_hzeta - 1) + 3 ) =  zeta0*zeta0*hzeta_right;
      // RHS
      B[ row ] = - Q( N_hzeta-1, 0, Theta ) * zeta0 * zeta0 * hzeta_right 
                 + Example::A * Base[ otheta ];   
    
      max_residual = B.norm_inf();
      cout << "***                                              Maximum residual = " 
           << B.norm_inf() << "\n";  

      Timer timer;
      timer.start();
      Vector<double> x;
      x = A.solve( B );
      B = x;      
      timer.print();  
      timer.stop();

      /* Update the known values using the corrections which we just found */  
        
      for ( std::size_t i = 0; i <= Example::N ; ++i )
      {
        for ( std::size_t j = 0; j <= Example::M ; ++j )
        {
          Q( i, j, Phi )    += B[ col( i, j, Phi ) ];
          Q( i, j, Psi )    += B[ col( i, j, Psi ) ];
          Q( i, j, U )      += B[ col( i, j, U ) ];
          Q( i, j, Theta )  += B[ col( i, j, Theta ) ];
        }
      }
      Example::A += B[ 4 * ( Example::N + 1 ) * ( Example::M + 1 ) ];

      cout << "***    Iteration = " << iteration 
           << "    Maximum correction = " << B.norm_inf() << "\n";  
      ++iteration;
    } while ( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) );
    
    /* End of solution iteration */

    if ( iteration >= max_iterations )
    {
      cout << "STOPPED AFTER TOO MANY ITERATIONS \n";   
      assert(false);
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
                              + base.solution().get_interpolated_vars( eta )[PhiB];
        Q_output( i, j, 5 ) = Q( i, j, Psi)   
                   + hzeta * base.solution().get_interpolated_vars( eta )[PsiB];
        Q_output( i, j, 6 ) =   Q( i, j, U)
                              + base.solution().get_interpolated_vars( eta )[UB];
        Q_output( i, j, 7 ) = Q( i, j, Theta) 
                   + hzeta * base.solution().get_interpolated_vars( eta )[ThetaB];
      }
    }
    // Convert to string //TODO we need a utility function to do this
    std::stringstream ss;
    ss << zeta0;
    std::string zeta0_str = ss.str(); 

    Q_output.dump_gnu( Example::output_path + "Qout_" + zeta0_str + ".dat" );

    Vector<double> Base( base.solution().get_interpolated_vars( 0.0 ) ); 
    U_eta = -( 3*Q_output(0,0,U+4) - 4*Q_output(0,1,U+4) 
            + Q_output(0,2,U+4) )*Example::Yd(0.0)/(2*dY);
  
    // Find value of eta on zeta=0 at which U=1/2
    std::size_t lower = 0;
    std::size_t upper = 1;
    for (std::size_t j=0; j < Example::M; ++j)
    {
    if ( Q_output(0,j,U+4) < 0.5 && Q_output(0,j+1,U+4) > 0.5 ) { lower = j; upper=j+1; } 
    }
    // linearly interpolate
    eta_half =  ( 0.5 - Q_output(0,lower,U+4) ) * ( eta_nodes[upper] - eta_nodes[lower] ) 
              / ( Q_output(0,upper,U+4) - Q_output(0,lower,U+4)  ) + eta_nodes[lower];

    metric.update();

    cout << " zeta0 = " << zeta0 << ", A = " << Example::A << endl;
    zeta0 += 0.5; // 0.5 is best for blowing

  }while(zeta0 < 4.5);


  cout << " FINISHED SOLUTION \n";

}
