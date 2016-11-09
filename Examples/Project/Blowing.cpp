#include <cassert>
#include <cmath>

#include "Core"


// PDE enumeration, this means that Phi=0, Psi=1, U=2, Theta=3
enum{ Phi, Psi, U, Theta };

// ODE enumeration -- similarly for the ODE
enum{ oU, oUd, oPhi, oTheta, oThetad, oPsi, ou, oud, ophi, otheta, othetad, opsi }; 

// either NONUNIFORM or UNIFORM meshes
#define NONUNIFORM // => X=zeta, Xd=1, Xdd=0, Y=eta, Yd=1, Ydd=0
namespace TSL
{
  namespace Example
  { 

    double C_w( -2.5 );	               // Wall transpiration factor
    double gamma( 20.0 );              // Phi_w steepness factor
    double A( 0.0 );                   // Displacement coefficient
    double zeta_hat_right( 16 );       // The size of the domain in the zeta_hat direction
    double eta_top( 128 );                 // The size of the domain in the eta direction
    unsigned N( 100 );                     // Number of intervals in the zeta_hat direction
    unsigned M( 100 );                     // Number of intervals in the eta direction
    unsigned Nvar( 4 );                    // Number of variables
    std::string output_path("./DATA/BLOWING_Cw_2.5_101x101_16_128/"); // Data output path
    
    unsigned col( const unsigned& i, const unsigned& j, const unsigned& k )           
    {        
        // Return the column number for kth variable at node (i,j)
        return Example::Nvar * ( i * (Example::M + 1) + j ) + k;   
    }

    double Phi_w( const double& zeta_hat )
    {
	      // Return a 'top-hat' function for the wall transpiration
	      return Example::C_w * 0.5 * (1.0 - tanh( Example::gamma * ( zeta_hat - 1.0 ) ) ) ;
    }

    double Phi_w_zeta_hat( const double& zeta_hat )
    {
        // Return the derivative of the wall transpiration
        return -Example::C_w * 0.5 * Example::gamma * std::pow( cosh( Example::gamma * ( zeta_hat - 1. ) ), -2 );
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

    const double gamma0( 2 );// bigger => more points around zeta_hat=1
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

    class Farfield_equation : public Equation<double>
    {
    public:

      Farfield_equation() : Equation<double>( 12 ) {}

      void residual_fn( const Vector<double> &z, Vector<double> &g ) const
      {
        // CORNER BL EQNS IN LARGE ZETA LIMIT
        g[ oU ] = z[ oUd ];
        g[ oUd ] = - z[ oPhi ] * z[ oUd ];
        //
        g[ oPhi ] = 2*z[ oU ] - z[ oPsi ];
        //
        g[ oTheta ] = z[ oThetad ];
        g[ oThetad ] = 2 * z[ oU ] * z[ oUd ] - z[ oPhi ] * z[ oThetad ] 
          - z[ oPsi ] * z[ oTheta ] - 2 * z[ oU ] * z[ oTheta ];
        //
        g[ oPsi ] = z[ oTheta ];
        // CORRECTION TO THE FAR-FIELD FOR THE POTENTIAL HALF LINE SINK OUTER FLOW
        g[ ou ] = z[ oud ];
        g[ oud ] = -z[oPhi]*z[oud]-z[ophi]*z[oUd]+2*z[oPsi]*z[ou];
        //
        g[ ophi ] = 2*z[ou]+z[opsi];
        //
        g[ otheta ] = z[othetad];
        g[ othetad ] = 2*z[oU]*z[oud]+2*z[ou]*z[oUd]-z[oPhi]*z[othetad]-z[ophi]*z[oThetad]+z[oPsi]*z[otheta]
             - z[opsi]*z[oTheta] - 2*z[oU]*z[otheta] - 2*z[ou]*z[oTheta];
        //
        g[ opsi ] = z[otheta];
      }
      
      /*void matrix0( const DenseVector<double> &x, DenseMatrix<double> &m ) const
      {
        Utility::fill_identity(m);
      }*/
    };

    class Farfield_plate_BC : public Residual<double>
    {
    public:
      Farfield_plate_BC() : Residual<double> ( 6, 12 ) {}

      void residual_fn( const Vector<double> &z, Vector<double> &B ) const
      {
        B[ 0 ] = z[ oU ];
        B[ 1 ] = z[ oPsi ];
        B[ 2 ] = z[ oPhi ];
        //
        B[ 3 ] = z[ ou ];
        B[ 4 ] = z[ ophi ];
        B[ 5 ] = z[ opsi ];
      }
    };

    class Farfield_free_BC : public Residual<double>
    {
    public:
      Farfield_free_BC() : Residual<double> ( 6, 12 ) {}

      void residual_fn( const Vector<double> &z, Vector<double> &B ) const
      {
        B[ 0 ] = z[ oU ] - 1.0;
        B[ 1 ] = z[ oPsi ] - 1.0;
        B[ 2 ] = z[ oTheta ];
        //
        B[ 3 ] = z[ ou ];
        B[ 4 ] = z[ otheta ];
        B[ 5 ] = z[ opsi ] - 1; 
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
  double zeta_hat_right = Example::zeta_hat_right; 
  double eta_top = Example::eta_top;                      
  //
  // Define the remapped (non-uniform mesh) domain
  double left   = Example::X( 0.0 );
  double right  = Example::X( zeta_hat_right );
  double bottom = Example::Y(0.0);
  double top    = Example::Y( eta_top );
  //
	// Number of points 
  unsigned N_zeta = Example::N + 1;
  unsigned N_X( N_zeta );
  unsigned N_eta = Example::M + 1;
  unsigned N_Y( N_eta );     
  // Nodal positions in the remapped domain (spanned by X,Y)
  Vector<double> X_nodes, Y_nodes;
  X_nodes.linspace( left, right, N_X );
  Y_nodes.linspace( bottom, top, N_Y );

  // Store the original coordinates for writing data on the original zeta-eta domain
  Vector<double> eta_nodes, zeta_hat_nodes;
  eta_nodes.linspace( 0.0, eta_top, N_eta );
  zeta_hat_nodes.linspace( 0.0, zeta_hat_right, N_zeta );
  //
  // To find eta=eta(Y) and zeta=zeta(X) we will use Newton iteration 
  Example::invert_eta find_eta;
  Newton<double> newton_eta( &find_eta, 50 );
  //
  for ( unsigned j = 0; j < N_Y; ++j )
  {
    unsigned kmin(0); double min(99);
    for ( unsigned k = 0; k < N_Y; ++k )
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
  for ( unsigned i = 0; i < N_X; ++i )
  {
    unsigned kmin(0); double min(99);
    for ( unsigned k = 0; k < N_X; ++k )
    {
      if ( std::abs( Example::X( zeta_hat_nodes[k] ) - X_nodes[i] ) < min )
      {
        min = std::abs( Example::X( zeta_hat_nodes[k] ) - X_nodes[i] );
        kmin = k;
      }
    }
    find_zeta.X0 = X_nodes[i];
    Vector<double> guess(1,1.0);
    guess[0] = zeta_hat_nodes[kmin];
    newton_zeta.iterate(guess);
    zeta_hat_nodes[i] = guess[0];
  } 
  //
  // Step sizes in the remapped domain : these should be constants
  const double dY( Y_nodes[ 1 ] - Y_nodes[ 0 ] );
  const double dX( X_nodes[ 1 ] - X_nodes[ 0 ] );
  //
  // Construct and solve the ODE for the far field in the BL  
  Example::Farfield_equation equation;
  Example::Farfield_free_BC freebc;
  Example::Farfield_plate_BC platebc;
  //

  unsigned ODE_N_eta( N_eta );                  // Number of nodes in the ODE solution
  Vector<double> ODE_eta_nodes( eta_nodes );    // ODE nodes
  cout << "# WE ARE SOLVING USING A " << N_zeta << " x " << N_eta << " mesh, and " 
       << ODE_N_eta << " points in the ODE solver.\n";  

  ODE_BVP<double> farfield( &equation, ODE_eta_nodes, &platebc, &freebc );
  farfield.max_iterations() = 50;
  // Initial guesstimate
  for ( unsigned j = 0; j < ODE_N_eta; ++j )
  {
    double eta( ODE_eta_nodes[ j ] );
    farfield.solution()( j, oU )      =   1 - exp( -eta );
    farfield.solution()( j, oUd )     =   exp(-eta);
    farfield.solution()( j, oPhi )    =   eta - 1.5 * tanh(eta);
    farfield.solution()( j, oPsi )    =   1 - exp( -eta );                              
    farfield.solution()( j, oTheta )  = - 5 * eta * exp( -2*eta );      
    farfield.solution()( j, oThetad ) = - 5 * exp( -eta ) + 10 * eta * exp( -eta );
  }
  // solve -- probably converges to the 2D solution -- user should check!
  farfield.solve_bvp();
  cout << "# We have solved the ODE problem, it is output to " + Example::output_path + "farfield.dat.\n";
  // initial guess  
  farfield.solution().output( Example::output_path + "farfield.dat" );

  //cout << "# THE NUMBER BELOW SHOULD BE CLOSE TO ZERO for the 2D ODE solution\n";
  //cout << farfield.solution().integral2(oU) - farfield.solution().integral2(oPsi) << endl;
  cout << "# We now have a solution for the far-field on-plate flow " <<  endl;  
  cout << "# UB'(0) = " << farfield.solution()( 0, oUd ) << endl;
  cout << "# Thetabar(0) = " << farfield.solution()( 0, otheta ) << endl;
  // set the current guess states  
  TwoD_node_mesh<double> Q( X_nodes, Y_nodes, 4 );
  // we use the mesh below to write the data on the original zeta-eta domain
  TwoD_node_mesh<double> Q_output( zeta_hat_nodes, eta_nodes, 8 );
  
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
  Vector<double> B( 4 * N_eta * N_zeta + 1, 0.0 );
  
  do                                                   // Iterate over values of zeta_0
  {
    /* Iterate to a solution */
    double max_residual( 0.0 );                        // Maximum residual
    unsigned max_iterations( 20 );                     // Maximum number of iterations  
    unsigned iteration( 0 );                           // Initialise iteration counter
    do 
    {
      /* N_eta x N_zeta mesh, with 4 unknowns at each node + 1 the coefficient "A". */
      SparseMatrix<double> A( 4 * N_eta * N_zeta + 1, 4 * N_eta * N_zeta + 1 );
      cout << "Assembling global sparse matrix problem.\n";
      
      using namespace Example;
      unsigned row( 0 );                               // Initialise row counter

      /* zeta_hat = 0 boundary ( left boundary ) */
      unsigned i( 0 );
      
      for ( unsigned j = 0; j < Example::M + 1 ; ++j )
      {
          double zeta_hat( zeta_hat_nodes[ 0 ] );
          double Xd( Example::Xd(zeta_hat) );           
          
          // phi_zeta = 0
          A( row, col( i, j, Phi ) )      = - 3.*Xd/(2*dX);
          A( row, col( i + 1, j, Phi ) )  =   4.*Xd/(2*dX);
          A( row, col( i + 2, j, Phi ) )  = - 1.*Xd/(2*dX);
          B[ row ]                        = - ( ( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi) 
                                            - Q(i+2,j,Phi) )*Xd/(2*dX) );
          ++row;
          // psi = 0
          A( row, col( i, j, Psi ) )      =   1;
          B[ row ]                        = - ( Q( i, j, Psi ) );
          ++row;
          // u_zeta = 0
          A( row, col( i, j, U ) )        = - 3.*Xd/(2*dX);
          A( row, col( i + 1, j, U ) )    =   4.*Xd/(2*dX);
          A( row, col( i + 2, j, U ) )    = - 1.*Xd/(2*dX);
          B[ row ]                        = - ( ( -3*Q(i,j,U) + 4*Q(i+1,j,U) 
                                            - Q(i+2,j,U) )*Xd/(2*dX) );
          ++row;
          //// theta = 0
          A( row, col( i, j, Theta ) )    =   1;
          B[ row ]                        = - Q( i, j, Theta );
          ++row;
      } // end for loop over LHS eta nodes

      /* Interior points between the zeta_hat boundaries */ 
   
      for ( unsigned i = 1; i < Example::N ; ++i )
      {
        /* eta = 0 boundary ( bottom boundary ) */
        unsigned j( 0 );						
        // zeta location
        double zeta_hat( zeta_hat_nodes[ i ] );
        double Xd( Example::Xd(zeta_hat) );
        double Xdd( Example::Xdd(zeta_hat) );
        // eta location
        double eta( eta_nodes[ j ] );
        double Yd( Example::Yd(eta) );
        // Wall transpiration
        double Phi_w( Example::Phi_w( zeta_hat ) );
        double Phi_w_zeta_hat( Example::Phi_w_zeta_hat( zeta_hat ) );
        // Get base flow solution from the ODE
        Vector<double> base( farfield.solution().get_interpolated_vars( eta ) );
        //
        // Boundary conditions	
        // Phi = Phi_w
        A( row, col( i, j, Phi ) )        =   1;
        B[ row ]                          = - Q( i, j, Phi ) + Phi_w;
        ++row;
        // Psi = 0
        A( row, col( i, j, Psi ) )        =   1;
        B[ row ]                          = - Q( i, j, Psi );
        ++row;
        // U = 0
        A( row, col( i, j, U ) )          =   1;
        B[ row ]                          = - Q( i, j, U );
        ++row;
        // Theta - Psi_eta = -(1/zeta0^2)*Phi_w_zeta_hat
        A( row, col( i, j, Theta ) )      =   1;
        A( row, col( i, j, Psi ) )        =   Yd*3./(2*dY);
        A( row, col( i, j + 1, Psi ) )    = - Yd*4./(2*dY);
        A( row, col( i, j + 2, Psi ) )    =   Yd/(2*dY);
        B[ row ]                          = - Q( i, j, Theta ) + Yd*( -3*Q(i, j, Psi) 
                                            + 4 * Q(i, j + 1, Psi) - Q(i, j + 2, Psi) )
                                            / (2*dY) - Phi_w_zeta_hat / (zeta0 * zeta0);
        ++row;     
      
        // Main interior grid points 

        for ( unsigned j = 1; j < Example::M ; ++j )
        {
          double eta( eta_nodes[ j ] );
          double Yd( Example::Yd(eta) );
          double Ydd( Example::Ydd(eta) );
          /* Get the underlying 2D base flow solution for this value of eta */
          Vector<double> base( farfield.solution().get_interpolated_vars( eta ) );
          /* Get the guessed/known components and find various derivative values */      
          Vector<double> Guess( Q.get_nodes_vars(i,j) );
          Vector<double> Guess_Y( ( Q.get_nodes_vars(i,j+1) 
                                  - Q.get_nodes_vars(i,j-1) ) / (2*dY) );
          Vector<double> Guess_X( ( Q.get_nodes_vars(i+1,j) 
                                  - Q.get_nodes_vars(i-1,j) ) / (2*dX) );
          Vector<double> Guess_laplace( ( Q.get_nodes_vars(i+1,j) 
            - Q.get_nodes_vars(i,j)*2 + Q.get_nodes_vars(i-1,j) ) * Xd * Xd 
            / (zeta0*zeta0*dX*dX)
            + ( Q.get_nodes_vars(i+1,j) - Q.get_nodes_vars(i-1,j) ) * Xdd 
            / (zeta0*zeta0*2*dX) 
            + ( Q.get_nodes_vars(i,j+1) - Q.get_nodes_vars(i,j)*2 
            +   Q.get_nodes_vars(i,j-1) ) * Yd * Yd / (dY*dY) 
            + ( Q.get_nodes_vars(i,j+1) - Q.get_nodes_vars(i,j-1) ) * Ydd / (2*dY) );
          //

          //////////////////
          // Phi equation //
          //////////////////
          
          // phi_i,j+1
          A( row, col( i, j + 1, Phi ) )    =   Yd*Yd/(dY*dY) + Ydd/(2*dY);       
          // phi_i,j
          A( row, col( i, j, Phi ) )        = - 2*Yd*Yd/(dY*dY) - 2*Xd*Xd
                                              / ( zeta0 * zeta0 * dX * dX );
          // phi_i,j-1
          A( row, col( i, j - 1, Phi ) )    =   Yd*Yd/(dY*dY) - Ydd/(2*dY);
          // phi_i+1,j
          A( row, col( i + 1, j, Phi ) )    =   Xd*Xd/(zeta0*zeta0*dX*dX) + Xdd
                                              / ( 2*zeta0*zeta0*dX );
          // phi_i-1,j
          A( row, col( i - 1, j, Phi ) )    =   Xd*Xd/(zeta0*zeta0*dX*dX) - Xdd
                                              / ( 2*zeta0*zeta0*dX );
          // u_i,j+1
          A( row, col( i, j + 1, U ) )      = - Yd/dY;
          // u_i,j-1
          A( row, col( i, j - 1, U ) )      =   Yd/dY;
          // theta_i+1,j
          A( row, col( i + 1, j, Theta ) )  =   Xd/(2*dX);
          // theta_i-1,j
          A( row, col( i - 1, j, Theta ) )  = - Xd/(2*dX);
          // RHS
          B[ row ]                          = - Guess_laplace[ Phi ] 
                                              + 2 * Yd * Guess_Y[ U ] 
                                              - Xd * Guess_X[ Theta ]; 
          ++row;

          //////////////////
          // Psi equation //
          //////////////////
          
          // psi_i,j+1
          A( row, col( i, j + 1, Psi ) )    =   Yd*Yd/(dY*dY) + Ydd/(2*dY);        
          // psi_i,j
          A( row, col( i, j, Psi ) )        = - 2*Yd*Yd/(dY*dY) - 2*Xd*Xd
                                              / ( zeta0*zeta0*dX*dX );
          // psi_i,j-1
          A( row, col( i, j - 1, Psi ) )    =   Yd*Yd/(dY*dY) - Ydd/(2*dY);
          // psi_i+1,j
          A( row, col( i + 1, j, Psi ) )    =   Xd*Xd/(zeta0*zeta0*dX*dX) + Xdd
                                              / ( 2*zeta0*zeta0*dX );
          // psi_i-1,j  
          A( row, col( i - 1, j, Psi ) )    =   Xd*Xd/(zeta0*zeta0*dX*dX) - Xdd
                                              / ( 2*zeta0*zeta0*dX );
          // u_i+1,j
          A( row, col( i + 1, j, U ) )      = - Xd/(zeta0*zeta0*dX);        
          // u_i-1,j  
          A( row, col( i - 1, j, U ) )      =   Xd/(zeta0*zeta0*dX);        
          // theta_i,j+1
          A( row, col( i, j + 1, Theta ) )  = - Yd/(2*dY);        
          // theta_i,j-1
          A( row, col( i, j - 1, Theta ) )  =   Yd/(2*dY);        
          // RHS      
          B[ row  ]                         = - Guess_laplace[ Psi ] 
                                              + 2 * Xd * Guess_X[ U ] / (zeta0*zeta0) 
                                              + Yd * Guess_Y[ Theta ];
          ++row;

          ////////////////    
          // U equation //
          ////////////////

          // u_i,j+1
          A( row, col( i, j + 1, U ) )      =   Yd*Yd/(dY*dY) + Ydd/(2*dY) + (Guess[Phi] 
                                              + base[oPhi])* Yd / (2*dY);
          // u_i,j
          A( row, col( i, j, U ) )          = - 2*Yd*Yd/(dY*dY) - 2*Xd*Xd
                                              / ( zeta0*zeta0*dX*dX );
          // u_i,j-1
          A( row, col( i, j - 1, U ) )      =   Yd*Yd/(dY*dY) - Ydd/(2*dY) - (Guess[Phi] 
                                              + base[oPhi])* Yd / (2*dY);
          // u_i+1,j
          A( row, col( i + 1, j, U ) )      =   Xd*Xd/(zeta0*zeta0*dX*dX) + Xdd
                                              / ( 2 * zeta0 * zeta0 * dX ) 
                                              + (Guess[Psi] + base[oPsi]*zeta_hat) * Xd 
                                              / ( 2 * dX );
          // u_i-1,j
          A( row, col( i - 1, j, U ) )      =   Xd*Xd/(zeta0*zeta0*dX*dX) - Xdd
                                              / ( 2*zeta0*zeta0*dX ) 
                                              - (Guess[Psi] + base[oPsi]*zeta_hat) * Xd 
                                              / ( 2*dX );
          // phi_i,j
          A( row, col( i, j, Phi ) )        =   ( Yd * Guess_Y[U] + base[oUd] );
          // psi_i,j
          A( row, col( i, j, Psi ) )        =   ( Xd * Guess_X[U] + 0);
          // RHS
          B[ row ]                          = - Guess_laplace[ U ] 
                                              - Yd * (base[oPhi] + Guess[Phi]) 
                                              * Guess_Y[ U ] - Xd 
                                              * (zeta_hat*base[oPsi] + Guess[Psi]) 
                                              * Guess_X[ U ] - base[oUd] * Guess[Phi];
          ++row;

          ////////////////////
          // Theta equation //
          ////////////////////

          // theta_i,j+1
          A( row, col( i, j + 1, Theta ) )  =   Yd*Yd/(dY*dY) + Ydd/(2*dY) 
                                              + (Guess[Phi]+base[oPhi]) * Yd / (2*dY);
          // theta_i,j  
          A( row, col( i, j, Theta ) )      = - 2 * Yd*Yd / (dY*dY) - 2*Xd*Xd 
                                              / (zeta0*zeta0*dX*dX) 
                                              +  2 * (Guess[U]+base[oU]);
          // theta_i,j-1
          A( row, col( i, j - 1, Theta ) )  =   Yd*Yd/(dY*dY) - Ydd/(2*dY) 
                                              - (Guess[Phi]+base[oPhi]) * Yd / (2*dY);
          // theta_i+1,j
          A( row, col( i + 1, j, Theta ) )  =   Xd*Xd/(zeta0*zeta0*dX*dX) + Xdd
                                              / ( 2*zeta0*zeta0*dX ) 
                                              + (Guess[Psi]+zeta_hat*base[oPsi]) * Xd 
                                              / (2*dX);
          // theta_i-1,j
          A( row, col( i - 1, j, Theta ) )  =   Xd*Xd/(zeta0*zeta0*dX*dX) - Xdd
                                              / ( 2*zeta0*zeta0*dX ) 
                                              - (Guess[Psi]+zeta_hat*base[oPsi]) * Xd 
                                              / (2*dX);
          // u_i,j+1
          A( row, col( i, j + 1, U ) )      = - zeta_hat*(Guess[U]+base[oU])*Yd / dY;
          // u_i,j
          A( row, col( i, j, U ) )          = - 2*zeta_hat*(Yd*Guess_Y[U]+base[oUd]) 
                                              + 2*eta*Xd*Guess_X[ U ]/(zeta0*zeta0) 
                                              + 2*(Guess[Theta]+zeta_hat*base[oTheta]) ;
          // u_i,j-1
          A( row, col( i, j - 1, U ) )      =   zeta_hat*(Guess[U]+base[oU])*Yd / dY;
          // u_i+1,j
          A( row, col( i + 1, j, U ) )      =   eta*(Guess[ U ]+base[oU])*Xd
                                              / ( zeta0*zeta0*dX );
          // u_i-1,j
          A( row, col( i - 1, j, U ) )      = - eta*(Guess[ U ]+base[oU])*Xd
                                              / ( zeta0*zeta0*dX );
          // phi_i,j
          A( row, col( i, j, Phi ) )        =   ( Yd*Guess_Y[ Theta ] 
                                              + zeta_hat*base[oThetad] );
          // psi_i,j
          A( row, col( i, j, Psi ) )        =   ( Xd*Guess_X[ Theta ] + base[oTheta] );
          // RHS
          B[ row ]    = - Guess_laplace[ Theta ] - Yd * (base[oPhi] + Guess[Phi]) 
                        * Guess_Y[Theta] - Xd * (zeta_hat * base[oPsi] + Guess[Psi]) 
                        * Guess_X[Theta] + 2 * zeta_hat * Yd * (base[oU] + Guess[U]) 
                        * Guess_Y[U] - 2 * (eta / (zeta0 * zeta0)) * Xd 
                        * (base[oU] + Guess[U]) * Guess_X[U] - 2 * (Guess[Theta] 
                        + zeta_hat * (base[oTheta] - base[oUd])) * Guess[U]
					              - zeta_hat * base[oThetad] * Guess[Phi] - base[oTheta] * Guess[Psi]
					              - 2 * base[oU] * Guess[Theta]; 
          ++row;
        }
        // eta = eta_inf boundary ( top boundary )
        j = Example::M ;
        //
        eta = eta_nodes[ j ];
        Yd = Example::Yd(eta);
        
        // phi_eta = A*(...)
        {
          A( row, col( i, j, Phi ) )        =   3.0 *Yd/ (2*dY);
          A( row, col( i, j - 1, Phi ) )    = - 4.0 *Yd/ (2*dY);
          A( row, col( i, j - 2, Phi ) )    =   1.0 *Yd/ (2*dY);
          A( row, 4 * N_zeta * N_eta )      = - (zeta0 * zeta0 * zeta_hat * zeta_hat 
                                              - eta * eta) / pow(eta*eta 
                                              + zeta0*zeta0*zeta_hat*zeta_hat,2);
          B[ row ]                          = - ( ( 3*Q(i,j,Phi) - 4*Q(i,j-1,Phi) 
                                              + Q(i,j-2,Phi) ) * Yd / (2*dY) 
                                              - Example::A*(zeta0*zeta0*zeta_hat*zeta_hat
                                              - eta*eta) / pow(eta*eta 
                                              + zeta0*zeta0*zeta_hat*zeta_hat,2) );
          ++row;
        }

        // psi_eta = A*(...)
        {
          A( row, col( i, j, Psi ) )        =   3.0 *Yd/ (2*dY);
          A( row, col( i, j - 1, Psi ) )    = - 4.0 *Yd/ (2*dY);
          A( row, col( i, j - 2, Psi ) )    =   1.0 *Yd/ (2*dY);
          A( row, 4 * N_zeta * N_eta )      = - (-2*eta*zeta_hat) 
                                              / pow(eta*eta
                                              + zeta0*zeta0*zeta_hat*zeta_hat,2);
          B[ row ]    = - ( ( 3*Q(i,j,Psi) - 4*Q(i,j-1,Psi) + Q(i,j-2,Psi) ) * Yd / (2*dY) 
                        - Example::A*(-2*eta*zeta_hat) / pow(eta*eta 
                        + zeta0*zeta0*zeta_hat*zeta_hat,2) );
          ++row; 
        }
          
        // u = 0
        A( row, col( i, j, U ) )            =   1;
        B[ row ]                            = - ( Q( i, j, U ) );
        ++row;
        
        // theta = 0
        A( row, col( i, j, Theta ) )        =   1;
        B[ row ]                            = - ( Q(i,j,Theta) );
        ++row;

      } // end of for loop over interior nodes
        
      // zeta_hat = zeta_hat_inf ( right boundary )
      for ( unsigned j = 0; j < Example::M + 1; ++j )
      {
        //offset for global problem
        unsigned i( N_zeta-1 );
        double zeta_hat( zeta_hat_nodes[ i ] );
        double Xd( Example::Xd(zeta_hat) );
        double eta( eta_nodes[j] );
        
        double zeta( zeta0*zeta_hat );
        
        Vector<double> base( farfield.solution().get_interpolated_vars( eta ) );
        
        // zeta_hat*phi_zeta_hat + 2*phi = A*(...)
        A( row, col( i, j, Phi ) )          =   zeta_hat_right*3.*Xd/(2*dX) + 2.0;
        A( row, col( i - 1, j, Phi ) )      = - zeta_hat_right*4.*Xd/(2*dX);
        A( row, col( i - 2, j, Phi ) )      =   zeta_hat_right*Xd/(2*dX);
        A( row, 4 * N_zeta * N_eta )        = - 2*pow(eta,3)/pow(zeta*zeta+eta*eta,2);
        B[ row ]                            = - ( zeta_hat_right*( 3*Q( i, j, Phi) 
                                              - 4*Q( i - 1, j, Phi) 
                                              + Q( i - 2, j, Phi) )*Xd/(2*dX) 
                                              + 2*Q( i, j, Phi) - 2*Example::A*pow(eta,3.)
                                              / pow(zeta*zeta+eta*eta,2) );
        ++row;
        
        // zeta_hat*psi_zeta + psi = A*(...)   
        A( row, col( i, j, Psi ) )          =   zeta_hat_right*3.*Xd/(2*dX) + 1.0;
        A( row, col( i - 1, j, Psi ) )      = - zeta_hat_right*4.*Xd/(2*dX);
        A( row, col( i - 2, j, Psi ) )      =   zeta_hat_right*Xd/(2*dX);
        A( row, 4 * N_zeta * N_eta )        = - 2.0*zeta_hat*pow(eta,2.0)
                                              / pow(zeta*zeta+eta*eta,2);
        B[ row ]                            = - ( zeta_hat_right*( 3*Q( i, j, Psi) 
                                              - 4*Q( i - 1, j, Psi) 
                                              + Q( i - 2, j, Psi) )*Xd/(2*dX) 
                                              + Q( i, j, Psi)
                                              - 2.0*Example::A*zeta_hat*pow(eta,2.0)
                                              / pow(zeta*zeta+eta*eta,2) );
        ++row;

        // zeta_hat*u_zeta_hat + 2*u = 0
        A( row, col( i, j, U ) )            =   zeta_hat_right*3.*Xd/(2*dX) + 2.0;
        A( row, col( i - 1, j, U ) )        = - zeta_hat_right*4.*Xd/(2*dX);
        A( row, col( i - 2, j, U ) )        =   zeta_hat_right*Xd/(2*dX);
        B[ row  ]                           = - ( zeta_hat_right*( 3*Q( i, j, U) 
                                              - 4*Q( i - 1, j, U) + Q( i - 2, j, U) )*Xd
                                              / (2*dX) + 2.0 * Q( i, j, U) );
        ++row;
        
        // zeta_hat*theta_zeta_hat + theta = 0
        A( row, col( i, j, Theta )  )       =   zeta_hat_right*3.*Xd/(2*dX) + 1.0;
        A( row, col( i - 1, j, Theta ) )    = - zeta_hat_right*4.*Xd/(2*dX);
        A( row, col( i - 2, j, Theta ) )    =   zeta_hat_right*Xd/(2*dX);
        B[ row ]                            = - ( zeta_hat_right*( 3*Q( i, j, Theta) 
                                              - 4*Q( i - 1, j, Theta) 
                                              + Q( i - 2, j, Theta) )*Xd/(2*dX) 
                                              + Q(i, j, Theta) );
        ++row;

      }
        
      // WE ALSO NEED A FINAL CLOSING CONDITION FOR THE DISPLACEMENT COEFFICIENT "A"
      
      // zeta0^2 * zeta_hat_max theta( zeta=zeta_max, eta=0 ) = A*otheta(0)
      Vector<double> base( farfield.solution().get_interpolated_vars( 0.0 ) );
      A( 4 * N_eta * N_zeta, 4 * N_eta * N_zeta ) = -base[otheta];
      A( 4 * N_eta * N_zeta, 4 * N_eta * (N_zeta - 1) + 3 ) =  zeta0*zeta0*zeta_hat_right;
      // RHS
      B[ row ] = -( Q( N_zeta-1, 0, Theta )*zeta0*zeta0*zeta_hat_right 
                 - Example::A*base[otheta] );   
    
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
        
      for ( unsigned i = 0; i <= Example::N ; ++i )
      {
        for ( unsigned j = 0; j <= Example::M ; ++j )
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
    for ( unsigned i = 0; i < N_zeta; ++i )
    {
      double zeta_hat=zeta_hat_nodes[i];
      for ( unsigned j = 0; j < N_eta; ++j )
      {
        double eta=eta_nodes[j];
        // first 4 values output are the without the underlying 2D base flow
        Q_output( i, j, 0 ) = Q( i, j, Phi);
        Q_output( i, j, 1 ) = Q( i, j, Psi);
        Q_output( i, j, 2 ) = Q( i, j, U);
        Q_output( i, j, 3 ) = Q( i, j, Theta);
        // second 4 values are the "full" solution, but still with the zeta0 scaling
        Q_output( i, j, 4 ) =   Q( i, j, Phi)   
                              + farfield.solution().get_interpolated_vars( eta )[oPhi];
        Q_output( i, j, 5 ) = Q( i, j, Psi)   
                   + zeta_hat * farfield.solution().get_interpolated_vars( eta )[oPsi];
        Q_output( i, j, 6 ) =   Q( i, j, U)
                              + farfield.solution().get_interpolated_vars( eta )[oU];
        Q_output( i, j, 7 ) = Q( i, j, Theta) 
                   + zeta_hat * farfield.solution().get_interpolated_vars( eta )[oTheta];
      }
    }
    // Convert to string //TODO need a utility function to do this
    std::stringstream ss;
    ss << zeta0;
    std::string zeta0_str = ss.str(); 

    Q_output.dump_gnu( Example::output_path + "Qout_" + zeta0_str + ".dat" );

    Vector<double> base( farfield.solution().get_interpolated_vars( 0.0 ) ); 
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
