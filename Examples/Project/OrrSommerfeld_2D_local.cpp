// Solve the 2D Orr-Sommerfeld equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"

//enum{ v, w, q, s }; enum already included in OrrSommerfeld_2D.h

namespace TSL
{
  std::size_t eta_intervals;

  std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
  {
    // Return the column number for the kth variable at node (i,j)
    return 4 * ( i * ( eta_intervals + 1 ) + j ) + k;
  }

}

using namespace std;
using namespace TSL;

class mySelfSimInjection : public SelfSimInjection {
public:
  // Define the injection function
  double Phi_w_func( const double& hzeta ){
    return - K * exp( - hzeta * hzeta );
  }
}; // End of class mySelfSimInjection


int main()
{
  cout << "*** ------- Solving the 2D OrrSommerfeld equation (Global) ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 32.0 ); // Size of the domain in the zeta_hat direction
  double eta_top( 32.0 );     // Size of the domain in the eta direction
  std::size_t N( 180 );       // Number of intervals in the zeta_hat direction
  std::size_t M( 180 );       // Number of intervals in the eta direction
  std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.5 );         // Hartree parameter
  double zeta0( 1.0 );        // Transpiration width
  double K( 9.0 );           // Transpiration parameter ( +ve = blowing )
  double alpha( 0.4 );        // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );   // Local Reynolds number

  std::complex<double> target( 0.76, 0.0 ); // Target for eigensolver

  // Solve the self similar injection flow
  mySelfSimInjection SSI;
  SSI.hzeta_right() = hzeta_right;
  SSI.eta_top() = eta_top;
  SSI.hzeta_intervals() = N;
  SSI.eta_intervals() = M;
  SSI.base_intervals() = MB;
  SSI.hartree() = beta;
  SSI.injection_width() = zeta0;
  SSI.injection() = K;
  SSI.wavenumber() = alpha;
  SSI.local_Reynolds() = Rx;
  SSI.set_mesh( "NONUNIFORM" );
  SSI.set_base_flow( "2D" );
  SSI.speed_up( false );


  cout << "*** Solving the self-similar flow on a coarse mesh ***" << endl;
  cout << "  * N = " << N << ", M = " << M << endl;
  Timer timer;
  timer.start();

  // Solve self-similar equations on a coarse grid
  SSI.set_output_path();
  SSI.mesh_setup();
  SSI.solve();
  /*Vector<double> ETA_NODES      = SSI.eta_nodes();
  Vector<double> HZETA_NODES    = SSI.hzeta_nodes();
  Vector<double> X_NODES        = SSI.x_nodes();
  Vector<double> Y_NODES        = SSI.y_nodes();
  Vector<double> BASE_ETA_NODES = SSI.base_eta_nodes();

  TwoD_node_mesh<double> sol( HZETA_NODES, ETA_NODES, 8 ); // Mesh for storing the solution
  OneD_node_mesh<double> base( BASE_ETA_NODES, 6 );

  // Don't bother solving it all again if the solution file already exists
  bool exists;
  exists = Utility::file_exists( SSI.output_path() + "Qout_" + Utility::stringify( zeta0 ) + ".dat" );
  bool base_exists;
  base_exists = Utility::file_exists( SSI.output_path() + "Base_soln.dat" );

  try
  {
    if ( !exists || !base_exists ){
      SSI.solve();
      SSI.output();
      SSI.output_base_solution();
      sol = SSI.solution();
      base = SSI.base_flow_solution();
      cout << "  * zeta0 = " << SSI.injection_width() << ", A = " << SSI.mass_flux() << endl;
     }
    if ( exists ){
      cout << "--- Reading solution from file" << endl;
      sol.read( SSI.output_path() + "Qout_" + Utility::stringify( zeta0 ) + ".dat" );
      cout << "--- Finished reading" << endl;
    }

    if ( base_exists ){
      base.read( SSI.output_path() + "Base_soln.dat" );
      std::cout << "  * UB'(eta=0) =" << base( 0, 1 ) << std::endl;
    }

    if ( exists && base_exists ){
      SSI.set_solved( true );
      SSI.set_solution( sol );
      SSI.set_base_solution( base );
    }
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED (SSI) \033[0m\n";
    assert( false );
  }*/

  // Setup the generalised eigenvalue problem A p = c B p (solved using SLEPc)
  cout << "*** Setting up the generalised eigenvalue problem ***" << endl;
  cout << "--- K = " << K << ", alpha = " << alpha << ", Rx^1/2 = " << sqrt(Rx) << endl;

  // Create the OrrSommerfeld_2D object
  std::size_t nev( 1 );
  OrrSommerfeld_2D orrsommerfeld_2D( SSI, alpha, Rx, nev );
  orrsommerfeld_2D.set_region(0.1,1.0,-1.0,1.0);
  orrsommerfeld_2D.set_target( target );
  orrsommerfeld_2D.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D.calc_eigenvectors() = true;
  orrsommerfeld_2D.solve_evp();
  orrsommerfeld_2D.output();

  // Return the eigenvalues and eigenvectors on the coarse mesh
  Vector< std::complex<double> > evals;
  evals = orrsommerfeld_2D.eigenvalues();
  TwoD_node_mesh< std::complex<double> > evecs;
  evecs = orrsommerfeld_2D.eigenvectors(); // v, w, q, s

  // Redefine the SSI mesh and resolve
  cout << "*** Solving the self-similar flow on a refined mesh ***" << endl;
  N = 200;
  M = 200;
  MB = M * 100;
  cout << "  * N = " << N << ", M = " << M << endl;
  SSI.hzeta_intervals() = N;
  SSI.eta_intervals() = M;
  SSI.base_intervals() = MB;
  SSI.mesh_setup();
  SSI.solve();
  TSL::eta_intervals = M;

  // Remesh the eigenvectors onto the refined mesh
  Vector<double> ETA_NODES      = SSI.eta_nodes();
  Vector<double> HZETA_NODES    = SSI.hzeta_nodes();
  Vector<double> X_NODES        = SSI.x_nodes();
  Vector<double> Y_NODES        = SSI.y_nodes();
  Vector<double> BASE_ETA_NODES = SSI.base_eta_nodes();

  TwoD_node_mesh<double> sol( HZETA_NODES, ETA_NODES, 8 ); // Mesh for storing the solution
  OneD_node_mesh<double> base( BASE_ETA_NODES, 6 );
  sol = SSI.solution();
  base = SSI.base_flow_solution();

  evecs.remesh1( HZETA_NODES, ETA_NODES );

  // Normalise the eigenvectors

  int v = static_cast<int>(OS_2D::v);
  int w = static_cast<int>(OS_2D::w);
  int q = static_cast<int>(OS_2D::q);
  int s = static_cast<int>(OS_2D::s);
  double norm;
  norm = real( evecs.square_integral2D( v ) + evecs.square_integral2D( w ) );
            //+ evecs.square_integral2D( q ) + evecs.square_integral2D( s ) );
  for ( std::size_t i=0; i<evecs.xnodes().size(); ++i )
  {
    for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
    {
      evecs( i, j, v ) = evecs( i, j, v ) / norm;
      evecs( i, j, w ) = evecs( i, j, w ) / norm;
      evecs( i, j, q ) = evecs( i, j, q ) / norm;
      evecs( i, j, s ) = evecs( i, j, s ) / norm;
    }
  }

  // Solve the nonlinear system
  cout << "*** ------- Solving the 2D OrrSommerfeld equation (Local) ------- ***" << endl;

  // Current guess g
  TwoD_node_mesh< std::complex<double> > Q( X_NODES, Y_NODES, 4 );
  std::complex<double> c_g( evals[0] );

  // Set the initial guess
  for ( std::size_t i = 0; i < N + 1; ++i )
  {
    for ( std::size_t j = 0; j < M + 1; ++j )
    {
      Q( i, j, v )  = evecs( i, j, v );
      Q( i, j, w )  = evecs( i, j, w );
      Q( i, j, q )  = evecs( i, j, q );
      Q( i, j, s )  = evecs( i, j, s );
    }
  }

  // Normalise based on the extra condition q(0,0) = 1
  std::complex<double> lambda( Q( 0, 0, q ) );
  for ( std::size_t i = 0; i < N + 1; ++i )
  {
    for ( std::size_t j = 0; j < M + 1; ++j )
    {
      Q( i, j, v )  = Q( i, j, v ) / lambda;
      Q( i, j, w )  = Q( i, j, w ) / lambda;
      Q( i, j, q )  = Q( i, j, q ) / lambda;
      Q( i, j, s )  = Q( i, j, s ) / lambda;
    }
  }

  //cout << "Q( 0, 0, q ) = " << Q( 0, 0, q ) << endl;

  // Step sizes
  const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
  const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );
  // Vector for the RHS of the matrix problem
  std::size_t size( 4 * ( M + 1 ) * ( N + 1 ) + 1 );
  Vector< std::complex<double> > B( size, 0.0 );

  /* Iterate to a solution */
  double max_residual( 0.0 );             // Maximum residual
  std::size_t iteration( 0 );             // Initialise iteration counter
  std::size_t max_iterations( 20 );       // Maximum number of iterations

  std::complex<double> iaR ( 0.0, 1.0 / ( alpha * sqrt( Rx ) ) );
  std::complex<double> Ralpha ( 1.0 / sqrt( Rx ), alpha );

  do {
      SparseMatrix< std::complex<double> > A( size, size );
      std::cout << "  * Assembling sparse matrix problem" << std::endl;

      Timer timer;                                        // Timer
      timer.start();

      std::size_t row( 0 );                               // Initialise row counter

      // hzeta = 0 boundary ( left boundary )
      std::size_t i( 0 );

      for ( std::size_t j = 0; j < M + 1 ; ++j )
      {
        double hzeta( HZETA_NODES[ 0 ] );
        double Xd( SSI.mesh_Xd( hzeta ) );
        double Xdd( SSI.mesh_Xdd( hzeta ) );
        // eta location
        double eta( ETA_NODES[ j ] );
        double Yd( SSI.mesh_Yd( eta ) );
        double Ydd( SSI.mesh_Ydd( eta ) );

        // Symmetry boundary conditions

        // w = 0
        A( row, col( i, j, w ) ) = 1.;

        B[ row ] = - Q( i, j, w );
        ++row;

        // v_hzeta = 0
        A( row, col( i, j, v ) )       = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, v ) )   =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, v ) )   = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, v ) + 4. * Q( i + 1, j, v )
                                - Q( i + 2, j, v ) ) / ( 2 * dX ) );
        ++row;

        // s = 0
        A( row, col( 0, j, s ) ) = 1.;

        B[ row ] = - Q( i, j, s );
        ++row;

        // q_hzeta = 0
        A( row, col( i, j, q ) )      = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, q ) )  =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, q ) )  = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, q ) + 4. * Q( i + 1, j, q )
                                - Q( i + 2, j, q ) ) / ( 2 * dX ) );
        ++row;

      } // End of for loop over LHS eta nodes

      // Interior points between the hzeta boundaries
      for ( std::size_t i = 1; i < N; ++i )
      {
        // hzeta location
        double hzeta( HZETA_NODES[ i ] );
        double Xd( SSI.mesh_Xd( hzeta ) );
        double Xdd( SSI.mesh_Xdd( hzeta ) );

        // eta = 0 boundary ( bottom boundary )
        std::size_t j( 0 );
        double eta( ETA_NODES[ j ] );
        double Yd( SSI.mesh_Yd( eta ) );
        double Ydd( SSI.mesh_Ydd( eta ) );

        // v = 0
        A( row, col( i, j, v ) ) = 1.;

        B[ row ] = - Q( i, j, v );
        ++row;

        // w = 0
        A( row, col( i, j, w ) ) = 1.;

        B[ row ] = - Q( i, j, w );
        ++row;

        //s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
        A( row, col( i, j, s ) )     =   1.;
        A( row, col( i, j + 1, w ) ) = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 2, w ) ) =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 3, w ) ) =   iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s )
                   + iaR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) ) * Q( i, j + 1, w)
                   - iaR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) ) * Q( i, j + 2, w )
                   - iaR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, w );
        ++row;

        //q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
        A( row, col( i, j, q ) )     =   1.;
        A( row, col( i, j + 1, v ) ) = - 6. * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 2, v ) ) =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 3, v ) ) = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q )
                   + iaR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j + 1, v )
                   - iaR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 2, v )
                   + iaR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, v );
        ++row;


        // Main interior grid points
        for ( std::size_t j = 1; j < M; ++j )
        {
          // eta location
          double eta( ETA_NODES[ j ] );
          double Yd( SSI.mesh_Yd( eta ) );
          double Ydd( SSI.mesh_Ydd( eta ) );

          // Laplacian coefficients for finite-differencing

          // X(i,j-1)
          double laplace_1 =  ( Yd*Yd / (dY*dY) - Ydd / (2.*dY) ) ;
          // X(i-1,j)
          double laplace_3 = ( Xd*Xd / (dX*dX) - Xdd / (2.*dX) ) / ( zeta0 * zeta0 );
          // X(i,j)
          double laplace_4 = - 2 *( Yd*Yd / (dY*dY) + Xd*Xd / ( zeta0 * zeta0 * dX * dX ) )
                             - alpha * alpha;
          // X(i+1,j)
          double laplace_5 = ( Xdd / (2.*dX) + Xd*Xd / (dX*dX) ) / ( zeta0 * zeta0 );
          // X(i,j+1)
          double laplace_7 = ( Yd*Yd / (dY*dY) + Ydd / (2.*dY) );

          // Self similar solution
          double U = sol( i, j, 6 ); // UB + U_pert
          double U_hzeta = Xd * ( sol( i + 1, j, 2 ) - sol( i - 1, j, 2 ) ) / ( 2 * dX ); // U_pert_hzeta
          double U_eta = Yd * ( sol( i, j + 1, 2 ) - sol( i, j - 1, 2 ) ) / ( 2 * dY ); // U_pert_eta
          double UBd = base.get_interpolated_vars( eta )[ 1 ];
          U_eta += UBd;

          double UB = base.get_interpolated_vars( eta )[ 0 ];
          double PhiB = base.get_interpolated_vars( eta )[ 2 ];
          double U_eta_eta = laplace_1 * sol( i, j - 1, 2 ) + laplace_7 * sol( i, j + 1, 2 );
          U_eta_eta += - 2 * Yd * Yd * sol( i, j, 2 ) / ( dY * dY );
          U_eta_eta += beta * ( UB * UB - 1 ) - PhiB * UBd; // + UBdd = beta * [UB^2 - 1 ]- PhiB * UBd

          double U_hzeta_hzeta = laplace_3 * sol( i - 1, j, 2 ) + laplace_5 * sol( i + 1, j, 2 );
          U_hzeta_hzeta += - 2 * Xd * Xd * sol( i, j, 2 ) / ( dX * dX );

          double U_eta_hzeta = Yd * Xd * ( sol( i + 1, j + 1, 2 ) + sol( i - 1, j - 1, 2 )
                                         - sol( i + 1, j - 1, 2 ) - sol( i - 1, j + 1, 2 ) )
                                         / ( 4 * dY * dX );

          ///////////////////////////////
          // 2D OrrSommerfeld equation //
          ///////////////////////////////

          ////////////////
          // v equation //
          ////////////////

          // ( i / alpha * Rx^(1/2) ) * Laplacian of v
          A( row, col( i, j - 1, v ) )     = iaR * laplace_1;
          A( row, col( i - 1, j, v ) )     = iaR * laplace_3;
          A( row, col( i, j, v ) )         = iaR * laplace_4;
          A( row, col( i + 1, j, v ) )     = iaR * laplace_5;
          A( row, col( i, j + 1, v ) )     = iaR * laplace_7;

          // + ( U - c_g ) * v
          A( row, col( i, j, v ) )        +=  U - c_g;

          // + ((1 - beta) / Rx) * v
          A( row, col( i, j, v ) )        +=  ( 1. - beta ) / Rx;

          // - v_g * c
          A( row, size - 1 ) = - Q( i, j, v );

          // - q
          A( row, col( i, j, q ) )         = - 1.;

          B[ row ] =  - iaR * ( Q( i, j - 1, v ) * laplace_1
                             +  Q( i - 1, j, v ) * laplace_3
                             +  Q( i, j, v )     * laplace_4
                             +  Q( i + 1, j, v ) * laplace_5
                             +  Q( i, j + 1, v ) * laplace_7 )
                      - ( U - c_g + ( ( 1. - beta ) / Rx ) ) * Q( i, j, v )
                      + Q( i, j, q );
          ++row;

          ////////////////
          // w equation //
          ////////////////

          // ( i / alpha * Rx^(1/2) ) * Laplacian of v
          A( row, col( i, j - 1, w ) )     = iaR * laplace_1;
          A( row, col( i - 1, j, w ) )     = iaR * laplace_3;
          A( row, col( i, j, w ) )         = iaR * laplace_4;
          A( row, col( i + 1, j, w ) )     = iaR * laplace_5;
          A( row, col( i, j + 1, w ) )     = iaR * laplace_7;

          // + ( U - c_g ) * w
          A( row, col( i, j, w ) )        +=  U - c_g;

          // + ((1 - beta) / Rx) * w
          A( row, col( i, j, w ) )        +=  ( 1. - beta ) / Rx;

          // - w_g * c
          A( row, size - 1 ) = - Q( i, j, w );

          // - s
          A( row, col( i, j, s ) )         = - 1.;

          B[ row ] =  - iaR * ( Q( i, j - 1, w ) * laplace_1
                             +  Q( i - 1, j, w ) * laplace_3
                             +  Q( i, j, w )     * laplace_4
                             +  Q( i + 1, j, w ) * laplace_5
                             +  Q( i, j + 1, w ) * laplace_7 )
                      - ( U - c_g + ( ( 1. - beta ) / Rx ) ) * Q( i, j, w )
                      + Q( i, j, s );
          ++row;

          ////////////////
          // q equation //
          ////////////////

          // q_{eta eta}
          A( row, col( i, j - 1, q ) )      =   Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ;
          A( row, col( i, j, q ) )          = - 2 * Yd * Yd / ( dY * dY );
          A( row, col( i, j + 1, q ) )      =   Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ;

          // - alpha^2 * q
          A( row, col( i, j, q ) )         += - alpha * alpha;

          // + ( 1 / zeta0 ) * s_{hzeta eta}
          A( row, col( i + 1, j + 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j - 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i + 1, j - 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j + 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );

          // + ( 1 / zeta0 ) * U_{eta} * w_{hzeta}
          A( row, col( i + 1, j, w ) )     =   U_eta * Xd / ( 2 * dX * zeta0 );
          A( row, col( i - 1, j, w ) )     = - U_eta * Xd / ( 2 * dX * zeta0 );

          // - ( 1 / zeta0 ) * U_{hzeta} * w_{eta}
          A( row, col( i, j + 1, w ) )     = - U_hzeta * Yd / ( 2 * dY * zeta0 );
          A( row, col( i, j - 1, w ) )     =   U_hzeta * Yd / ( 2 * dY * zeta0 );

          // - ( 1 / zeta0 ) * U_{hzeta eta} * w
          A( row, col( i, j, w ) )         = - U_eta_hzeta / zeta0;

          // - U_{eta eta} * v
          A( row, col( i, j, v ) )          = - U_eta_eta;

          B[ row ] = - ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) * Q( i, j - 1, q )
                     - ( - 2 * Yd * Yd / ( dY * dY ) ) * Q( i, j, q )
                     - ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) * Q( i, j + 1, q )
                     + alpha * alpha * Q( i, j, q )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j + 1, s )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j - 1, s )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j - 1, s )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j + 1, s )
                     - ( U_eta * Xd / ( 2 * dX * zeta0 ) ) * Q( i + 1, j, w )
                     + ( U_eta * Xd / ( 2 * dX * zeta0 ) ) * Q( i - 1, j, w )
                     + ( U_hzeta * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j + 1, w )
                     - ( U_hzeta * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j - 1, w )
                     + ( U_eta_hzeta / zeta0 ) * Q( i, j, w )
                     + ( U_eta_eta ) * Q( i, j, v );
          ++row;



          ////////////////
          // s equation //
          ////////////////

          // ( 1 / zeta0^2 ) * s_{hzeta hzeta}
          A( row, col( i - 1, j, s ) )      =   ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) )
                                              / ( zeta0 * zeta0 );
          A( row, col( i, j, s ) )          = - 2 * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
          A( row, col( i + 1, j, s ) )      =   ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) )
                                              / ( zeta0 * zeta0 );

          // - alpha^2 * s
          A( row, col( i, j, s ) )         += - alpha * alpha;

          // + ( 1 / zeta0 ) * q_{hzeta eta}
          A( row, col( i + 1, j + 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j - 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i + 1, j - 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j + 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );

          // + ( 1 / zeta0 ) * U_{hzeta} * v_{eta}
          A( row, col( i, j + 1, v ) )     =   U_hzeta * Yd / ( 2 * dY * zeta0 );
          A( row, col( i, j - 1, v ) )     = - U_hzeta * Yd / ( 2 * dY * zeta0 );

          // - ( 1 / zeta0 ) * U_{eta} * v_{hzeta}
          A( row, col( i + 1, j, v ) )     = - U_eta * Xd / ( 2 * dX * zeta0 );
          A( row, col( i - 1, j, v ) )     =   U_eta * Xd / ( 2 * dX * zeta0 );

          // - ( 1 / zeta0 ) * U_{eta hzeta} * v
          A( row, col( i, j, v ) )         = - U_eta_hzeta / zeta0;

          // - ( 1 / zeta0^2 ) * U_{hzeta hzeta} * w
          A( row, col( i, j, w ) )         = - U_hzeta_hzeta / ( zeta0 * zeta0 );

          B[ row ] = - ( ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) ) / ( zeta0 * zeta0 ) ) * Q( i - 1, j, s )
                     - ( - 2 * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i, j, s )
                     - ( ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) ) / ( zeta0 * zeta0 ) ) * Q( i + 1, j, s )
                     + alpha * alpha * Q( i, j, s )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j + 1, q )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j - 1, q )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j - 1, q )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j + 1, q )
                     - ( U_hzeta * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j + 1, v )
                     + ( U_hzeta * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j - 1, v )
                     + ( U_eta * Xd / ( 2 * dX * zeta0 ) ) * Q( i + 1, j, v )
                     - ( U_eta * Xd / ( 2 * dX * zeta0 ) ) * Q( i - 1, j, v )
                     + ( U_eta_hzeta / zeta0 ) * Q( i, j, v )
                     + ( U_hzeta_hzeta / ( zeta0 * zeta0 ) ) * Q( i, j, w );
          ++row;

      }

      // eta = eta_inf boundary ( top boundary )
      j = M ;
      eta = ETA_NODES[ j ];
      Yd = SSI.mesh_Yd( eta );

      // v = 0
      A( row, col( i, j, v ) ) = 1.;

      B[ row ] = - Q( i, j, v );
      ++row;

      // w = 0
      A( row, col( i, j, w ) ) = 1.;

      B[ row ] = - Q( i, j, w );
      ++row;

      // s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
      A( row, col( i, j, s ) )        =   1.;
      A( row, col( i, j - 1, w ) )    = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) );
      A( row, col( i, j - 2, w ) )    =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) );
      A( row, col( i, j - 3, w ) )    =   iaR * Yd * Yd / ( dY * dY );

      B[ row ] = - Q( i, j, s )
                 + iaR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) ) * Q( i, j - 1, w)
                 - iaR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) ) * Q( i, j - 2, w )
                 - iaR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, w );
      ++row;

      // q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
      A( row, col( i, j, q ) )        =   1.;
      A( row, col( i, j - 1, v ) )    = - 6. * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, j - 2, v ) )    =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, j - 3, v ) )    = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

      B[ row ] = - Q( i, j, q )
                 + iaR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j - 1, v )
                 - iaR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 2, v )
                 + iaR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, v );
      ++row;

    } // End of loop over interior nodes

    // hzeta = hzeta_inf boundary ( right boundary )
    for ( std::size_t j = 0; j < M + 1; ++j )
    {
      std::size_t i( N );
      double hzeta( HZETA_NODES[ i ] );
      double Xd( SSI.mesh_Xd( hzeta ) );
      double Xdd( SSI.mesh_Xdd( hzeta ) );

      // eta location
      double eta( ETA_NODES[ j ] );
      double Yd( SSI.mesh_Yd( eta ) );
      double Ydd( SSI.mesh_Ydd( eta ) );

      // w = 0
      A( row, col( i, j, w ) ) = 1.;

      B[ row ] = - Q( i, j, w );
      ++row;

      // v = 0
      A( row, col( i, j, v ) ) = 1.;

      B[ row ] = - Q( i, j, v );
      ++row;

      // q - (1/zeta0^2) * (i / (alpha*Rx^(1/2))) * v_{hzeta hzeta} = 0
      A( row, col( i, j, q ) )        =   1.;
      A( row, col( i - 1, j, v ) )    = - iaR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) )
                                    / ( zeta0 * zeta0 );
      A( row, col( i - 2, j, v ) )    =   iaR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) )
                                    / ( zeta0 * zeta0 );
      A( row, col( i - 3, j, v ) )    =   iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );

      B[ row ] = - Q( i, j, q )
                 + ( iaR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) ) / ( zeta0 * zeta0 ) ) * Q( i - 1, j, v )
                 - ( iaR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) ) / ( zeta0 * zeta0 ) ) * Q( i - 2, j, v )
                 - ( iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i - 3, j, v );
      ++row;

      // s - (1/zeta0^2) * (i / (alpha*Rx^(1/2))) * w_{hzeta hzeta} = 0
      A( row, col( i, j, s ) )        =   1.;
      A( row, col( i - 1, j, w ) )    = - 6. * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
      A( row, col( i - 2, j, w ) )    =   ( 3. / 2. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
      A( row, col( i - 3, j, w ) )    = - ( 2. / 9. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );

      B[ row ] = - Q( i, j, s )
                 + ( 6. * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i - 1, j, w )
                 - ( ( 3. / 2. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i - 2, j, w )
                 + ( ( 2. / 9. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i - 3, j, w );
      ++row;

    } // End of loop over RHS eta nodes


    // TODO extra condition for c ??? q(0,0) = 1 ???
    A( row, col( 0, 0, q ) ) = 1.;
    B[ row ] = 1. - Q( 0, 0, q );
    ++row;

    // v_eta(0,0) = 1 ???
    /*double eta( ETA_NODES[ 0 ] );
    double Yd( SSI.mesh_Yd( eta ) );
    double Ydd( SSI.mesh_Ydd( eta ) );
    A( row, col( 0, 0, v ) ) = - 3 * Yd / ( 2 * dY );
    A( row, col( 0, 1, v ) ) =   4 * Yd / ( 2 * dY );
    A( row, col( 0, 2, v ) ) = - 1 * Yd / ( 2 * dY );
    B[ row ] = 1.0 + ( 3 * Yd / ( 2 * dY ) ) * Q( 0, 0, v )
                   - ( 4 * Yd / ( 2 * dY ) ) * Q( 0, 1, v )
                   + ( 1 * Yd / ( 2 * dY ) ) * Q( 0, 2, v );
    ++row;*/

    max_residual = B.norm_inf();
    std::cout << "***                                              Maximum residual = "
              << max_residual << std::endl;

    // Solve the sparse system
    Vector< std::complex<double> > x;
    x = A.solve( B );
    B = x;

    timer.print();
    timer.stop();

    // Update the known values using the correction which we just found
    for ( std::size_t i = 0; i < N + 1; ++i )
    {
      for ( std::size_t j = 0; j < M + 1; ++j )
      {
        Q( i, j, v ) += B[ col( i, j, v ) ];
        Q( i, j, w ) += B[ col( i, j, w ) ];
        Q( i, j, q ) += B[ col( i, j, q ) ];
        Q( i, j, s ) += B[ col( i, j, s ) ];
      }
    }

    c_g += B[ size - 1 ];

    std::cout << "***    Iteration = " << iteration
              << "    Maximum correction = " << B.norm_inf() << std::endl;

    ++iteration;
  }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) );


  cout << "  * c_refined = " << c_g << endl;


  timer.print();
  timer.stop();

  Q.dump( "./DATA/OS2D_local.dat" );

	cout << "FINISHED" << endl;

}
