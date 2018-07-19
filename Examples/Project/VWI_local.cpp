// Solve the 2D Orr-Sommerfeld equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"

enum{ v, w, q, s, Phi, Psi, U, Theta };

namespace TSL
{
  std::size_t eta_intervals;

  std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
  {
    // Return the column number for the kth variable at node (i,j)
    return 8 * ( i * ( eta_intervals + 1 ) + j ) + k;
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
  double hzeta_right( 32.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 32.0 );           // Size of the domain in the eta direction
  std::size_t N( 150 );       // Number of intervals in the zeta_hat direction
  std::size_t M( 150 );       // Number of intervals in the eta direction
  std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 7.5 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.1 );              // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );         // Local Reynolds number
  double Sigma( 0.0 );              // Wave amplitude

  std::complex<double> target(0.8,-0.01); // Target for eigensolver

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
  //SSI.wave_amplitude() = Sigma;
  //SSI.local_Reynolds() = Rx;
  //SSI.wavenumber() = alpha;


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
  TwoD_node_mesh< std::complex<double> > Q( X_NODES, Y_NODES, 8 );
  std::complex<double> c_g( evals[0] );

  // Set the initial guess
  for ( std::size_t i = 0; i < N + 1; ++i )
  {
    for ( std::size_t j = 0; j < M + 1; ++j )
    {
      Q( i, j, v )     = evecs( i, j, v );
      Q( i, j, w )     = evecs( i, j, w );
      Q( i, j, q )     = evecs( i, j, q );
      Q( i, j, s )     = evecs( i, j, s );
      Q( i, j, Phi )   = sol( i, j, 0 );
      Q( i, j, Psi )   = sol( i, j, 1 );
      Q( i, j, U )     = sol( i, j, 2);
      Q( i, j, Theta ) = sol( i, j, 3 );
    }
  }
  // Step sizes
  const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
  const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );
  // Vector for the RHS of the matrix problem
  std::size_t size( 8 * ( M + 1 ) * ( N + 1 ) + 1 );
  Vector< std::complex<double> > B( size, 0.0 );
  //cout << " B.size() = " << B.size() << endl;
  /* Iterate to a solution */
  double max_residual( 0.0 );             // Maximum residual
  std::size_t iteration( 0 );             // Initialise iteration counter
  std::size_t max_iterations( 20 );       // Maximum number of iterations

  std::complex<double> iaR ( 0.0, 1.0 / ( alpha * sqrt( Rx ) ) );
  std::complex<double> Ralpha ( 1.0 / sqrt( Rx ), alpha );

  do {
      SparseMatrix< std::complex<double> > A( size, size );
      //cout << "A.rows() = " << A.rows() << endl;
      //cout << "A.cols() = " << A.cols() << endl;
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
        A( row, col( i, j, w ) ) = 1;

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
        A( row, col( 0, j, s ) ) = 1;

        B[ row ] = - Q( i, j, s );
        ++row;

        // q_hzeta = 0
        A( row, col( i, j, q ) )      = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, q ) )  =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, q ) )  = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, q ) + 4. * Q( i + 1, j, q )
                                - Q( i + 2, j, q ) ) / ( 2 * dX ) );
        ++row;

        // Phi_hzeta = 0
        A( row, col( i, j, Phi ) )      = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, Phi ) )  =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, Phi ) )  = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, Phi ) + 4. * Q( i + 1, j, Phi )
                             - Q( i + 2, j, Phi ) ) / ( 2 * dX ) );
        ++row;

        // Psi = 0
        A( row, col( i, j, Psi ) )      =   1;

        B[ row ]                        = - Q( i, j, Psi );
        ++row;

        // U_hzeta = 0
        A( row, col( i, j, U ) )        = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, U ) )    =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, U ) )    = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, U ) + 4. * Q( i + 1, j, U )
                             - Q( i + 2, j, U ) ) / ( 2 * dX ) );
        ++row;

        // Theta = 0
        A( row, col( i, j, Theta ) )    =   1;

        B[ row ]                        = - Q( i, j, Theta );
        ++row;

      } // End of for loop over LHS eta nodes

      // Interior points between the hzeta boundaries
      for ( std::size_t i = 1; i < N; ++i )
      {
        // hzeta location
        double hzeta( HZETA_NODES[ i ] );
        double Xd( SSI.mesh_Xd( hzeta ) );
        double Xdd( SSI.mesh_Xdd( hzeta ) );

        // Wall transpiration
        double Phi_w( SSI.Phi_w_func( hzeta ) );

        // eta = 0 boundary ( bottom boundary )
        std::size_t j( 0 );
        double eta( ETA_NODES[ j ] );
        double Yd( SSI.mesh_Yd( eta ) );
        double Ydd( SSI.mesh_Ydd( eta ) );

        // v = 0
        A( row, col( i, j, v ) ) = 1;

        B[ row ] = - Q( i, j, v );
        ++row;

        // w = 0
        A( row, col( i, j, w ) ) = 1;

        B[ row ] = - Q( i, j, w );
        ++row;

        //s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
        A( row, col( i, j, s ) )     =   1;
        A( row, col( i, j + 1, w ) ) = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 2, w ) ) =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 3, w ) ) =   iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s )
                   + iaR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) ) * Q( i, j + 1, w)
                   - iaR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) ) * Q( i, j + 2, w )
                   - iaR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, w );
        ++row;

        //q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
        A( row, col( i, j, q ) )     =   1;
        A( row, col( i, j + 1, v ) ) = - 6. * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 2, v ) ) =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 3, v ) ) = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q )
                   + iaR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j + 1, v )
                   - iaR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 2, v )
                   + iaR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, v );
        ++row;

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

        // Theta = Psi_eta  - Phi_zeta
        A( row, col( i, j, Theta ) )   =   1;
        A( row, col( i, j, Psi ) )     =   3 * Yd / ( 2 * dY );
        A( row, col( i, j + 1, Psi ) ) = - 4 * Yd / ( 2 * dY );
        A( row, col( i, j + 2, Psi ) ) =   1 * Yd / ( 2 * dY );
        A( row, col( i + 1, j, Phi ) ) =   1 * Xd / ( 2 * dX );
        A( row, col( i - 1, j, Phi ) ) = - 1 * Xd / ( 2 * dX );

        B[ row ] = - Q( i, j, Theta )
                   + Yd * (( -3. * Q( i, j, Psi ) + 4. * Q( i, j + 1, Psi )
                            - Q( i, j + 2, Psi ) ) / ( 2 * dY ))
                   - Xd * ( Q( i + 1, j, Phi ) - Q( i - 1, j, Phi ) ) / ( 2 * dX );
        ++row;

        // Main interior grid points
        for ( std::size_t j = 1; j < M; ++j )
        {
          // eta location
          double eta( ETA_NODES[ j ] );
          double Yd( SSI.mesh_Yd( eta ) );
          double Ydd( SSI.mesh_Ydd( eta ) );
          Vector<double> Base( base.get_interpolated_vars( eta ) );

          // Laplacian coefficients for finite-differencing

          // OS
          // X(i,j-1)
          double laplace_1_os =  ( Yd*Yd / (dY*dY) - Ydd / (2.*dY) ) ;
          // X(i-1,j)
          double laplace_3_os = ( Xd*Xd / (dX*dX) - Xdd / (2.*dX) ) / ( zeta0 * zeta0 );
          // X(i,j)
          double laplace_4_os = - 2 *( Yd*Yd / (dY*dY) + Xd*Xd / ( zeta0 * zeta0 * dX * dX ) )
                             - alpha * alpha;
          // X(i+1,j)
          double laplace_5_os = ( Xdd / (2.*dX) + Xd*Xd / (dX*dX) ) / ( zeta0 * zeta0 );
          // X(i,j+1)
          double laplace_7_os = ( Yd*Yd / (dY*dY) + Ydd / (2.*dY) );

          // Self-sim
          // X(i,j-1)
          double laplace_1 =  ( Yd*Yd/(dY*dY) - Ydd/ (2.*dY) ) ;
          // X(i-1,j)
          double laplace_3 = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) ) / ( zeta0 * zeta0 );
          // X(i,j)
          double laplace_4 = -2.*( Yd*Yd / (dY*dY) + Xd*Xd/( zeta0 * zeta0 * dX * dX ) );
          // X(i+1,j)
          double laplace_5 = ( Xdd/(2.*dX) + Xd*Xd/(dX*dX) ) / ( zeta0 * zeta0 );
          // X(i,j+1)
          double laplace_7 = ( Yd*Yd/(dY*dY) + Ydd/ (2.*dY) );

          // Guessed/known components and various derivative values
          Vector< std::complex<double> > Guess( Q.get_nodes_vars( i, j ) );
          Vector< std::complex<double> > Guess_eta( ( Q.get_nodes_vars( i, j + 1 )
                                    - Q.get_nodes_vars( i, j - 1 ) ) * ( Yd /( 2 * dY )) );
          Vector< std::complex<double> > Guess_hzeta( ( Q.get_nodes_vars( i + 1, j )
                                      - Q.get_nodes_vars( i - 1, j ) )
                                      * ( Xd /( 2 * dX )) );
          Vector< std::complex<double> > Guess_laplace( Q.get_nodes_vars( i, j - 1 ) * laplace_1
                                     +  Q.get_nodes_vars( i - 1, j ) * laplace_3
                                     +  Q.get_nodes_vars( i, j ) * laplace_4
                                     +  Q.get_nodes_vars( i + 1, j ) * laplace_5
                                     +  Q.get_nodes_vars( i, j + 1 ) * laplace_7 );

          Vector< std::complex<double> > Guess_eta_eta( ( Q.get_nodes_vars( i, j + 1 ) * laplace_7
                                                        + Q.get_nodes_vars( i, j ) * ( - 2. * Yd * Yd / ( dY * dY ) )
                                                        + Q.get_nodes_vars( i, j - 1 ) * laplace_1 ) );
          Vector< std::complex<double> > Guess_hzeta_hzeta( ( Q.get_nodes_vars( i + 1, j ) * laplace_5
                                                        + Q.get_nodes_vars( i, j ) * ( - 2. * Xd * Xd / ( zeta0 * zeta0 * dX * dX ))
                                                        + Q.get_nodes_vars( i - 1, j ) * laplace_3 ) );

          Vector< std::complex<double> > Guess_eta_hzeta( ( Q.get_nodes_vars( i + 1, j + 1 )
                                                          + Q.get_nodes_vars( i - 1, j - 1 )
                                                          - Q.get_nodes_vars( i + 1, j - 1 )
                                                          - Q.get_nodes_vars( i - 1, j + 1 ) ) * ( Xd * Yd / ( 4 * dX * dY) ) );

          double UBdd = beta * ( Base[ UB ] * Base[ UB ] - 1. ) - Base[ PhiB ] * Base[ UBd ];

          // Self similar solution
          //double U = sol( i, j, 6 ); // UB + U_pert
          //double U_hzeta = Xd * ( sol( i + 1, j, 2 ) - sol( i - 1, j, 2 ) ) / ( 2 * dX ); // U_pert_hzeta
          //double U_eta = Yd * ( sol( i, j + 1, 2 ) - sol( i, j - 1, 2 ) ) / ( 2 * dY ); // U_pert_eta
          //double UBd = base.get_interpolated_vars( eta )[ 1 ];
          //U_eta += UBd;

          //double UB = base.get_interpolated_vars( eta )[ 0 ];
          //double PhiB = base.get_interpolated_vars( eta )[ 2 ];
          //double U_eta_eta = laplace_1_os * sol( i, j - 1, 2 ) + laplace_7_os * sol( i, j + 1, 2 );
          //U_eta_eta += - 2 * Yd * Yd * sol( i, j, 2 ) / ( dY * dY );
          //U_eta_eta += beta * ( UB * UB - 1 ) - PhiB * UBd; // + UBdd = beta * [UB^2 - 1 ]- PhiB * UBd

          //double U_hzeta_hzeta = laplace_3_os * sol( i - 1, j, 2 ) + laplace_5_os * sol( i + 1, j, 2 );
          //U_hzeta_hzeta += - 2 * Xd * Xd * sol( i, j, 2 ) / ( dX * dX );

          //double U_eta_hzeta = Yd * Xd * ( sol( i + 1, j + 1, 2 ) + sol( i - 1, j - 1, 2 )
          //                               - sol( i + 1, j - 1, 2 ) - sol( i - 1, j + 1, 2 ) )
          //                               / ( 4 * dY * dX );

          ///////////////////////////////
          // 2D OrrSommerfeld equation //
          ///////////////////////////////

          ////////////////
          // v equation //
          ////////////////

          // ( i / alpha * Rx^(1/2) ) * Laplacian of v
          A( row, col( i, j - 1, v ) )     = iaR * laplace_1_os;
          A( row, col( i - 1, j, v ) )     = iaR * laplace_3_os;
          A( row, col( i, j, v ) )         = iaR * laplace_4_os;
          A( row, col( i + 1, j, v ) )     = iaR * laplace_5_os;
          A( row, col( i, j + 1, v ) )     = iaR * laplace_7_os;

          // + ( UB + UG - c_g ) * v
          A( row, col( i, j, v ) )        +=  Base[ UB ] + Guess[ U ] - c_g;

          // + ((1 - beta) / Rx) * v
          A( row, col( i, j, v ) )        +=  ( 1.0 - beta ) / Rx;

          // + v_g * U
          A( row, col( i, j, U ) )         =  Q( i, j, v);

          // - v_g * c
          A( row, size - 1 )               = - Q( i, j, v );

          // - q
          A( row, col( i, j, q ) )         = - 1.;

          B[ row ] =  - iaR * ( Q( i, j - 1, v ) * laplace_1_os
                             +  Q( i - 1, j, v ) * laplace_3_os
                             +  Q( i, j, v )     * laplace_4_os
                             +  Q( i + 1, j, v ) * laplace_5_os
                             +  Q( i, j + 1, v ) * laplace_7_os )
                      - ( Base[ UB ] + Guess[ U ] - c_g
                             + ( ( 1.0 - beta ) / Rx ) ) * Q( i, j, v )
                      + Q( i, j, q );
          ++row;

          ////////////////
          // w equation //
          ////////////////

          // ( i / alpha * Rx^(1/2) ) * Laplacian of v
          A( row, col( i, j - 1, w ) )     = iaR * laplace_1_os;
          A( row, col( i - 1, j, w ) )     = iaR * laplace_3_os;
          A( row, col( i, j, w ) )         = iaR * laplace_4_os;
          A( row, col( i + 1, j, w ) )     = iaR * laplace_5_os;
          A( row, col( i, j + 1, w ) )     = iaR * laplace_7_os;

          // + ( U - c_g ) * w
          A( row, col( i, j, w ) )        +=  Base[ UB ] + Guess[ U ] - c_g;

          // + ((1 - beta) / Rx) * w
          A( row, col( i, j, w ) )        +=  ( 1.0 - beta ) / Rx;

          // + w_g * U
          A( row, col( i, j, U ) )         =  Q( i, j, w);

          // - w_g * c
          A( row, size - 1 )               = - Q( i, j, w );

          // - s
          A( row, col( i, j, s ) )         = - 1.;

          B[ row ] =  - iaR * ( Q( i, j - 1, w ) * laplace_1_os
                             +  Q( i - 1, j, w ) * laplace_3_os
                             +  Q( i, j, w )     * laplace_4_os
                             +  Q( i + 1, j, w ) * laplace_5_os
                             +  Q( i, j + 1, w ) * laplace_7_os )
                      - ( Base[ UB ] + Guess[ U ] - c_g
                             + ( ( 1.0 - beta ) / Rx ) ) * Q( i, j, w )
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

          // + ( 1 / zeta0 ) * (UB' + UG_{eta}) * w_{hzeta}
          A( row, col( i + 1, j, w ) )      =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 );
          A( row, col( i - 1, j, w ) )      = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 );

          // + ( 1 / zeta0 ) * w_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      =   Guess_hzeta[ w ] * Yd / ( 2 * dY * zeta0 );
          A( row, col( i, j - 1, U ) )      = - Guess_hzeta[ w ] * Yd / ( 2 * dY * zeta0 );

          // - ( 1 / zeta0 ) * U_{hzeta} * w_{eta}
          A( row, col( i, j + 1, w ) )      = - Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 );
          A( row, col( i, j - 1, w ) )      =   Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 );

          // - ( 1 / zeta0 ) * w_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      = - Guess_eta[ w ] * Xd / ( 2 * dX * zeta0 );
          A( row, col( i - 1, j, U ) )      =   Guess_eta[ w ] * Xd / ( 2 * dX * zeta0 );

          // - ( 1 / zeta0 ) * UG_{hzeta eta} * w
          A( row, col( i, j, w ) )          = - Guess_eta_hzeta[ U ] / zeta0;

          // - w_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  =   Guess[ w ] * Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j - 1, U ) )  =   Guess[ w ] * Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i + 1, j - 1, U ) )  = - Guess[ w ] * Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j + 1, U ) )  = - Guess[ w ] * Xd * Yd / ( 4 * dX * dY * zeta0 );

          // - ( UB'' + UG_{eta eta} ) * v
          A( row, col( i, j, v ) )          = - ( UBdd + Guess_eta_eta[ U ] );

          // - v_g * U_{eta eta}
          A( row, col( i, j + 1, U ) )      = - Guess[ v ] * laplace_7;
          A( row, col( i, j, U ) )          =   Guess[ v ] * 2. * Yd * Yd / ( dY * dY );
          A( row, col( i, j - 1, U ) )      = - Guess[ v ] * laplace_1;

          B[ row ] = - ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) * Q( i, j - 1, q )
                     - ( - 2 * Yd * Yd / ( dY * dY ) ) * Q( i, j, q )
                     - ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) * Q( i, j + 1, q )
                     + alpha * alpha * Q( i, j, q )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j + 1, s )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j - 1, s )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j - 1, s )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j + 1, s )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 ) ) * Q( i + 1, j, w )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 ) ) * Q( i - 1, j, w )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j + 1, w )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j - 1, w )
                     + ( Guess_eta_hzeta[ U ] / zeta0 ) * Q( i, j, w )
                     + ( UBdd + Guess_eta_eta[ U ] ) * Q( i, j, v );
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

          // + ( 1 / zeta0 ) * UG_{hzeta} * v_{eta}
          A( row, col( i, j + 1, v ) )     =   Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 );
          A( row, col( i, j - 1, v ) )     = - Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 );

          // + ( 1 / zeta0 ) * v_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      = - Guess_eta[ v ] * Xd / ( 2 * dX * zeta0 );
          A( row, col( i - 1, j, U ) )      =   Guess_eta[ v ] * Xd / ( 2 * dX * zeta0 );

          // - ( 1 / zeta0 ) * (UBd + UG_{eta}) * v_{hzeta}
          A( row, col( i + 1, j, v ) )     = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 );
          A( row, col( i - 1, j, v ) )     =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 );

          // - ( 1 / zeta0 ) * v_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      = - Guess_hzeta[ v ] * Yd / ( 2 * dY * zeta0 );
          A( row, col( i, j - 1, U ) )      =   Guess_hzeta[ v ] * Yd / ( 2 * dY * zeta0 );

          // - ( 1 / zeta0 ) * UG_{eta hzeta} * v
          A( row, col( i, j, v ) )         = - Guess_eta_hzeta[ U ] / zeta0;

          // - v_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  =   Guess[ v ] * Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j - 1, U ) )  =   Guess[ v ] * Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i + 1, j - 1, U ) )  = - Guess[ v ] * Xd * Yd / ( 4 * dX * dY * zeta0 );
          A( row, col( i - 1, j + 1, U ) )  = - Guess[ v ] * Xd * Yd / ( 4 * dX * dY * zeta0 );

          // - ( 1 / zeta0^2 ) * UG_{hzeta hzeta} * w
          A( row, col( i, j, w ) )         = - Guess_eta_hzeta[ U ] / ( zeta0 * zeta0 );

          // - w_g * U_{hzeta hzeta} / zeta0^2
          A( row, col( i + 1, j, U ) )      = - Guess[ w ] * laplace_5;
          A( row, col( i, j, U ) )          =   Guess[ w ] * 2. * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
          A( row, col( i - 1, j, U ) )      = - Guess[ w ] * laplace_3;

          B[ row ] = - ( ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) ) / ( zeta0 * zeta0 ) ) * Q( i - 1, j, s )
                     - ( - 2 * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i, j, s )
                     - ( ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) ) / ( zeta0 * zeta0 ) ) * Q( i + 1, j, s )
                     + alpha * alpha * Q( i, j, s )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j + 1, q )
                     - ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j - 1, q )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i + 1, j - 1, q )
                     + ( Xd * Yd / ( 4 * dX * dY * zeta0 ) ) * Q( i - 1, j + 1, q )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j + 1, v )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * zeta0 ) ) * Q( i, j - 1, v )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 ) ) * Q( i + 1, j, v )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * zeta0 ) ) * Q( i - 1, j, v )
                     + ( Guess_eta_hzeta[ U ] / zeta0 ) * Q( i, j, v )
                     + ( Guess_hzeta_hzeta[ U ] / ( zeta0 * zeta0 ) ) * Q( i, j, w );
                     //TODO replace these with Guess[ ... ] instead of Q( i, j, ...)
          ++row;

          ///////////////////////////////
          // Self-sim streak equations //
          ///////////////////////////////

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
          A( row, col( i, j + 1, U ) )        = - ( 2. - beta ) * Yd / ( 2 * dY );
          A( row, col( i, j - 1, U ) )        =   ( 2. - beta ) * Yd / ( 2 * dY );
          // Theta_hzeta
          A( row, col( i + 1, j, Theta ) )    =   Xd / ( 2 * dX );
          A( row, col( i - 1, j, Theta ) )    = - Xd / ( 2 * dX );

          // Residual
          B[ row ]      = - Guess_laplace[ Phi ] + ( 2. - beta ) * Guess_eta[ U ]
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
          A( row, col( i + 1, j, U ) )        = - ( 2. - beta ) * Xd
                                                / ( 2. * dX * zeta0 * zeta0 );
          A( row, col( i - 1, j, U ) )        =   ( 2. - beta ) * Xd
                                                / ( 2. * dX * zeta0 * zeta0 );

          // -Theta_eta
          A( row, col( i, j + 1, Theta ) )    = - Yd / ( 2 * dY );
          A( row, col( i, j - 1, Theta ) )    =   Yd / ( 2 * dY );

          // Residual
          B[ row ]      = - Guess_laplace[ Psi ] + ( 2. - beta )
                          * ( Guess_hzeta[ U ] )
                          / ( zeta0 * zeta0 )
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
          A( row, col( i, j, U ) )           += - 2.* beta * ( Base[ UB ]
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
          A( row, col( i, j, Phi ) )          =   Base[ UBd ] + Guess_eta[ U ];

          //TODO Forcing (LHS)

          // Forcing (RHS)
          std::complex<double> F_1;
          F_1 = Guess[ v ] * ( conj( Guess_eta_eta[ v ] ) + conj( Guess_eta_hzeta[ w ] ) )
             + Guess[ w ] * ( conj( Guess_eta_hzeta[ v ] ) + conj( Guess_hzeta_hzeta[ w ] ) );
          F_1 *= std::complex<double> (0.0, - 1.0 / alpha);
          F_1 = F_1 + conj( F_1 );

          // Residual
          B[ row ]        = - Guess_laplace[ U ]
                            + beta * ( 2. * Base[ UB ] + Guess[ U ] ) * Guess[ U ]
                            - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                            * ( Guess_hzeta[ U ] )
                            - ( Base[ PhiB ] + Guess[ Phi ] ) * Guess_eta[ U ]
                            - Base[UBd] * Guess[Phi];
                            //+ std::pow( Rx, -2.0/3.0 ) * Sigma * Sigma * F_1;
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
          A( row, col( i, j + 1, U ) )         = - 2. * ( 1. - beta )
                                                   * ( Base[ UB ] + Guess[ U ] )
                                                   * ( hzeta )
                                                   * Yd / ( 2 * dY );
          A( row, col( i, j - 1, U ) )         =   2. * ( 1. - beta )
                                                   * ( Base[ UB ] + Guess[ U ] )
                                                   * ( hzeta )
                                                   * Yd / ( 2 * dY );

          // -2 * (1-beta) * (UB' + UG) * ( hzeta ) * U
          A( row, col( i, j, U ) )             = - 2. * ( 1. - beta )
                                                   * ( Base[ UBd ] + Guess_eta[ U ] )
                                                   * ( hzeta );

          // (2 * (1-beta) * eta * UG_hzeta / (zeta0^2)) * U
          A( row, col( i, j, U ) )            +=  2. * ( 1. - beta )
                                                  * eta * Guess_hzeta[ U ]
                                                  / ( zeta0 * zeta0 );

          // 2 * (1-beta) * eta * (UB + UG) * U_hzeta / ( zeta0^2 )
          A( row, col( i + 1, j, U ) )         =  2. * ( 1. - beta ) * eta
                                                  * ( Base[ UB ] + Guess[ U ] )
                                                  * Xd / ( 2 * dX * zeta0 * zeta0 );
          A( row, col( i - 1, j, U ) )         = -2. * ( 1. - beta ) * eta
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

          // [ThetaB + ThetaG_hzeta] * Psi
          A( row, col( i, j, Psi ) )           =  Base[ ThetaB ] + Guess_hzeta[ Theta ];

          // (2-beta) * ( UB + UG ) * Theta
          A( row, col( i, j, Theta ) )        +=   ( 2. - beta ) * ( Base[ UB ]
                                                 + Guess[ U ] );

          // (2-beta) * ( hzeta * ThetaB + ThetaG ) * U
          A( row, col( i, j, U ) )            +=   ( 2. - beta ) * ( hzeta * Base[ ThetaB ]
                                                 + Guess[ Theta ] );

          //TODO Forcing (LHS)

          // Forcing (RHS)
          std::complex<double> F_2;
          F_2 =   2. * ( Guess_hzeta[ v ] - Guess_eta[ w ] )
                     * ( conj( Guess_eta[ v ] ) + conj( Guess_hzeta[ w ] ) )
                + Guess[ v ] * ( 2. * conj( Guess_eta_hzeta[ v ] ) + conj( Guess_hzeta_hzeta[ w ] ) - conj( Guess_eta_eta[ w ] ) )
                + Guess[ w ] * ( conj( Guess_hzeta_hzeta[ v ] ) - conj( Guess_eta_eta[ v ] ) - 2. * conj( Guess_eta_hzeta[ w ] ) );
          F_2 = F_2 + conj( F_2 );

          // Residual
          B[ row ]      = - Guess_laplace[ Theta ]
                          + 2.*( 1. - beta )
                          * ( hzeta * ( Base[ UB ] + Guess[ U ] )
                          * Guess_eta[ U ] + hzeta * Base[ UBd ] * Guess[ U ]
                          - eta * ( Base[ UB ] + Guess[ U ] )
                          * ( Guess_hzeta[ U ] )
                          / ( zeta0 * zeta0 ) )
                          - ( Base[ PhiB ] + Guess[ Phi ] ) * Guess_eta[ Theta ]
                          - hzeta * Base[ ThetaBd ] * Guess[ Phi ]
                          - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                          * ( Guess_hzeta[ Theta ] ) - Guess[ Psi ] * Base[ ThetaB ]
                          - ( 2. - beta ) * ( ( Base[ UB ] + Guess[ U ] )
                          * Guess[ Theta ] + hzeta * Base[ ThetaB ] * Guess[ U ] );
                          //+ std::pow( Rx, -1.0/6.0 ) * Sigma * Sigma * F_2;
          ++row;

      }

      // eta = eta_inf boundary ( top boundary )
      j = M ;
      eta = ETA_NODES[ j ];
      Yd = SSI.mesh_Yd( eta );

      // v = 0
      A( row, col( i, j, v ) ) = 1;

      B[ row ] = - Q( i, j, v );
      ++row;

      // w = 0
      A( row, col( i, j, w ) ) = 1;

      B[ row ] = - Q( i, j, w );
      ++row;

      // s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
      A( row, col( i, j, s ) )        =   1;
      A( row, col( i, j - 1, w ) )    = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) );
      A( row, col( i, j - 2, w ) )    =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) );
      A( row, col( i, j - 3, w ) )    =   iaR * Yd * Yd / ( dY * dY );

      B[ row ] = - Q( i, j, s )
                 + iaR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) ) * Q( i, j - 1, w)
                 - iaR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) ) * Q( i, j - 2, w )
                 - iaR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, w );
      ++row;

      // q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
      A( row, col( i, j, q ) )        =   1;
      A( row, col( i, j - 1, v ) )    = - 6. * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, j - 2, v ) )    =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, j - 3, v ) )    = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

      B[ row ] = - Q( i, j, q )
                 + iaR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j - 1, v )
                 - iaR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 2, v )
                 + iaR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, v );
      ++row;

      const double rsq =  eta * eta + zeta0 * zeta0 * hzeta * hzeta ;
      // Phi_eta*( eta^2 + zeta_0^2*hzeta^2) + [ 2*eta - (eta^2 + zeta_0^2*hzeta^2)/eta ]*Phi = 0
      A( row, col( i, j, Phi ) )        =   3 * Yd * rsq / ( 2 * dY );
      A( row, col( i, j - 1, Phi ) )    = - 4 * Yd * rsq / ( 2 * dY );
      A( row, col( i, j - 2, Phi ) )    =   1 * Yd * rsq / ( 2 * dY );
      A( row, col( i, j, Phi ) )       +=   2 * eta - ( rsq / eta );

      B[ row ] = - ( 3. * Q( i, j, Phi ) - 4. * Q( i, j-1, Phi )
                     + Q( i, j-2, Phi ) ) * Yd * rsq / ( 2 * dY )
                 + ( ( rsq / eta ) - 2. * eta ) * Q( i, j, Phi );
      ++row;

      // Psi_eta*( eta^2 + zeta_0^2*hzeta^2) + 2 * eta * Psi = 0
      A( row, col( i, j, Psi ) )        =   3 * Yd * rsq / (2*dY);
      A( row, col( i, j - 1, Psi ) )    = - 4 * Yd * rsq / (2*dY);
      A( row, col( i, j - 2, Psi ) )    =   1 * Yd * rsq / (2*dY);
      A( row, col( i, j, Psi ) )       +=   2 * eta;

      B[ row ] =  - (( 3. * Q( i, j, Psi ) - 4. * Q( i, j-1, Psi )
                      + Q( i, j-2, Psi ) ) * Yd * rsq / ( 2 * dY ))
                      - 2. * eta * Q( i, j, Psi );
      ++row;

      // U = 0
      A( row, col( i, j, U ) )            =   1;

      B[ row ]                            = - Q( i, j, U );
      ++row;

      // Theta = 0
      A( row, col( i, j, Theta ) )        =   1;

      B[ row ]                            = - Q( i, j, Theta );
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

      const double rsq =  eta * eta + zeta0 * zeta0 * hzeta * hzeta ;

      // w = 0
      A( row, col( i, j, w ) ) = 1;

      B[ row ] = - Q( i, j, w );
      ++row;

      // v = 0
      A( row, col( i, j, v ) ) = 1;

      B[ row ] = - Q( i, j, v );
      ++row;

      // q - (1/zeta0^2) * (i / (alpha*Rx^(1/2))) * v_{hzeta hzeta} = 0
      A( row, col( i, j, q ) )        =   1;
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
      A( row, col( i, j, s ) )        =   1;
      A( row, col( i - 1, j, w ) )    = - 6. * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
      A( row, col( i - 2, j, w ) )    =   ( 3. / 2. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
      A( row, col( i - 3, j, w ) )    = - ( 2. / 9. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );

      B[ row ] = - Q( i, j, s )
                 + ( 6. * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i - 1, j, w )
                 - ( ( 3. / 2. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i - 2, j, w )
                 + ( ( 2. / 9. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 ) ) * Q( i - 3, j, w );
      ++row;

      // (eta^2 + zeta_0^2 * hzeta^2) * Phi_hzeta + 2 * zeta_0^2 * hzeta * Phi = 0
      A( row, col( i, j, Phi ) )          =   3 * Xd * rsq / ( 2 * dX );
      A( row, col( i - 1, j, Phi ) )      = - 4 * Xd * rsq / ( 2 * dX );
      A( row, col( i - 2, j, Phi ) )      =   1 * Xd * rsq / ( 2 * dX );
      A( row, col( i, j, Phi ) )         +=   2 * zeta0 * zeta0 * hzeta;

      B[ row ]        = - rsq * ( 3. * Q( i, j, Phi) - 4. * Q( i - 1, j, Phi)
                        + Q( i - 2, j, Phi) ) * Xd / ( 2 * dX )
                        - 2. * zeta0 * zeta0 * hzeta * Q( i, j, Phi );
      ++row;

      // (eta^2 + zeta_0^2 * hzeta^2)*Psi_hzeta + (2*zeta_0^2*hzeta-(eta^2 + zeta_0^2*hzeta^2)/hzeta)*Psi = 0
      A( row, col( i, j, Psi ) )          =   3 * Xd * rsq / ( 2 * dX );
      A( row, col( i - 1, j, Psi ) )      = - 4 * Xd * rsq / ( 2 * dX );
      A( row, col( i - 2, j, Psi ) )      =   1 * Xd * rsq / ( 2 * dX );
      A( row, col( i, j, Psi ) )         +=   2 * zeta0 * zeta0 * hzeta - (rsq / hzeta);


      B[ row ]  = - (rsq * ( 3. * Q( i, j, Psi ) - 4. * Q( i - 1, j, Psi )
                    + Q( i - 2, j, Psi) ) * Xd / ( 2 * dX ))
                    - 2. * zeta0 * zeta0 * hzeta  * Q( i, j, Psi)
                    + (rsq / hzeta)  * Q( i, j, Psi);
      ++row;

      // hzeta * U_hzeta + 2 * U = 0
      A( row, col( i, j, U ) )            =   hzeta * 3. * Xd / ( 2 * dX ) + 2.;
      A( row, col( i - 1, j, U ) )        = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, U ) )        =   hzeta * 1. * Xd / ( 2 * dX );

      B[ row  ] = - ( hzeta * Xd * ( 3. * Q( i, j, U ) - 4. * Q( i - 1, j, U )
                    + Q( i - 2, j, U) ) / ( 2 * dX ) ) - 2. * Q( i, j, U );
      ++row;

      // hzeta * Theta_hzeta + Theta = 0
      A( row, col( i, j, Theta )  )       =   hzeta * 3. * Xd / ( 2 * dX ) + 1.;
      A( row, col( i - 1, j, Theta ) )    = - hzeta * 4. * Xd / ( 2 * dX );
      A( row, col( i - 2, j, Theta ) )    =   hzeta * 1. * Xd / ( 2 * dX );

      B[ row ] = - ( hzeta * Xd * ( 3. * Q( i, j, Theta ) - 4. * Q( i - 1, j, Theta )
                        + Q( i - 2, j, Theta ) ) / ( 2 * dX ) ) - Q( i, j, Theta ) ;
      ++row;

    } // End of loop over RHS eta nodes


    // TODO extra condition for c ??? q(0,0) = 1 ???
    A( row, col( 0, 0, q ) ) = 1;
    B[ row ] = 1.0 - Q( 0, 0, q );
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
        Q( i, j, v )     += B[ col( i, j, v ) ];
        Q( i, j, w )     += B[ col( i, j, w ) ];
        Q( i, j, q )     += B[ col( i, j, q ) ];
        Q( i, j, s )     += B[ col( i, j, s ) ];
        Q( i, j, Phi )   += B[ col( i, j, Phi ) ];
        Q( i, j, Psi )   += B[ col( i, j, Psi ) ];
        Q( i, j, U )     += B[ col( i, j, U ) ];
        Q( i, j, Theta ) += B[ col( i, j, Theta ) ];
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

  Q.dump( "./DATA/VWI_local.dat" );

	cout << "FINISHED" << endl;

}
