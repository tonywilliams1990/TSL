// Solve the 2D Orr-Sommerfeld equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"

enum{ Phi, Psi, U, Theta };
enum{ v, w, q, s };

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
  cout << "*** ------- Solving the VWI equations ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 20.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 20.0 );           // Size of the domain in the eta direction
  std::size_t N( 200 );             // Number of intervals in the zeta_hat direction
  std::size_t M( 200 );             // Number of intervals in the eta direction
  std::size_t MB( M * 100 );        // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 8.9 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.4 );              // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );         // Local Reynolds number
  double Sigma( 0.05 );              // Wave amplitude
  double tol( 1e-3 );               // Tolerance for c_i = 0
  double conv_tol( 1e-4 );          // Convergence tolerance for inner loop

  TSL::eta_intervals = M;

  double K_min( 0.0 );
  double K_step( 0.05 );
  double Sigma_step_initial( 0.05 );
  double Sigma_step( Sigma_step_initial );

  std::complex<double> target( 0.782, 0.0 ); // Target for eigensolver

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

  cout << "*** Solving the self-similar equations ( without forcing ) ***" << endl;
  cout << "  * N = " << N << ", M = " << M << endl;
  Timer timer;
  timer.start();

  // Solve self-similar equations on a coarse grid
  SSI.set_output_path();
  SSI.mesh_setup();

  Vector<double> ETA_NODES      = SSI.eta_nodes();
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


  /* Solve the 2D Orr-Sommerfeld equation (Global) */
  cout << "--- K = " << K << ", alpha = " << alpha << ", Rx^1/2 = " << sqrt(Rx) << endl;

  // Create the OrrSommerfeld_2D object
  std::size_t nev( 1 );
  OrrSommerfeld_2D orrsommerfeld_2D( SSI, alpha, Rx, nev );
  orrsommerfeld_2D.set_region(0.1,1.0,-1.0,1.0);
  orrsommerfeld_2D.set_target( target );
  orrsommerfeld_2D.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D.calc_eigenvectors() = true;
  double c_i( 0.0 ); // Imaginary part of eigenvalue
  double c_r( 0.0 ); // Real part of eigenvalue

  // Output K and Sigma to file
  TrackerFile metric( "./DATA/K_vs_Sigma_alpha_" + Utility::stringify( alpha, 3 ) + ".dat" );
  metric.push_ptr( &SSI.injection(), "K" );
  metric.push_ptr( &Sigma, "Sigma" );
  metric.push_ptr( &c_i, "c_i" );
  metric.header();

  do{

    double c_i_minus, c_i_plus, Sigma_minus, Sigma_plus;

    do{

        cout << "*** Solving the stability equations for v and w ***" << endl;
        orrsommerfeld_2D.update_SSI( SSI );
        Timer timer_OS;
        timer_OS.start();
        orrsommerfeld_2D.solve_evp();
        orrsommerfeld_2D.output();
        timer_OS.print();
        timer_OS.stop();

        c_i = orrsommerfeld_2D.eigenvalues()[0].imag();

        // Return the eigenvalues and eigenvectors
        Vector< std::complex<double> > evals;
        evals = orrsommerfeld_2D.eigenvalues();
        TwoD_node_mesh< std::complex<double> > evecs;
        evecs = orrsommerfeld_2D.eigenvectors(); // v, w, q, s

        // Normalise the eigenvectors
        double norm;
        norm = real( evecs.square_integral2D( v ) + evecs.square_integral2D( w ) );

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

        /* Resolve self-sim equations (with forcing) */
        cout << "*** Solving the self-similar equations ( with forcing ) ***" << endl;

        // Current guess g
        TwoD_node_mesh<double> Q( X_NODES, Y_NODES, 4 );
        // Set the initial guess
        for ( std::size_t i = 0; i < N + 1; ++i )
        {
          for ( std::size_t j = 0; j < M + 1; ++j )
          {
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
        std::size_t size( 4 * ( M + 1 ) * ( N + 1 ) );
        Vector<double> B( size, 0.0 );
        // Iterate to a solution
        double max_residual( 0.0 );             // Maximum residual
        std::size_t iteration( 0 );             // Initialise iteration counter
        std::size_t max_iterations( 20 );       // Maximum number of iterations

        do {
            SparseMatrix<double> A( size, size );
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

                Vector<double> Guess_eta_eta( ( Q.get_nodes_vars( i, j + 1 ) * laplace_7
                                                              + Q.get_nodes_vars( i, j ) * ( - 2. * Yd * Yd / ( dY * dY ) )
                                                              + Q.get_nodes_vars( i, j - 1 ) * laplace_1 ) );
                Vector<double> Guess_hzeta_hzeta( ( Q.get_nodes_vars( i + 1, j ) * laplace_5
                                                              + Q.get_nodes_vars( i, j ) * ( - 2. * Xd * Xd / ( zeta0 * zeta0 * dX * dX ))
                                                              + Q.get_nodes_vars( i - 1, j ) * laplace_3 ) );

                Vector<double> Guess_eta_hzeta( ( Q.get_nodes_vars( i + 1, j + 1 )
                                                                + Q.get_nodes_vars( i - 1, j - 1 )
                                                                - Q.get_nodes_vars( i + 1, j - 1 )
                                                                - Q.get_nodes_vars( i - 1, j + 1 ) ) * ( Xd * Yd / ( 4 * dX * dY) ) );

                double UBdd = beta * ( Base[ UB ] * Base[ UB ] - 1. ) - Base[ PhiB ] * Base[ UBd ];

                // Wave eigenvectors for forcing terms
                Vector< std::complex<double> > Wave( evecs.get_nodes_vars( i, j ) );
                Vector< std::complex<double> > Wave_eta( ( evecs.get_nodes_vars( i, j + 1 )
                                          - evecs.get_nodes_vars( i, j - 1 ) ) * ( Yd /( 2 * dY )) );
                Vector< std::complex<double> > Wave_hzeta( ( evecs.get_nodes_vars( i + 1, j )
                                            - evecs.get_nodes_vars( i - 1, j ) )
                                            * ( Xd /( 2 * dX )) );
                Vector< std::complex<double> > Wave_eta_eta( ( evecs.get_nodes_vars( i, j + 1 ) * laplace_7
                                              + evecs.get_nodes_vars( i, j ) * ( - 2. * Yd * Yd / ( dY * dY ) )
                                              + evecs.get_nodes_vars( i, j - 1 ) * laplace_1 ) );
                Vector< std::complex<double> > Wave_hzeta_hzeta( ( evecs.get_nodes_vars( i + 1, j ) * laplace_5
                                              + evecs.get_nodes_vars( i, j ) * ( - 2. * Xd * Xd / ( zeta0 * zeta0 * dX * dX ))
                                              + evecs.get_nodes_vars( i - 1, j ) * laplace_3 ) );

                Vector< std::complex<double> > Wave_eta_hzeta( ( evecs.get_nodes_vars( i + 1, j + 1 )
                                              + evecs.get_nodes_vars( i - 1, j - 1 )
                                              - evecs.get_nodes_vars( i + 1, j - 1 )
                                              - evecs.get_nodes_vars( i - 1, j + 1 ) ) * ( Xd * Yd / ( 4 * dX * dY) ) );

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
                /*std::complex<double> F_1;
                F_1 = Guess[ v ] * ( conj( Guess_eta_eta[ v ] ) + conj( Guess_eta_hzeta[ w ] ) )
                   + Guess[ w ] * ( conj( Guess_eta_hzeta[ v ] ) + conj( Guess_hzeta_hzeta[ w ] ) );
                F_1 *= std::complex<double> (0.0, - 1.0 / alpha);
                F_1 = F_1 + conj( F_1 );*/

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

                // Forcing
                std::complex<double> F_2;
                F_2 =   2. * ( Wave_hzeta[ v ] - Wave_eta[ w ] )
                           * ( conj( Wave_eta[ v ] ) + conj( Wave_hzeta[ w ] ) )
                      + Wave[ v ] * ( 2. * conj( Wave_eta_hzeta[ v ] ) + conj( Wave_hzeta_hzeta[ w ] ) - conj( Wave_eta_eta[ w ] ) )
                      + Wave[ w ] * ( conj( Wave_hzeta_hzeta[ v ] ) - conj( Wave_eta_eta[ v ] ) - 2. * conj( Wave_eta_hzeta[ w ] ) );
                F_2 = F_2 + conj( F_2 );

                double Forcing = real( F_2 );

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
                                * Guess[ Theta ] + hzeta * Base[ ThetaB ] * Guess[ U ] )
                                + std::pow( Rx, -1.0/6.0 ) * Sigma * Sigma * Forcing;
                ++row;

            }

            // eta = eta_inf boundary ( top boundary )
            j = M ;
            eta = ETA_NODES[ j ];
            Yd = SSI.mesh_Yd( eta );

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

          max_residual = B.norm_inf();
          std::cout << "***                                              Maximum residual = "
                    << max_residual << std::endl;

          // Solve the sparse system
          Vector<double> x;
          x = A.solve( B );
          B = x;

          timer.print();
          timer.stop();

          // Update the known values using the correction which we just found
          for ( std::size_t i = 0; i < N + 1; ++i )
          {
            for ( std::size_t j = 0; j < M + 1; ++j )
            {
              Q( i, j, Phi )   += B[ col( i, j, Phi ) ];
              Q( i, j, Psi )   += B[ col( i, j, Psi ) ];
              Q( i, j, U )     += B[ col( i, j, U ) ];
              Q( i, j, Theta ) += B[ col( i, j, Theta ) ];
            }
          }

          std::cout << "***    Iteration = " << iteration
                    << "    Maximum correction = " << B.norm_inf() << std::endl;

          ++iteration;
        }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) );


        for ( std::size_t i = 0; i < N + 1; ++i )
        {
          double hzeta=HZETA_NODES[i];
          for ( std::size_t j = 0; j < M + 1; ++j )
          {
            double eta=ETA_NODES[j];
            sol( i, j, 0 ) = Q( i, j, Phi);
            sol( i, j, 1 ) = Q( i, j, Psi);
            sol( i, j, 2 ) = Q( i, j, U);
            sol( i, j, 3 ) = Q( i, j, Theta);
            sol( i, j, 4 ) =   Q( i, j, Phi)
                                  + base.get_interpolated_vars( eta )[ PhiB ];
            sol( i, j, 5 ) =   Q( i, j, Psi)
                                  + hzeta * base.get_interpolated_vars( eta )[ PsiB ];
            sol( i, j, 6 ) =   Q( i, j, U)
                                  + base.get_interpolated_vars( eta )[ UB ];
            sol( i, j, 7 ) =   Q( i, j, Theta)
                                  + base.get_interpolated_vars( eta )[ ThetaB ];
            }
        }

        SSI.set_solution( sol );

      c_i = orrsommerfeld_2D.eigenvalues()[0].imag();
      c_r = orrsommerfeld_2D.eigenvalues()[0].real();

      // Decide how to vary Sigma
      if ( c_i > 0.0 && std::abs( c_i ) > tol )
      {
        cout << "  * c_i = " << c_i << endl;
        cout << "  * Sigma = " << Sigma << endl;
        c_i_plus = c_i;
        Sigma_plus = Sigma;
        Sigma = Sigma_minus
           + ( c_i_minus * ( Sigma_minus - Sigma_plus ) / ( c_i_plus - c_i_minus ) );
        Sigma_step *= 0.1;
        cout << "*** Stepping in sigma ***" << endl;
        cout << "  * Step = " << Sigma_step << endl;
      }
      else if ( c_i < 0.0 && std::abs( c_i ) > tol )
      {
        cout << "  * c_i = " << c_i << endl;
        cout << "  * Sigma = " << Sigma << endl;
        c_i_minus = c_i;
        Sigma_minus = Sigma;
        Sigma += Sigma_step;
        cout << "*** Stepping in sigma ***" << endl;
        cout << "  * Step = " << Sigma_step << endl;
      }

      if ( c_i < 0.0 )
      {
        target.imag( c_i + 0.005 );
      } else {
        target.imag( 0.0 );
      }
      target.real( c_r );
      orrsommerfeld_2D.set_target( target );
      cout << "  * target = " << target << endl;

    }while( std::abs( c_i ) > tol );

    c_i = orrsommerfeld_2D.eigenvalues()[0].imag();
    cout << "  * c_i = " << c_i << endl;
    cout << "  * Sigma = " << Sigma << endl;
    cout << "  * K = " << SSI.injection() << endl;
    metric.update();
    sol.dump( "./DATA/VWI_try_K_" + Utility::stringify( SSI.injection(), 3 )
              + "_alpha_" + Utility::stringify( alpha, 3 ) + ".dat" );
    // Step in K
    SSI.injection() -= K_step;
    cout << "*** Stepping in K ***" << endl;
    Sigma_step = Sigma_step_initial; // Reset sigma step

  }while( SSI.injection() >= K_min );

  timer.print();
  timer.stop();

  sol.dump( "./DATA/VWI_try_K_" + Utility::stringify( SSI.injection(), 3 )
            + "_alpha_" + Utility::stringify( alpha, 3 ) + ".dat" );

	cout << "FINISHED" << endl;

}
