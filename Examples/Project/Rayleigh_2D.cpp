// Solve the 2D Rayleigh pressure equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"

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
  cout << "*** ------- Solving the 2D Rayleigh equation (EVP) ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 30.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 30.0 );           // Size of the domain in the eta direction
  const std::size_t N( 100 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 100 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.0 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 4.0 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.3 );              // Wavenumber (alpha hat)

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
  SSI.set_mesh( "NONUNIFORM" );
  SSI.set_base_flow( "2D" );
  SSI.speed_up( false );

  cout << "*** Solving the self similar injection base flow ***" << endl;
  Timer timer;
  timer.start();

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
  }

  const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
  const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );

  // Setup the generalised eigenvalue problem A p = c B p (solved using SLEPc)
  cout << "*** Setting up the generalised eigenvalue problem ***" << endl;
  cout << "--- K = " << K << ", alpha = " << alpha << endl;
  SlepcInitialize(NULL,NULL,(char*)0,(char*)0);

  std::size_t N_hzeta = N + 1;
  std::size_t N_eta = M + 1;

  SparseMatrix< std::complex<double> > A( N_eta * N_hzeta, N_eta * N_hzeta );
  SparseMatrix< std::complex<double> > B( N_eta * N_hzeta, N_eta * N_hzeta );

  SparseEigenSystem< std::complex<double> > system( &A, &B );
  system.set_nev(5);
  system.set_region(0.1,1.0,-1.0,1.0);
  system.set_target( std::complex<double>(0.56,0.18) );
  //system.set_target( std::complex<double>(0.45,0.001) );
  system.set_order( "EPS_TARGET_IMAGINARY" );
  system.calc_eigenvectors() = true;

  // Fill the sparse matrices
  std::size_t row( 0 );                               // Initialise row counter

  // hzeta = 0 boundary ( left boundary )
  std::size_t i( 0 );

  for ( std::size_t j = 0; j < M + 1 ; ++j )
  {
    double hzeta( HZETA_NODES[ 0 ] );
    double Xd( SSI.mesh_Xd( hzeta ) );
    double eta( ETA_NODES[ j ] );

    // P_hzeta = 0
    A( row, i * N_eta + j )           = -3 * Xd / ( 2 * dX );
    A( row, ( i + 1 ) * N_eta + j )   =  4 * Xd / ( 2 * dX );
    A( row, ( i + 2 ) * N_eta + j )   = -1 * Xd / ( 2 * dX );

    ++row;

  } // end for loop over LHS eta nodes

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

    // P_eta = 0
    A( row, i * N_eta + j )        = -3 * Yd / ( 2 * dY );
    A( row, i * N_eta + j + 1 )    =  4 * Yd / ( 2 * dY );
    A( row, i * N_eta + j + 2 )    = -1 * Yd / ( 2 * dY );

    ++row;

    // Main interior grid points
    for ( std::size_t j = 1; j < M; ++j )
    {
      // eta location
      double eta( ETA_NODES[ j ] );
      double Yd( SSI.mesh_Yd( eta ) );
      double Ydd( SSI.mesh_Ydd( eta ) );

      // Self similar solution
      double U = sol( i, j, 6 ); // UB + U_pert
      double U_hzeta = Xd * ( sol( i + 1, j, 2 ) - sol( i - 1, j, 2 ) ) / ( 2 * dX ); // U_pert_hzeta
      double U_eta = Yd * ( sol( i, j + 1, 2 ) - sol( i, j - 1, 2 ) ) / ( 2 * dY ); // U_pert_eta
      U_eta += base.get_interpolated_vars( eta )[ 1 ]; // + UBd

      // Laplacian coefficients for finite-differencing

      // X(i,j-1)
      double laplace_1 =  ( Yd*Yd / (dY*dY) - Ydd / (2.*dY) ) ;
      // X(i-1,j)
      double laplace_3 = ( Xd*Xd / (dX*dX) - Xdd / (2.*dX) ) / ( zeta0 * zeta0 );
      // X(i,j)
      double laplace_4 = -2.*( Yd*Yd / (dY*dY)
                         + Xd*Xd / ( zeta0 * zeta0 * dX * dX ) );
      // X(i+1,j)
      double laplace_5 = ( Xdd / (2.*dX) + Xd*Xd / (dX*dX) ) / ( zeta0 * zeta0 );
      // X(i,j+1)
      double laplace_7 = ( Yd*Yd / (dY*dY) + Ydd / (2.*dY) );

      //////////////////////////
      // 2D Rayleigh equation //
      //////////////////////////

      // A matrix

      // Laplacian of P * U
      A( row, i * N_eta + j - 1 )    = laplace_1 * U;
      A( row, (i - 1) * N_eta + j )  = laplace_3 * U;
      A( row, i * N_eta + j )        = laplace_4 * U;
      A( row, (i + 1) * N_eta + j )  = laplace_5 * U;
      A( row, i * N_eta + j + 1 )    = laplace_7 * U;

      // - alpha^2 * U * P
      A( row, i * N_eta + j )       += - alpha * alpha * U;

      // - 2 * U_eta * P_eta
      A( row, i * N_eta + j + 1 )   += -2 * U_eta * ( Yd / ( 2 * dY ) );
      A( row, i * N_eta + j - 1 )   +=  2 * U_eta * ( Yd / ( 2 * dY ) );

      // - 2 * U_hzeta * P_hzeta / zeta0^2
      A( row, (i + 1) * N_eta + j ) += -2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );
      A( row, (i - 1) * N_eta + j ) +=  2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );

      // B matrix

      // Laplacian of P
      B( row, i * N_eta + j - 1 )    = laplace_1;
      B( row, (i - 1) * N_eta + j )  = laplace_3;
      B( row, i * N_eta + j )        = laplace_4;
      B( row, (i + 1) * N_eta + j )  = laplace_5;
      B( row, i * N_eta + j + 1 )    = laplace_7;

      // - alpha^2 * P
      B( row, i * N_eta + j )       += - alpha * alpha;

      ++row;

    }

    // eta = eta_inf boundary ( top boundary )
    j = M ;
    eta = ETA_NODES[ j ];
    Yd = SSI.mesh_Yd( eta );

    // P = 0 (Dirichlet)
    //A( row, i * N_eta + j ) = 1;

    // Robin condition (P_eta + alpha*P = 0)
    // P_eta
    A( row, i * N_eta + j )        =  3 * Yd / ( 2 * dY );
    A( row, i * N_eta + j - 1 )    = -4 * Yd / ( 2 * dY );
    A( row, i * N_eta + j - 2 )    =  1 * Yd / ( 2 * dY );

    // + alpha * P
    A( row, i * N_eta + j )       += alpha;

    ++row;

  } // End of for loop over interior nodes

  // hzeta = hzeta_inf boundary ( right boundary )

  for ( std::size_t j = 0; j < M + 1; ++j )
  {
    std::size_t i( N );
    // P = 0
    A( row, i * N_eta + j ) = 1;
    ++row;

  } // End of loop over nodes

  //A.output( "./A_mat.dat", 4 );
  //B.output( "./B_mat.dat", 4 );

  // Solve the sparse eigenvalue problem
  try
  {
    system.eigensolve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED (EIGENSYS) \033[0m\n";
    assert( false );
  }

  Vector< std::complex<double> > lambdas;
  lambdas = system.eigenvalues();
  //cout << "lambdas = " << lambdas << endl;
  lambdas.output( SSI.output_path() + "alpha_"
                + Utility::stringify( alpha, 3 ) + "_evals.dat", 15 );


  Matrix< std::complex<double> > evecs;
  evecs = system.eigenvectors();

  std::size_t nev( system.get_nconv() ); // Number of converged eigenvalues

  // Convert matrix to 2D node mesh (complex)
  TwoD_node_mesh< std::complex<double> > output( HZETA_NODES, ETA_NODES, nev );
  // Push the data back into the unmapped domain
  for ( std::size_t n = 0; n < nev; ++n )
  {
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      for ( std::size_t j = 0; j < N_eta; ++j )
      {
        output( i, j, n ) = evecs( n, i * N_eta + j );
      }
    }
  }

  // Output mesh to file
  output.dump_gnu( SSI.output_path() + "alpha_"
                 + Utility::stringify( alpha, 3 ) + "_evecs.dat" );

  SlepcFinalize();

  timer.print();
  timer.stop();

	cout << "FINISHED" << endl;

}
