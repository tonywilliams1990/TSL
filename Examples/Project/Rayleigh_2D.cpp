// Solve the 2D Rayleigh pressure equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"

using namespace std;
using namespace TSL;


int main()
{
  cout << "----- Solving the 2D Rayleigh equation (EVP) -----" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 30.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 30.0 );          // Size of the domain in the eta direction
  const std::size_t N( 200 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 200 );       // Number of intervals in the eta direction
  double beta( 0.0 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 8.0 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.0 );              // Wavenumber (alpha hat)

  // Solve the self similar injection flow
  SelfSimInjection SSI;
  SSI.hzeta_right() = hzeta_right;
  SSI.eta_top() = eta_top;
  SSI.hzeta_intervals() = N;
  SSI.eta_intervals() = M;
  SSI.hartree() = beta;
  SSI.injection_width() = zeta0;
  SSI.injection() = K;
  SSI.set_mesh( "UNIFORM" );
  SSI.set_base_flow( "2D" );
  SSI.speed_up( false );

  cout << "*** Solving the self similar injection base flow ***" << endl;
  SSI.solve();
  SSI.output();
  TwoD_node_mesh<double> sol = SSI.solution();

  OneD_node_mesh<double> base = SSI.base_flow_solution();
  cout << "  * zeta0 = " << SSI.injection_width() << ", A = " << SSI.mass_flux() << endl;

  Vector<double> ETA_NODES   = SSI.eta_nodes();
  Vector<double> HZETA_NODES = SSI.hzeta_nodes();
  Vector<double> X_NODES     = SSI.x_nodes();
  Vector<double> Y_NODES     = SSI.y_nodes();

  const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
  const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );

  //cout << "dX = " << dX << endl;
  //cout << "dY = " << dY << endl;

  // Read from file
  TwoD_node_mesh<double> Rich_sol( X_NODES, Y_NODES, 4 );
  cout << "--- Reading from file" << endl;

  ifstream inFile;
  inFile.open( "./DATA/Q_8_1.dat" );
  if (!inFile) {
    cerr << "Unable to open file.";
    exit(1);
  }
  double x, y;
  Vector<double> vals( 4 );

  for ( std::size_t i=0; i<X_NODES.size(); ++i )
  {
    for ( std::size_t j=0; j<Y_NODES.size(); ++j )
    {
      inFile >> x;
      inFile >> y;
      for ( std::size_t v=0; v<vals.size(); ++v )
      {
        inFile >> vals[v];
      }

      Rich_sol.set_nodes_vars( i, j, vals );
      //cout << "x = " << x << ", y = " << y << endl;
      //cout << "vals = " << vals << endl;
    }
  }

  inFile.close();

  //Rich_sol.read( "./DATA/Q_8_1.dat", true );
  cout << "--- Finished reading" << endl;

  // Check data read from file
  Rich_sol.dump_gnu( "Test_out.dat" );

  // Setup the generalised eigenvalue problem A p = c B p (solved using SLEPc)
  cout << "*** Setting up the generalised eigenvalue problem ***" << endl;

  SlepcInitialize(NULL,NULL,(char*)0,(char*)0);

  std::size_t N_hzeta = N + 1;
  std::size_t N_eta = M + 1;

  SparseMatrix< std::complex<double> > A( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta );
  SparseMatrix< std::complex<double> > B( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta );

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
    A( row, i * N_eta + j )      = -3 * Xd / ( 2 * dX );
    A( row, (i+1) * N_eta + j )  =  4 * Xd / ( 2 * dX );
    A( row, (i+2) * N_eta + j )  = -1 * Xd / ( 2 * dX );

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
      // Self similar solution ( U = UB + U_pert solution is stored in column 7 )
      // U_pert stored in column 3 -> index 2
      //double U = sol( i, j, 6 );
      //double U_eta = Yd * ( sol( i, j + 1, 2 ) - sol( i, j - 1, 2 ) ) / ( 2 * dY );
      //U_eta += base.get_interpolated_vars( eta )[ 1 ];
      //double U_hzeta = Xd * ( sol( i + 1, j, 6 ) - sol( i - 1, j, 6 ) ) / ( 2 * dX );
      //double U_hzeta = Xd * ( sol( i + 1, j, 2 ) - sol( i - 1, j, 2 ) ) / ( 2 * dX );

      // Use Rich's solution

      double U = Rich_sol( i, j, 2 );
      //double U_eta = Yd * ( Rich_sol( i, j + 1, 2 ) - Rich_sol( i, j - 1, 2 ) ) / ( 2 * dY );
      double U_hzeta = Xd * ( Rich_sol( i + 1, j, 2 ) - Rich_sol( i - 1, j, 2 ) ) / ( 2 * dX );

      double U_pert = U - base.get_interpolated_vars( eta )[ 0 ];
      double eta_jp1( ETA_NODES[ j + 1 ] );
      double eta_jm1( ETA_NODES[ j - 1 ] );
      double U_pert_jp1 = Rich_sol( i, j + 1, 2 ) - base.get_interpolated_vars( eta_jp1 )[ 0 ];
      double U_pert_jm1 = Rich_sol( i, j - 1, 2 ) - base.get_interpolated_vars( eta_jm1 )[ 0 ];

      double U_eta = Yd * ( U_pert_jp1 - U_pert_jm1 ) / ( 2 * dY );
      U_eta += base.get_interpolated_vars( eta )[ 1 ];


      /*double U = Rich_sol( i, j, 2 ) + base.get_interpolated_vars( eta )[ 0 ]; // + UB
      double U_eta = Yd * ( Rich_sol( i, j + 1, 2 ) - Rich_sol( i, j - 1, 2 ) ) / ( 2 * dY );
      U_eta += base.get_interpolated_vars( eta )[ 1 ]; // + UBd
      double U_hzeta = Xd * ( Rich_sol( i + 1, j, 2 ) - Rich_sol( i - 1, j, 2 ) ) / ( 2 * dX );*/


      // Laplacian coefficients for finite-differencing

      // X(i,j-1)
      double laplace_1 =  ( Yd*Yd/(dY*dY) - Ydd/ (2.*dY) ) ;
      // X(i-1,j)
      double laplace_3 = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) ) / ( zeta0 * zeta0 );
      // X(i,j)
      double laplace_4 = -2.*( Yd*Yd / (dY*dY)
                         + Xd*Xd/( zeta0 * zeta0 * dX * dX ) );
      // X(i+1,j)
      double laplace_5 = ( Xdd/(2.*dX) + Xd*Xd/(dX*dX) ) / ( zeta0 * zeta0 );
      // X(i,j+1)
      double laplace_7 = ( Yd*Yd/(dY*dY) + Ydd/ (2.*dY) );

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

      //cout << "2 * U_eta * Yd / ( 2 * dY ) = " << 2 * U_eta * Yd / ( 2 * dY ) << endl;

      // - 2 * U_hzeta * P_hzeta / zeta0^2
      A( row, (i + 1) * N_eta + j ) += -2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );
      A( row, (i - 1) * N_eta + j ) +=  2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );

      //cout << "2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 ) = " << 2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 ) << endl;

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
    //eta = ETA_NODES[ j ];
    //Yd = SSI.mesh_Yd( eta );

    // P = 0
    A( row, i * N_eta + j ) = 1;

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

  A.output( "./A_mat.dat", 4 );
  B.output( "./B_mat.dat", 4 );

  // Create the sparse eigenvalue problem
  Vector< std::complex<double> > lambdas;
  SparseEigenSystem< std::complex<double> > system( &A, &B );

  system.set_nev(1);
  system.set_region(0.1,1.0,-1.0,1.0);
  system.set_target( std::complex<double>(0.56,0.18) );
  system.set_order( "EPS_TARGET_IMAGINARY" );
  //system.set_order( "EPS_TARGET_MAGNITUDE" );

  try
  {
    //system.eigensolve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  lambdas = system.eigenvalues();
  cout << "lambdas = " << lambdas << endl;

  SlepcFinalize();

	cout << "FINISHED" << endl;

}
