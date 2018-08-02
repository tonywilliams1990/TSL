// Solve the 2D Orr-Sommerfeld equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"
#include "VWI.h"

enum{ v, w, q, s, Phi, Psi, U, Theta };

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

class myVWI : public VWI {
public:
  using VWI::VWI;
  // Define the injection function (should be same as above)
  double Phi_w_func( const double& hzeta ){
    return - K * exp( - hzeta * hzeta );
  }
}; // End of class myVWI

int main()
{
  cout << "*** ------- Solving the VWI equations ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 20.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 20.0 );           // Size of the domain in the eta direction
  std::size_t N( 130 );             // Number of intervals in the zeta_hat direction
  std::size_t M( 130 );             // Number of intervals in the eta direction
  std::size_t MB( M * 100 );        // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 8.5 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.4 );              // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );         // Local Reynolds number
  double Sigma( 0.1 );             // Wave amplitude
  double tol( 1e-3 );               // Tolerance for c_i = 0


  TSL::eta_intervals = M;

  double K_min( 0.0 );
  double K_step( 0.5 );
  double Sigma_step_initial( 1.0 );
  double Sigma_step( Sigma_step_initial );

  std::complex<double> target( 0.77, 0.0 ); // Target for eigensolver

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

  // Solve self-similar equations
  SSI.set_output_path();
  SSI.mesh_setup();
  SSI.solve_check_exists();

  Vector<double> ETA_NODES      = SSI.eta_nodes();
  Vector<double> HZETA_NODES    = SSI.hzeta_nodes();
  Vector<double> X_NODES        = SSI.x_nodes();
  Vector<double> Y_NODES        = SSI.y_nodes();
  Vector<double> BASE_ETA_NODES = SSI.base_eta_nodes();

  TwoD_node_mesh<double> sol( HZETA_NODES, ETA_NODES, 8 ); // Mesh for storing the solution
  OneD_node_mesh<double> base( BASE_ETA_NODES, 6 );
  sol = SSI.solution();
  base = SSI.base_flow_solution();

  cout << "*** Solving the 2D OrrSommerfeld equation ( global ) ***" << endl;
  cout << "--- K = " << K << ", alpha = " << alpha << ", Rx^1/2 = " << sqrt(Rx) << endl;

  // Create the OrrSommerfeld_2D object
  std::size_t nev( 1 );
  OrrSommerfeld_2D orrsommerfeld_2D( SSI, alpha, Rx, nev );
  orrsommerfeld_2D.set_region(0.1,1.0,-1.0,1.0);
  orrsommerfeld_2D.set_target( target );
  orrsommerfeld_2D.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D.calc_eigenvectors() = true;

  // Solve the global eigenvalue problem
  Timer timer_OS;
  timer_OS.start();
  orrsommerfeld_2D.solve_evp();
  orrsommerfeld_2D.output();
  timer_OS.print();
  timer_OS.stop();

  // Return the eigenvalues and eigenvectors
  Vector< std::complex<double> > eigenvalues;
  eigenvalues = orrsommerfeld_2D.eigenvalues();
  TwoD_node_mesh< std::complex<double> > eigenvectors;
  eigenvectors = orrsommerfeld_2D.eigenvectors(); // v, w, q, s
  std::complex<double> c_guess( eigenvalues[0] );

  // Normalise the eigenvectors
  double normalisation;
  normalisation = real( eigenvectors.square_integral2D( v )
                      + eigenvectors.square_integral2D( w ) );

  for ( std::size_t i = 0; i < eigenvectors.xnodes().size(); ++i )
  {
    for ( std::size_t j = 0; j < eigenvectors.ynodes().size(); ++j )
    {
      eigenvectors( i, j, v ) = eigenvectors( i, j, v ) / sqrt( normalisation );
      eigenvectors( i, j, w ) = eigenvectors( i, j, w ) / sqrt( normalisation );
      eigenvectors( i, j, q ) = eigenvectors( i, j, q ) / sqrt( normalisation );
      eigenvectors( i, j, s ) = eigenvectors( i, j, s ) / sqrt( normalisation );
    }
  }


  cout << "*** Solving VWI equations ( local ) ***" << endl;

  // Create the VWI object
  myVWI vwi( SSI, alpha, Rx, Sigma );
  vwi.c_guess() = c_guess;

  // Set the initial guess
  TwoD_node_mesh< std::complex<double> > guess( X_NODES, Y_NODES, 8 );
  for ( std::size_t i = 0; i < X_NODES.size(); ++i )
  {
    for ( std::size_t j = 0; j < Y_NODES.size(); ++j )
    {
      guess( i, j, v )      = eigenvectors( i, j, v );
      guess( i, j, w )      = eigenvectors( i, j, w );
      guess( i, j, q )      = eigenvectors( i, j, q );
      guess( i, j, s )      = eigenvectors( i, j, s );
      guess( i, j, Phi )    = sol( i, j, 0 );
      guess( i, j, Psi )    = sol( i, j, 1 );
      guess( i, j, U )      = sol( i, j, 2 );
      guess( i, j, Theta )  = sol( i, j, 3 );
    }
  }
  vwi.set_guess( guess );

  vwi.solve_local();

  TwoD_node_mesh< std::complex<double> > solution;
  solution = vwi.solution();
  solution.dump_gnu("./DATA/VWI_local_output.dat");

  cout << "  * c = " << vwi.c_guess() << endl;










	cout << "FINISHED" << endl;

}
