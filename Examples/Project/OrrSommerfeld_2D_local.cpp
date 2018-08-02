// Solve the 2D Orr-Sommerfeld equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"

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
  cout << "*** ------- 2D OrrSommerfeld equation (Local) ------- ***" << endl;

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

  Timer timer;
  timer.start();

  cout << "*** Solving the self-similar flow on a coarse mesh ***" << endl;
  cout << "  * N = " << N << ", M = " << M << endl;
  SSI.set_output_path();
  SSI.mesh_setup();
  SSI.solve_check_exists();

  cout << "*** Solving the 2D OrrSommerfeld equation (Global) ***" << endl;
  cout << "--- K = " << K << ", alpha = " << alpha << ", Rx^1/2 = " << sqrt(Rx) << endl;

  // Create the OrrSommerfeld_2D object
  std::size_t nev( 1 );
  OrrSommerfeld_2D orrsommerfeld_2D( SSI, alpha, Rx, nev );
  orrsommerfeld_2D.set_region(0.1,1.0,-1.0,1.0);
  orrsommerfeld_2D.set_target( target );
  orrsommerfeld_2D.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D.calc_eigenvectors() = true;
  Timer timer_OS;
  timer_OS.start();
  orrsommerfeld_2D.solve_evp();
  orrsommerfeld_2D.output();
  timer_OS.print();
  timer_OS.stop();

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
  SSI.set_output_path();
  SSI.mesh_setup();
  SSI.solve_check_exists();

  // Remesh the eigenvectors onto the refined mesh
  evecs.remesh1( SSI.hzeta_nodes(), SSI.eta_nodes() );

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
      evecs( i, j, v ) = evecs( i, j, v ) / sqrt(norm);
      evecs( i, j, w ) = evecs( i, j, w ) / sqrt(norm);
      evecs( i, j, q ) = evecs( i, j, q ) / sqrt(norm);
      evecs( i, j, s ) = evecs( i, j, s ) / sqrt(norm);
    }
  }



  // Solve the nonlinear system
  cout << "*** ------- Solving the 2D OrrSommerfeld equation (Local) ------- ***" << endl;

  // Current guess g
  TwoD_node_mesh< std::complex<double> > Q_out( SSI.hzeta_nodes(), SSI.eta_nodes(), 4 );
  std::complex<double> c_g( evals[0] );

  // Set the initial guess
  for ( std::size_t i = 0; i < N + 1; ++i )
  {
    for ( std::size_t j = 0; j < M + 1; ++j )
    {
      Q_out( i, j, v )  = evecs( i, j, v );
      Q_out( i, j, w )  = evecs( i, j, w );
      Q_out( i, j, q )  = evecs( i, j, q );
      Q_out( i, j, s )  = evecs( i, j, s );
    }
  }

  Q_out.dump_gnu("./DATA/OS2D_global.dat");

  // Solve
  orrsommerfeld_2D.update_SSI( SSI );
  TwoD_node_mesh<double> sol = SSI.solution();
  orrsommerfeld_2D.solve_local( c_g, Q_out, sol );

  cout << "  * c_refined = " << c_g << endl;
  Q_out.dump_gnu( "./DATA/OS2D_local.dat" );

  timer.print();
  timer.stop();

	cout << "FINISHED" << endl;
}
