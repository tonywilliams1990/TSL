// Solve the 2D Orr-Sommerfeld equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"

using namespace std;
using namespace TSL;

enum {v, w, q, s};

class mySelfSimInjection : public SelfSimInjection {
public:
  // Define the injection function
  double Phi_w_func( const double& hzeta ){
    return - K * exp( - hzeta * hzeta );
    // Rich's function
    //return - K * exp( - 0.1 * hzeta * hzeta ) * ( 1. - 2 * 0.1 * hzeta * hzeta );

  }
}; // End of class mySelfSimInjection

int main()
{
  cout << "*** ------- Solving the 2D OrrSommerfeld equation (EVP) ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 20.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 20.0 );           // Size of the domain in the eta direction
  const std::size_t N( 180 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 180 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 11.2 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.8 );              // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );           // Local Reynolds number

  std::complex<double> target(0.76,0.0); // Target for eigensolver
  double alpha_max( 0.5 );
  double alpha_min( 0.01 );

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


  cout << "*** Solving the self similar injection base flow ***" << endl;
  Timer timer;
  timer.start();

  SSI.set_output_path();
  SSI.mesh_setup();
  SSI.solve_check_exists();
  //SSI.solve();

  // Setup the generalised eigenvalue problem A p = c B p (solved using SLEPc)
  cout << "*** Setting up the generalised eigenvalue problem ***" << endl;
  cout << "--- K = " << K << ", alpha = " << alpha << ", Rx^1/2 = " << sqrt(Rx) << endl;

  // Create the OrrSommerfeld_2D object
  std::size_t nev( 1 );
  OrrSommerfeld_2D orrsommerfeld_2D( SSI, alpha, Rx, nev );

  // Setup
  orrsommerfeld_2D.set_region(0.1,1.0,-1.0,1.0);
  orrsommerfeld_2D.set_target( target );
  orrsommerfeld_2D.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D.calc_eigenvectors() = true;

  // Solve
  Timer timer_OS;
  timer_OS.start();
  orrsommerfeld_2D.solve_evp();
  orrsommerfeld_2D.output();
  timer_OS.print();
  timer_OS.stop();

  TwoD_node_mesh< std::complex<double> > evecs;
  evecs = orrsommerfeld_2D.eigenvectors(); // v, w, q, s
  double norm;
  norm = real( evecs.square_integral2D( v ) + evecs.square_integral2D( w ) );
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

  OneD_node_mesh< std::complex<double> > evec_centreline( evecs.ynodes(), 4 );

  std::size_t i( 0 );
  for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
  {
    evec_centreline( j, v ) = evecs( i, j, v );
    evec_centreline( j, w ) = evecs( i, j, w );
    evec_centreline( j, q ) = evecs( i, j, q );
    evec_centreline( j, s ) = evecs( i, j, s );
  }

  evecs.dump_gnu("./DATA/OS2D_global.dat");

  cout << "int v on centreline = " << evec_centreline.integral2( v ) << endl;
  cout << "int w on centreline = " << evec_centreline.integral2( w ) << endl;
  cout << "int q on centreline = " << evec_centreline.integral2( q ) << endl;
  cout << "int s on centreline = " << evec_centreline.integral2( s ) << endl;

  std::complex<double> lambda;
  double eta( SSI.eta_nodes()[ 0 ] );
  double dY( SSI.y_nodes()[ 1 ] - SSI.y_nodes()[ 0 ] );
  double Yd( SSI.mesh_Yd( eta ) );
  lambda = ( 3 * Yd / ( 2 * dY ) ) * evecs( 0, 0, q )
           - ( 4 * Yd / ( 2 * dY ) ) * evecs( 0, 1, q )
           + ( 1 * Yd / ( 2 * dY ) ) * evecs( 0, 2, q );

  cout << "q_eta( 0, 0 ) = " << lambda << endl;

  std::complex<double> sum( 0.0 );
  double dx( 0.0 );
  // Sum
  for ( std::size_t j = 0; j < SSI.y_nodes().size() - 1; ++j )
  {
    //dx = ( NODES[ node + 1 ] - NODES[ node ] );
    eta = SSI.eta_nodes()[ j ];
    Yd = SSI.mesh_Yd( eta );
    sum += 0.5 * dY * ( evecs( 0, j, v ) + evecs( 0, j + 1, v ) ) / Yd;
  }
  cout << "sum = " << sum << endl;

  //TODO set an initial guess somehow = faster convergence ???
  //orrsommerfeld_2D.set_initial_guess_from_evec( 0 );
  //target = orrsommerfeld_2D.eigenvalues()[ 0 ];
  //orrsommerfeld_2D.set_target( target );

  // Solve again (to see if it's any quicker)
  //timer_OS.start();
  //orrsommerfeld_2D.solve_evp();
  //orrsommerfeld_2D.output();
  //timer_OS.print();
  //timer_OS.stop();

  //orrsommerfeld_2D.step_in_alpha( 0.01, alpha_max );

  // Step backwards in alpha
  /*OrrSommerfeld_2D orrsommerfeld_2D_back( SSI, alpha, Rx, nev );
  orrsommerfeld_2D_back.set_region(0.1,1.0,-1.0,1.0);
  orrsommerfeld_2D_back.set_target( target );
  orrsommerfeld_2D_back.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D_back.calc_eigenvectors() = true;
  orrsommerfeld_2D_back.step_back_in_alpha( 0.01, alpha_min );*/

  timer.print();
  timer.stop();

	cout << "FINISHED" << endl;

}
