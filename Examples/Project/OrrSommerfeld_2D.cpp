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
    // Rich's function
    //return - K * exp( - 0.1 * hzeta * hzeta ) * ( 1. - 2 * 0.1 * hzeta * hzeta );

  }
}; // End of class mySelfSimInjection

int main()
{
  cout << "*** ------- Solving the 2D OrrSommerfeld equation (EVP) ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 32.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 32.0 );           // Size of the domain in the eta direction
  const std::size_t N( 300 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 300 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 9.0 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.4 );              // Wavenumber (alpha hat)
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

  //TODO set an initial guess somehow = faster convergence ???
  orrsommerfeld_2D.set_initial_guess_from_evec( 0 );
  target = orrsommerfeld_2D.eigenvalues()[ 0 ];
  orrsommerfeld_2D.set_target( target );

  // Solve again (to see if it's any quicker)
  timer_OS.start();
  orrsommerfeld_2D.solve_evp();
  orrsommerfeld_2D.output();
  timer_OS.print();
  timer_OS.stop();

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
