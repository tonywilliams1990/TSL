#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <sstream>

#include "Core"
#include "SelfSimInjection.h"

using namespace std;
using namespace TSL;

class mySelfSimInjection : public SelfSimInjection {
public:
  // Define the injection function
  double Phi_w_func( const double& hzeta ){
    // Top-hat injection
    //return - K * 0.5 * ( tanh( GAMMA * ( hzeta - 1. ) )
    //       - tanh( GAMMA * ( hzeta - 2. ) ) );
    // Gaussian
    return - K * exp( - hzeta * hzeta );
  }
}; // End of class mySelfSimInjection

int main()
{
  // Define the domain + short scale injection parameters
  double hzeta_right( 16.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 64.0 );          // Size of the domain in the eta direction
  const std::size_t N( 200 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 200 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.0 );               // Hartree parameter
  double zeta0( 1.0 );              // Injection width ( initial )
  double K( 2.5 );                  // Injection parameter ( +ve = blowing )
  double zeta0_max( 1.0 );          // Maximum value of the injection width
  double zeta0_step( 1.0 );         // Increment when iterating

  // Setup the problem
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

  cout << "*** ---------- Injection with pressure gradient ---------- ***" << endl;
  cout << "  * We are solving using a " << N + 1 << " x " << M + 1
       << " mesh with zeta_hat_inf = " << hzeta_right << " and eta_inf = "
       << eta_top << "." << endl;

  // Solve the system and iterate on zeta0
  try
  {
    SSI.iterate_on_zeta0( zeta0_step, zeta0_max );
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  cout << "FINISHED" << endl;
}
