#include "Core"
#include "SelfSimInjection.h"

using namespace std;
using namespace TSL;

class mySelfSimInjection : public SelfSimInjection {
public:
  // Define the injection function
  double Phi_w_func( const double& hzeta ){
    // N_TRANSP and GAMMA defined the constructor in SelfSimInjection.h
    if ( N_TRANSP < 1 )
    {
      return 0.0;
    }
    else
    {
      double sum( 0.0 );
      int sign;
      for (std::size_t i=1; i<N_TRANSP; ++i)
      {
        sign = i % 2 ? -1 : 1; // equivalent to (-1)^i since i % 2 = 0 = false and i % 2 = 1 = true
        sum += sign * tanh( GAMMA * ( hzeta - ((1.*i)/N_TRANSP) ) );
      }
      sign = N_TRANSP % 2 ? -1 : 1; // (-1)^N
      return - K * 0.5 *( 1 + 2 * sum + sign * tanh( GAMMA * ( hzeta - 1. ) ) );
    }
  }
}; // End of class mySelfSimInjection

int main()
{
  // Define the domain + short scale injection parameters
  double hzeta_right( 16.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 128.0 );          // Size of the domain in the eta direction
  const std::size_t N( 200 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 200 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.0 );               // Hartree parameter
  double zeta0( 1.0 );              // Injection width ( initial )
  double K( 0.0 );                  // Injection parameter ( +ve = blowing )
  double K_max( 4.0 );              // Maximum value of the injection parameter
  double K_step( 0.1 );             // Increment when iterating

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


  cout << "*** ---------- Injection K step ---------- ***" << endl;
  cout << "  * We are solving using a " << N + 1 << " x " << M + 1
       << " mesh with zeta_hat_inf = " << hzeta_right << " and eta_inf = "
       << eta_top << "." << endl;

  // Solve the system and iterate on K
  try
  {
    SSI.iterate_on_K( K_step, K_max );
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  cout << "FINISHED" << endl;
}
