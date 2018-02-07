// Solve the 2D Rayleigh pressure equation
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "Rayleigh_2D.h"

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
  cout << "*** ------- Approximating the critical alpha at which c_i=0 ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 30.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 30.0 );           // Size of the domain in the eta direction
  const std::size_t N( 300 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 300 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.0 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 3.0 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.4 );              // Wavenumber (alpha hat)
  double epsilon( 1e-3 );           // Tolerance (keep iterating while c_i > epsilon)

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

  // Setup the generalised eigenvalue problem A p = c B p (solved using SLEPc)
  cout << "*** Setting up the generalised eigenvalue problem ***" << endl;
  cout << "--- K = " << K << ", alpha = " << alpha << endl;

  // Create the Rayleigh_2D object
  std::size_t nev( 1 );
  Rayleigh_2D rayleigh_2D( SSI, alpha, nev );

  // Setup
  rayleigh_2D.set_region(0.1,1.0,-1.0,1.0);
  rayleigh_2D.set_target( std::complex<double>(0.45,0.05) );
  rayleigh_2D.set_order( "EPS_TARGET_IMAGINARY" );
  rayleigh_2D.calc_eigenvectors() = false;

  // Solve
  rayleigh_2D.iterate_to_neutral( epsilon );
  double alpha_crit( rayleigh_2D.alpha() );
  cout << "Critical wavenumber = " << alpha_crit << endl;

  timer.print();
  timer.stop();

	cout << "FINISHED" << endl;

}
