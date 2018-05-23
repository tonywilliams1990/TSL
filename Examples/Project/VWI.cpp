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
    //return - K * exp( - hzeta * hzeta );
    // Rich's function
    return - K * exp( - 0.1 * hzeta * hzeta ) * ( 1. - 2 * 0.1 * hzeta * hzeta );

  }
}; // End of class mySelfSimInjection

int main()
{
  cout << "*** ------- Vortex-wave interaction ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 30.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 30.0 );           // Size of the domain in the eta direction
  const std::size_t N( 200 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 200 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 8.0 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.05 );              // Wavenumber (alpha hat)
  double Rx( 500 * 500 );            // Local Reynolds number
  double Sigma( 0.0 );              // Wave amplitude

  double K_min( 8.0 );
  double K_step( 0.1 );
  double Sigma_step( 10.0 );

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
  SSI.wave_amplitude() = Sigma;
  SSI.local_Reynolds() = Rx;
  SSI.wavenumber() = alpha;
  SSI.set_output( false );

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

  /* Solve streak equations (no forcing) */
  cout << "*** Solving the streak equations (no forcing) ***" << endl;

  SSI.solve();
  //SSI.output();
  //SSI.output_base_solution();
  sol = SSI.solution();
  base = SSI.base_flow_solution();
  cout << "  * zeta0 = " << SSI.injection_width() << ", A = " << SSI.mass_flux() << endl;

  TwoD_node_mesh<double> new_sol( HZETA_NODES, ETA_NODES, 8 );
  TwoD_node_mesh<double> diff( HZETA_NODES, ETA_NODES, 8);

  // Turn on forcing
  SSI.forcing( true );
  SSI.wave_amplitude() = 1;

  /* Setup the stability equations */
  // Create the OrrSommerfeld_2D object
  std::size_t nev( 1 );
  OrrSommerfeld_2D orrsommerfeld_2D( SSI, alpha, Rx, nev );

  // Setup
  orrsommerfeld_2D.set_region(0.1,1.0,-1.0,1.0);
  orrsommerfeld_2D.set_target( std::complex<double>(0.46,0.01) );
  orrsommerfeld_2D.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D.calc_eigenvectors() = true;
  double c_i( 0.0 ); // Imaginary part of eigenvalue

  do {

    /* Iterate */
    double max_residual( 0.0 );
    std::size_t iteration( 0 );
    std::size_t max_iterations( 10 );

    do {

      /* Solve the stability equations */
      cout << "*** Solving the stability equations for v and w ***" << endl;
      orrsommerfeld_2D.update_SSI( SSI );
      orrsommerfeld_2D.solve_evp();

      // Return eigenvectors
      TwoD_node_mesh< std::complex<double> > evecs;
      evecs = orrsommerfeld_2D.eigenvectors(); // v, w, q, s

      // Put into separate v and w meshes
      TwoD_node_mesh< std::complex<double> > v( evecs.xnodes(), evecs.ynodes(), 1 );
      TwoD_node_mesh< std::complex<double> > w( evecs.xnodes(), evecs.ynodes(), 1 );

      //TODO normalise the eigenvectors ?

      for ( std::size_t i=0; i<evecs.xnodes().size(); ++i )
      {
        for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
        {
          v( i, j, 0 ) = evecs( i, j, 0 );
          w( i, j, 0 ) = evecs( i, j, 1 );
        }
      }

      // Pass to SSI object to create forcing terms
      SSI.set_v_wave( v );
      SSI.set_w_wave( w );

      /* Solve the streak equations (with forcing) */
      cout << "*** Solving the streak equations (with forcing) ***" << endl;

      SSI.solve();
      new_sol = SSI.solution();
      cout << "  * A = " << SSI.mass_flux() << endl;

      // Calculate the difference
      diff = sol - new_sol;
      //diff.dump( "./diff_dump.dat" );
      Vector<double> diff_vars;
      diff_vars = diff.get_vars();
      max_residual = diff_vars.norm_inf();
      cout << "  * max_residual = " << max_residual << endl;

      sol = new_sol;

      ++iteration;
    }while( ( max_residual > 1e-3 ) && ( iteration < max_iterations ) );

    c_i = orrsommerfeld_2D.eigenvalues()[0].imag();
    cout << "  * c_i = " << c_i << endl;
    cout << "  * Sigma = " << SSI.wave_amplitude() << endl;
    cout << "  * K = " << SSI.injection() << endl;
    // Decide how to vary K and Sigma and then resolve self-similar eqns
    if ( c_i > 0.01 )
    {
      SSI.injection() -= K_step;
      cout << "*** Stepping in K ***" << endl;
      SSI.solve();
    }
    else if ( c_i < 0.01 && c_i > 0.0 )
    {
      SSI.wave_amplitude() += Sigma_step;
      cout << "*** Stepping in sigma ***" << endl;
      SSI.solve();
    }
    else
    {
      break;
    }

  }while( SSI.injection() > K_min );

  // Now output the solution (after solving once more)
  SSI.set_output( true );
  SSI.solve();
  SSI.output();
  SSI.output_base_solution();

  timer.print();
  timer.stop();

	cout << "FINISHED" << endl;

}
