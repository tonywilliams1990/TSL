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
  cout << "*** ------- Vortex-wave interaction ------- ***" << endl;

  // Define the domain + short scale injection parameters
  double hzeta_right( 32.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 32.0 );           // Size of the domain in the eta direction
  const std::size_t N( 300 );       // Number of intervals in the zeta_hat direction
  const std::size_t M( 300 );       // Number of intervals in the eta direction
  const std::size_t MB( M * 100 );  // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 7.5 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.1 );              // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );         // Local Reynolds number
  double Sigma( 0.0 );              // Wave amplitude
  double Relax( 1.0 );              // Relaxation parameter

  double K_min( 0.0 );
  double K_step( 0.1 );
  double Sigma_step( 0.05 );

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
  cout << "  * K = " << SSI.injection() << endl;
  SSI.solve();
  //SSI.output();
  //SSI.output_base_solution();
  sol = SSI.solution();
  base = SSI.base_flow_solution();
  cout << "  * zeta0 = " << SSI.injection_width() << ", A = " << SSI.mass_flux() << endl;

  //TwoD_node_mesh<double> new_sol( HZETA_NODES, ETA_NODES, 8 );
  //TwoD_node_mesh<double> diff( HZETA_NODES, ETA_NODES, 8);

  // Turn on forcing
  SSI.forcing( true );
  SSI.wave_amplitude() = 0.1;

  // Setup the stability equations
  std::size_t nev( 1 );
  OrrSommerfeld_2D orrsommerfeld_2D( SSI, alpha, Rx, nev );
  orrsommerfeld_2D.set_region(0.1,1.0,-1.0,1.0);
  std::complex<double> target(0.76,0.0);
  //orrsommerfeld_2D.set_target( std::complex<double>(0.76,0.0) );
  orrsommerfeld_2D.set_target( target );
  orrsommerfeld_2D.set_order( "EPS_TARGET_IMAGINARY" );
  orrsommerfeld_2D.calc_eigenvectors() = true;
  double c_i( 0.0 ); // Imaginary part of eigenvalue

  // Step in sigma
  /*do {
      // Solve the stability equations
      cout << "*** Solving the stability equations for v and w ***" << endl;
      orrsommerfeld_2D.update_SSI( SSI );
      orrsommerfeld_2D.solve_evp();

      // Return eigenvectors
      TwoD_node_mesh< std::complex<double> > evecs;
      evecs = orrsommerfeld_2D.eigenvectors(); // v, w, q, s

      // Put into separate v and w meshes
      TwoD_node_mesh< std::complex<double> > v( evecs.xnodes(), evecs.ynodes(), 1 );
      TwoD_node_mesh< std::complex<double> > w( evecs.xnodes(), evecs.ynodes(), 1 );

      for ( std::size_t i=0; i<evecs.xnodes().size(); ++i )
      {
        for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
        {
          v( i, j, 0 ) = evecs( i, j, 0 );
          w( i, j, 0 ) = evecs( i, j, 1 );
        }
      }

      // Normalise the eigenvectors
      double norm;
      norm= real( v.square_integral2D() + w.square_integral2D() );

      for ( std::size_t i=0; i<evecs.xnodes().size(); ++i )
      {
        for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
        {
          v( i, j, 0 ) = v( i, j, 0 ) / norm;
          w( i, j, 0 ) = w( i, j, 0 ) / norm;
        }
      }

      // Pass to SSI object to create forcing terms
      SSI.set_v_wave( v );
      SSI.set_w_wave( w );

      // Solve the streak equations (with forcing)
      cout << "*** Solving the streak equations (with forcing) ***" << endl;

      SSI.solve();
      new_sol = SSI.solution();
      cout << "  * A = " << SSI.mass_flux() << endl;
      sol = new_sol;

      c_i = orrsommerfeld_2D.eigenvalues()[0].imag();
      orrsommerfeld_2D.set_target( orrsommerfeld_2D.eigenvalues()[0] ); //TODO
      cout << "  * c_i = " << c_i << endl;
      cout << "  * Sigma = " << SSI.wave_amplitude() << endl;
      cout << "  * K = " << SSI.injection() << endl;
      // Decide how to vary K and Sigma and then resolve self-similar eqns

      SSI.wave_amplitude() += Sigma_step;
      cout << "*** Stepping in sigma ***" << endl;
      SSI.solve();

  }while( c_i < 0.0 );*/

  double c_i_minus, c_i_plus, Sigma_minus, Sigma_plus;

  do {

    // Iterate
    std::size_t iteration( 0 );
    std::size_t max_iterations( 20 );
    double c_i_old( 0.0 );
    double c_i_diff( 0.0 );

    TwoD_node_mesh< std::complex<double> > v( HZETA_NODES, ETA_NODES, 1 );
    TwoD_node_mesh< std::complex<double> > w( HZETA_NODES, ETA_NODES, 1 );
    TwoD_node_mesh< std::complex<double> > v_old( HZETA_NODES, ETA_NODES, 1 );
    TwoD_node_mesh< std::complex<double> > w_old( HZETA_NODES, ETA_NODES, 1 );
    TwoD_node_mesh< std::complex<double> > v_relax( HZETA_NODES, ETA_NODES, 1 );
    TwoD_node_mesh< std::complex<double> > w_relax( HZETA_NODES, ETA_NODES, 1 );

    do {

      // Solve the stability equations
      cout << "*** Solving the stability equations for v and w ***" << endl;
      orrsommerfeld_2D.update_SSI( SSI );
      Timer timer_OS;
      timer_OS.start();
      orrsommerfeld_2D.solve_evp();
      timer_OS.print();
      timer_OS.stop();

      c_i = orrsommerfeld_2D.eigenvalues()[0].imag();
      c_i_diff = c_i - c_i_old;
      c_i_old = c_i;

      // Return eigenvectors
      TwoD_node_mesh< std::complex<double> > evecs;
      evecs = orrsommerfeld_2D.eigenvectors(); // v, w, q, s

      // Put into separate v and w meshes
      for ( std::size_t i=0; i<evecs.xnodes().size(); ++i )
      {
        for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
        {
          v( i, j, 0 ) = evecs( i, j, 0 );
          w( i, j, 0 ) = evecs( i, j, 1 );
        }
      }

      // Normalise the eigenvectors
      double norm;
      norm= real( v.square_integral2D() + w.square_integral2D() );

      for ( std::size_t i=0; i<evecs.xnodes().size(); ++i )
      {
        for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
        {
          v( i, j, 0 ) = v( i, j, 0 ) / norm;
          w( i, j, 0 ) = w( i, j, 0 ) / norm;
        }
      }

      // Set relaxation wave perturbation
      for ( std::size_t i=0; i<evecs.xnodes().size(); ++i )
      {
        for ( std::size_t j=0; j<evecs.ynodes().size(); ++j )
        {
          v_relax( i, j, 0 ) = v_old( i, j, 0 ) + Relax * ( v( i, j, 0 ) - v_old( i, j, 0 ));
          w_relax( i, j, 0 ) = w_old( i, j, 0 ) + Relax * ( w( i, j, 0 ) - w_old( i, j, 0 ));
        }
      }

      // Pass to SSI object to create forcing terms
      SSI.set_v_wave( v_relax );
      SSI.set_w_wave( w_relax );

      // Update old v and w
      v_old = v;
      w_old = w;

      // Solve the streak equations (with forcing)
      cout << "*** Solving the streak equations (with forcing) ***" << endl;
      SSI.solve();
      cout << "  * A = " << SSI.mass_flux() << endl;
      cout << "  * iter = " << iteration << endl;
      cout << "  * c_i_diff = " << c_i_diff << endl;
      ++iteration;
    }while( ( std::abs( c_i_diff ) > 1e-4 ) && ( iteration < max_iterations ) );


    //TODO how to step in sigma to converge to c_i = 0?
    c_i = orrsommerfeld_2D.eigenvalues()[0].imag();
    cout << "  * c_i = " << c_i << endl;
    cout << "  * Sigma = " << SSI.wave_amplitude() << endl;
    // Decide how to vary Sigma and then resolve self-similar eqns
    if ( c_i > 0.0 )
    {
      c_i_plus = c_i;
      Sigma_plus = SSI.wave_amplitude();
      SSI.wave_amplitude() = ( c_i_plus * Sigma_minus + c_i_minus * Sigma_plus )
                            / ( c_i_plus - c_i_minus );
      Sigma_step *= 0.1;
      cout << "*** Stepping in sigma ***" << endl;
      SSI.solve();
    }
    else
    {
      c_i_minus = c_i;
      Sigma_minus = SSI.wave_amplitude();
      SSI.wave_amplitude() += Sigma_step;
      cout << "*** Stepping in sigma ***" << endl;
      SSI.solve();
    }

    /*cout << "  * c_i = " << c_i << endl;
    cout << "  * Sigma = " << SSI.wave_amplitude() << endl;

    // Newton step in Sigma ??
    cout << "*** Calculating Sigma step ***" << endl;
    double delta_Sigma( 1e-3 );
    SSI.wave_amplitude() += delta_Sigma;
    double c_i_star( 0.0 );
    SSI.solve();
    orrsommerfeld_2D.update_SSI( SSI );
    Timer timer_OS;
    timer_OS.start();
    orrsommerfeld_2D.solve_evp();
    timer_OS.print();
    timer_OS.stop();
    c_i_star = orrsommerfeld_2D.eigenvalues()[0].imag();
    Sigma_step = - c_i * delta_Sigma / ( c_i_star - c_i );
    cout << "  * Sigma_step = " << Sigma_step << endl;
    cout << "*** Stepping in sigma ***" << endl;
    SSI.wave_amplitude() += Sigma_step - delta_Sigma;
    SSI.solve();

    //TODO ???
    target.imag( Sigma_step * ( c_i_star - c_i ) / delta_Sigma );
    cout << "  * target = " << target << endl;
    orrsommerfeld_2D.set_target( target );*/

  }while( std::abs( c_i ) > 1e-3 );

  c_i = orrsommerfeld_2D.eigenvalues()[0].imag();
  cout << "  * c_i = " << c_i << endl;
  cout << "  * Sigma = " << SSI.wave_amplitude() << endl;
  cout << "  * K = " << SSI.injection() << endl;


  // Now output the solution (after solving once more)
  //SSI.set_output( true );
  //SSI.solve();
  //SSI.output();
  //SSI.output_base_solution();
  //orrsommerfeld_2D.update_SSI( SSI );
  //orrsommerfeld_2D.solve_evp();

  timer.print();
  timer.stop();

	cout << "FINISHED" << endl;

}
