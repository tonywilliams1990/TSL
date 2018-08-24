// Solve the VWI interaction equations
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"
#include "VWI.h"

enum{ v, w, q, s };
enum{ v_r, v_i, w_r, w_i, q_r, q_i, s_r, s_i, Phi, Psi, U, Theta };

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
  std::size_t N( 120 );             // Number of intervals in the zeta_hat direction
  std::size_t M( 120 );             // Number of intervals in the eta direction
  std::size_t MB( M * 100 );        // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 9.0 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.4 );              // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );         // Local Reynolds number
  double Sigma( 0.0 );              // Wave amplitude
  double tol( 1e-3 );               // Tolerance for c_i = 0
  bool read_from_file( true );

  double K_min( 0.0 );
  double K_step( 0.5 );
  double Sigma_step_initial( 1e-4 );
  double Sigma_step( Sigma_step_initial );

  std::complex<double> target( 0.75, 0.0 ); // Target for eigensolver

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

  Vector< std::complex<double> > eigenvalues;
  TwoD_node_mesh< std::complex<double> > eigenvectors;
  std::complex<double> c_guess;

if ( !read_from_file )
{
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
  eigenvalues = orrsommerfeld_2D.eigenvalues();
  eigenvectors = orrsommerfeld_2D.eigenvectors(); // v, w, q, s
  c_guess = eigenvalues[0];

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
}

  cout << "*** Solving VWI equations ( local ) ***" << endl;

  // Create the VWI object
  myVWI vwi( SSI, alpha, Rx, Sigma );
if ( !read_from_file )
{
  vwi.c_guess() = c_guess;

  // Set the initial guess
  TwoD_node_mesh<double> guess( X_NODES, Y_NODES, 12 );
  for ( std::size_t i = 0; i < X_NODES.size(); ++i )
  {
    for ( std::size_t j = 0; j < Y_NODES.size(); ++j )
    {
      guess( i, j, v_r )    = eigenvectors( i, j, v ).real();
      guess( i, j, v_i )    = eigenvectors( i, j, v ).imag();
      guess( i, j, w_r )    = eigenvectors( i, j, w ).real();
      guess( i, j, w_i )    = eigenvectors( i, j, w ).imag();
      guess( i, j, q_r )    = eigenvectors( i, j, q ).real();
      guess( i, j, q_i )    = eigenvectors( i, j, q ).imag();
      guess( i, j, s_r )    = eigenvectors( i, j, s ).real();
      guess( i, j, s_i )    = eigenvectors( i, j, s ).imag();
      guess( i, j, Phi )    = sol( i, j, 0 );
      guess( i, j, Psi )    = sol( i, j, 1 );
      guess( i, j, U )      = sol( i, j, 2 );
      guess( i, j, Theta )  = sol( i, j, 3 );
    }
  }
  vwi.set_guess( guess );
}

  //TwoD_node_mesh<double> solution;
  vwi.speed_up() = true;

  vwi.set_output_path();
  vwi.make_output_directory();
  vwi.make_eigenvalues_directory();

  double c_i;

  //do{
    do{

      vwi.solve_check_exists();
      //vwi.solve_local();
      //vwi.output();
      //vwi.output_eigenvalue();
      //solution = vwi.solution();
      //solution.dump_gnu("./DATA/VWI_local_output.dat");
      c_i = imag( vwi.c_guess() );

      cout << "  * K = " << vwi.injection() << endl;
      cout << "  * Rx^1/2 = " << sqrt( vwi.Reynolds() )
           << "  =>  Rx = " << vwi.Reynolds() << endl;
      cout << "  * Sigma = " << vwi.Sigma() << endl;
      cout << "  * c = " << vwi.c_guess() << endl;

      //vwi.Reynolds() -= 50;
      vwi.Sigma() += 0.1;
      //if ( c_i > 0 ){ vwi.injection() -= 0.1; }
      //if ( c_i <= 0 ){ vwi.Sigma() += 1e-4; }
      //vwi.injection() -= 0.1;

    //}while( c_i < 0.01 && c_i > 0.0 );
    //vwi.injection() -= 0.05;
  }while( vwi.injection() >= 0.0 );

  // mesh refinement
  /*cout << "  * Refining the mesh " << endl;
  vwi.solve_check_exists();
  N = 80;
  M = 80;
  MB = M * 100;
  vwi.refine_mesh( N, M, MB );
  vwi.output();
  vwi.output_eigenvalue();*/

	cout << "FINISHED" << endl;

}
