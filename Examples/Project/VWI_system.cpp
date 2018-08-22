// Solve the VWI interaction equations
#include <cassert>
#include <fstream>

#include "Core"
#include "Eigenvalue"
#include "SelfSimInjection.h"
#include "OrrSommerfeld_2D.h"
#include "VWI.h"

enum{ v, w, q, s };

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

  double hzeta_right( 20.0 );       // Size of the domain in the zeta_hat direction
  double eta_top( 20.0 );           // Size of the domain in the eta direction
  std::size_t N( 200 );             // Number of intervals in the zeta_hat direction
  std::size_t M( 200 );             // Number of intervals in the eta direction
  std::size_t MB( M * 100 );        // Number of eta intervals in the base flow ODE
  double beta( 0.5 );               // Hartree parameter
  double zeta0( 1.0 );              // Transpiration width
  double K( 10.75 );                  // Transpiration parameter ( +ve = blowing )
  double alpha( 0.8 );              // Wavenumber (alpha hat)
  double Rx( 5000 * 5000 );         // Local Reynolds number
  double Sigma( 1.0 );              // Wave amplitude
  double tol( 1e-4 );               // Tolerance for c_i = 0
  double conv_tol( 1e-8 );          // Convergence tolerance for inner loop

  double K_min( 0.0 );
  double K_step( 0.25 );
  double Sigma_step_initial( 0.5 );
  double Sigma_step( Sigma_step_initial );

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

  cout << "*** Solving the self-similar equations ( without forcing ) ***" << endl;
  cout << "  * N = " << N << ", M = " << M << endl;
  Timer timer;
  timer.start();

  // Solve self-similar equations
  SSI.set_output_path();
  SSI.mesh_setup();
  //SSI.solve_check_exists();
  SSI.solve();

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
  //std::complex<double> c_guess;

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
  //c_guess = eigenvalues[0];

  // Set initial guess for next time
  //orrsommerfeld_2D.set_initial_guess_from_evec( 0 );
  target = eigenvalues[0];
  orrsommerfeld_2D.set_target( target );
  double c_i( eigenvalues[0].imag() ); // Imaginary part of eigenvalue

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


  cout << "*** Solving VWI equations ***" << endl;

  // Create the VWI object
  myVWI vwi( SSI, alpha, Rx, Sigma );

  // Set the initial guess (perturbation variables Phi, Psi, U, Theta)
  TwoD_node_mesh<double> guess( X_NODES, Y_NODES, 4 );
  for ( std::size_t i = 0; i < X_NODES.size(); ++i )
  {
    for ( std::size_t j = 0; j < Y_NODES.size(); ++j )
    {
      guess( i, j, 0 )  = sol( i, j, 0 );
      guess( i, j, 1 )  = sol( i, j, 1 );
      guess( i, j, 2 )  = sol( i, j, 2 );
      guess( i, j, 3 )  = sol( i, j, 3 );
    }
  }

  vwi.set_output_path();
  vwi.make_output_directory();
  vwi.make_eigenvalues_directory();

  //TrackerFile metric( "./DATA/ci_vs_Sigma_K_" + Utility::stringify( SSI.injection(), 3 ) + ".dat" );
  //metric.push_ptr( &vwi.Sigma(), "Sigma" );
  //metric.push_ptr( &c_i, "c_i" );
  //metric.header();

  TrackerFile metric( "./DATA/K_vs_Sigma_alpha_" + Utility::stringify( alpha, 3 ) + ".dat" );
  metric.push_ptr( &vwi.injection(), "K" );
  metric.push_ptr( &vwi.Sigma(), "Sigma" );
  metric.push_ptr( &c_i, "c_i" );
  metric.header();

  do{ // Step in K

    double c_i_minus, c_i_plus, Sigma_minus, Sigma_plus;

    do{ // Step in Sigma

        // Iterate to converge inner loop
        std::size_t iteration( 0 );
        std::size_t max_iterations( 50 );
        double c_i_old( c_i );
        double c_i_diff( 0.0 );

        do{


            cout << "*** Solving the streak equations (with forcing) ***" << endl;
            vwi.solve_streak( guess, eigenvectors );
            // update solution in SSI object
            for ( std::size_t i = 0; i < N + 1; ++i )
            {
              double hzeta=HZETA_NODES[i];
              for ( std::size_t j = 0; j < M + 1; ++j )
              {
                double eta=ETA_NODES[j];
                sol( i, j, 0 )  = guess( i, j, 0 );
                sol( i, j, 1 )  = guess( i, j, 1 );
                sol( i, j, 2 )  = guess( i, j, 2 );
                sol( i, j, 3 )  = guess( i, j, 3 );
                sol( i, j, 4 )  = guess( i, j, 0 )
                                  + base.get_interpolated_vars( eta )[PhiB];
                sol( i, j, 5 )  = guess( i, j, 1 )
                                  + hzeta * base.get_interpolated_vars( eta )[PsiB];
                sol( i, j, 6 )  = guess( i, j, 2 )
                                  + base.get_interpolated_vars( eta )[UB];
                sol( i, j, 7 )  = guess( i, j, 3 )
                                  + hzeta * base.get_interpolated_vars( eta )[ThetaB];
              }
            }
            SSI.set_solution( sol );

            // Solve the OS equations with updated streak
            cout << "*** Solving the stability equations for v and w ***" << endl;
            orrsommerfeld_2D.update_SSI( SSI );
            timer_OS.start();
            orrsommerfeld_2D.solve_evp();
            orrsommerfeld_2D.output();
            timer_OS.print();
            timer_OS.stop();

            eigenvalues = orrsommerfeld_2D.eigenvalues();
            eigenvectors = orrsommerfeld_2D.eigenvectors(); // v, w, q, s
            // Normalise the eigenvectors
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
            target = eigenvalues[0];
            orrsommerfeld_2D.set_target( target );

            c_i = eigenvalues[0].imag();
            c_i_diff = c_i - c_i_old;
            c_i_old = c_i;

            cout << "  * iter = " << iteration << endl;
            cout << "  * c_i_diff = " << c_i_diff << endl;
            ++iteration;

        }while( ( std::abs( c_i_diff ) > conv_tol ) && ( iteration < max_iterations ) );

        if ( iteration > 10 ){ Sigma_step *= 0.9; }
        if ( iteration < 5  ){ Sigma_step *= 1.1; }

        // Decide how to vary Sigma
        if ( c_i > 0.0 && std::abs( c_i ) > tol )
        {
          cout << "  * c_i = " << c_i << endl;
          cout << "  * Sigma = " << vwi.Sigma() << endl;
          c_i_plus = c_i;
          Sigma_plus = vwi.Sigma();
          vwi.Sigma() = Sigma_minus
  			     + ( c_i_minus * ( Sigma_minus - Sigma_plus ) / ( c_i_plus - c_i_minus ) );
          //Sigma_step *= 0.1;
          cout << "*** Stepping in sigma ***" << endl;
          //cout << "  * Step = " << Sigma_step << endl;
        }
        else if ( c_i < 0.0 && std::abs( c_i ) > tol )
        {
          cout << "  * c_i = " << c_i << endl;
          cout << "  * Sigma = " << vwi.Sigma() << endl;
          c_i_minus = c_i;
          Sigma_minus = vwi.Sigma();
          vwi.Sigma() += Sigma_step;
          cout << "*** Stepping in sigma ***" << endl;
          cout << "  * Step = " << Sigma_step << endl;
        }

      }while( std::abs( c_i ) > tol );

      cout << "  * ------------------------------------------ * " << endl;
      cout << "  * K = " << vwi.injection() << endl;
      cout << "  * Rx^1/2 = " << sqrt( vwi.Reynolds() )
           << "  =>  Rx = " << vwi.Reynolds() << endl;
      cout << "  * Sigma = " << vwi.Sigma() << endl;
      cout << "  * c = " << eigenvalues[0] << endl;
      cout << "  * ------------------------------------------ * " << endl;

      metric.update();
      sol.dump_gnu("./DATA/VWI/Streak_K_" + Utility::stringify( vwi.injection(), 5 ) + ".dat");
      eigenvectors.dump_gnu("./DATA/VWI/Evecs_K_" + Utility::stringify( vwi.injection(), 5 ) + ".dat");
      vwi.injection() -= K_step;
      cout << "*** Stepping in K ***" << endl;
      //Sigma_step = Sigma_step_initial; // Reset sigma step

  }while( vwi.injection() >= 0.0 );

	cout << "FINISHED" << endl;
}
