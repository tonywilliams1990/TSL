// Test the SelfSimInjection class
#include "Core"
#include "SelfSimInjection.h"

using namespace std;


int main()
{
  cout << "----- TESTING SelfSimInjection -----" << endl;

  // Test constructor
  TSL::SelfSimInjection SSI;

  // Setup
  SSI.hzeta_right() = 16;
  SSI.eta_top() = 128;
  SSI.injection() = 2;
  SSI.set_mesh( "NONUNIFORM" );
  SSI.set_base_flow( "2D" );
  SSI.speed_up( false );

  // Output setup
  cout << "hzeta_right = " << SSI.hzeta_right() << endl;
  cout << "eta_top = " << SSI.eta_top() << endl;
  cout << "N = " << SSI.hzeta_intervals() << endl;
  cout << "M = " << SSI.eta_intervals() << endl;
  cout << "beta = " << SSI.hartree() << endl;
  cout << "KB = " << SSI.base_injection() << endl;
  cout << "zeta0 = " << SSI.injection_width() << endl;
  cout << "K = " << SSI.injection() << endl;
  cout << "The mesh is " << SSI.mesh() << " and the base flow is "
       << SSI.base_flow() << endl;

  // Check mesh functions
  cout << "X( 0.1 ) = " << SSI.mesh_X( 0.1 ) << endl;
  cout << "Xd( 0.1 ) = " << SSI.mesh_Xd( 0.1 ) << endl;
  cout << "Xdd( 0.1 ) = " << SSI.mesh_Xdd( 0.1 ) << endl;
  cout << "Y( 0.1 ) = " << SSI.mesh_Y( 0.1 ) << endl;
  cout << "Yd( 0.1 ) = " << SSI.mesh_Yd( 0.1 ) << endl;
  cout << "Ydd( 0.1 ) = " << SSI.mesh_Ydd( 0.1 ) << endl;

  SSI.solve();
  TSL::TwoD_node_mesh<double> sol = SSI.solution();
  //SSI.output();
  //SSI.iterate_on_zeta0( 1.0, 4.0 );
  cout << "zeta0 = " << SSI.injection_width() << ", A = " << SSI.mass_flux() << endl;
  cout << "U_eta(0,0) = " << SSI.shear_at_origin() << endl;
  cout << "eta_half = " << SSI.eta_half() << endl;

	cout << "FINISHED" << endl;

}
