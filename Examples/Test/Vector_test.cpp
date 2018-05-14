// Test the Vector class
#include "Core"

using namespace std;

int main()
{
  cout << "----- TESTING Vector -----" << endl;

  /* ----- TESTING Vector class ----- */

  TSL::Vector<double> v1( 10, 1.0 );            // Create a vector object
  v1[2] = 0.4; v1[7] = 0.3;                     // Test indexing
  cout << "v1 = " << endl << v1 << endl;        // Test output operator

  TSL::Vector<double> v2;
  v2 = v1;                                      // Test assignment
  cout << "v2 = " << endl << +v2 << endl;       // Test unary +
  cout << "-v2 = " << endl << -v2 << endl;      // Test unary -
  TSL::Vector<double> v3;
  v3 = v1 + v2;                                 // Test binary +
  cout << "v3 = " << endl << v3 << endl;
  TSL::Vector<double> v4;
  v4 = v1 - v2;                                 // Test binary -
  cout << "v4 = " << endl << v4 << endl;

  v3 = v1 * 2.5;                                // Test scalar multiplication
  cout << "v3 = " << endl << v3 << endl;

  v4 = 2.5 * v2;                                // On either side
  cout << "v4 = " << endl << v4 << endl;

  v3 = v1 / 2.0;                                // Test scalar division
  cout << "v3 = " << endl << v3 << endl;

  v3 += v1;                                     // Test addition assignment
  cout << "v3 = " << endl << v3 << endl;

  v3 -= v1;                                     // Test subtraction assignment
  cout << "v3 = " << endl << v3 << endl;

  cout << "v3.size() = " << v3.size() << endl;  // Test size method
  v3.push_back( 3.14 );                         // Test push_back method
  cout << "v3 = " << endl << v3 << endl;
  cout << "v3.size() = " << v3.size() << endl;

  v3.resize( 4 );
  cout << "v3 = " << endl << v3 << endl;
  cout << "v3.size() = " << v3.size() << endl;

  v3.clear();
  cout << "v3 = " << endl << v3 << endl;
  cout << "v3.size() = " << v3.size() << endl;

  v4[ 0 ] = -1.1;
  v3 = v4.abs();                                // Test abs method
  cout << "v3 = " << endl << v3 << endl;

  v3.swap(0,1);                                 // Test swap method
  cout << "v3 = " << endl << v3 << endl;

  v3.assign( 3, 1.0 );                          // Test assign method
  cout << "v3 = " << endl << v3 << endl;

  v3.linspace( 0.0, 1.0, 11 );                  // Test linspace method
  cout << "v3 = " << endl << v3 << endl;

  cout << "v3.sum() = " << v3.sum() << endl;          // Test sum method
  cout << "v3.product() = " << v3.product() << endl;  // Test product method

  v3.scale( 10 );                                     // Test scale method
  cout << "v3 = " << endl << v3 << endl;

  cout << "v1.dot(v2) = " << v1.dot(v2) << endl;      // Test dot method

  cout << "v3.norm_1() = " << v3.norm_1() << endl;    // Test norm_1 method
  cout << "v3.norm_2() = " << v3.norm_2() << endl;    // Test norm_2 method
  cout << "v3.norm_p(3) = " << v3.norm_p(3) << endl;  // Test norm_p method
  cout << "v3.norm_inf() = " << v3.norm_inf() << endl;// Test norm_inf method

  cout << "v3.conjugate() = " << endl << v3.conjugate() << endl;

  std::complex<double> elem(0.3, 0.5);
  TSL::Vector< std::complex<double> > v_cmplx( 4, elem );
  cout << "v_cmplx = " << endl << v_cmplx << endl;
  cout << "v_cmplx.conjugate() = " << endl << v_cmplx.conjugate() << endl;

  // Test TwoD_node_mesh
  TSL::Vector<double> x_nodes;
  x_nodes.linspace( 0.0, 1.0, 11 );
  TSL::Vector<double> y_nodes( x_nodes );
  TSL::TwoD_node_mesh< std::complex<double> > mesh_2d( x_nodes, y_nodes, 1 );
  for ( std::size_t i=0; i<x_nodes.size(); ++i )
  {
    for ( std::size_t j=0; j<y_nodes.size(); ++j )
    {
      mesh_2d( i, j, 0 ) = std::complex<double>( 0.1*i, 0.2*j );
    }
  }

  //mesh_2d.dump();

  //mesh_2d.conjugate().dump();


	cout << "FINISHED" << endl;

}
