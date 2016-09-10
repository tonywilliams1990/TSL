// Test the Residual class
#include "Core"

using namespace std;

int main()
{
  cout << "----- TESTING OneD_node_mesh -----" << endl;

  /* ----- TESTING OneD_node_mesh class ----- */	

  TSL::OneD_node_mesh<double,double> default_mesh;        // Test default constructor

  size_t N( 11 );                                         // Number of nodes
  size_t nvars( 2 );                                      // Number of variables

  TSL::Vector<double> uniform;
  uniform.linspace(0,2,N);                                // Uniform vector of nodes

  TSL::OneD_node_mesh<double> mesh( uniform, nvars );     // Create a 1D mesh object
  
  for ( std::size_t i=0; i<N; ++i )
  {
    // Set the first variable equal to 2 * x
    mesh( i, 0 ) = 2 * mesh.coord( i );
    // Set the second variable equal to x^2
    mesh( i, 1 ) = mesh.coord( i ) * mesh.coord( i );
  }

  TSL::Vector<double> vars(2,3.14);
  mesh.set_nodes_vars( 0, vars );                         // Test set_nodes_vars
  cout << "mesh.get_nodes_vars( 0 ) = " << endl
  << mesh.get_nodes_vars( 0 ) << endl;                    // Test get_nodes_vars

  cout << "No nodes = " << mesh.get_nnodes() << endl;     // Test get_nnodes
  cout << "No vars = " << mesh.get_nvars() << endl;       // Test get_nvars
  cout << "nodes = " << endl << mesh.nodes() << endl;     // Test nodes

  cout << "Vars( 1.1 ) = " << endl << mesh.get_interpolated_vars( 1.1 ) << endl;

  mesh.output( "./DATA/mesh_test.dat" , 10 );

	cout << "FINISHED" << endl;

}
