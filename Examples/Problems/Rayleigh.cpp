/// Solving the Rayleigh equation //TODO - description

#include <cassert>

#include "Core"
#include "Eigenvalue"
#include "Rayleigh.h"

// enumerate the variables in the ODE
enum {phi, phid };

namespace TSL
{
  namespace Problem
  {
    // Base flow in the complex plane
    OneD_node_mesh<std::complex<double>, std::complex<double>> baseflow;
    // Rayleigh wavenumber
    double alpha;
    // Base flow profile
    std::complex<double> U( const std::complex<double>& y )
    {
      return std::sin( y );
    }
    // Curvature of the base flow
    std::complex<double> Udd( const std::complex<double>& y )
    {
      return - std::sin( y );
    }

  } // End of namespace Problem
} // End of namespace TSL

using namespace TSL;
using namespace std;

int main()
{
  cout << "=== Rayleigh equation eigenvalue problem (EVP) ===" << endl;

  // Define the problem
  Problem::alpha = 0.8;               // wavenumber
  double tol( 1.e-5 );                // tolerance
  double left( 0.0 );                 // y1 = 0
  double right(  2 * M_PI );          // y2 = 2*pi
  unsigned N( 201 );

  // Real distribution of nodes
  Vector<double> r_nodes;
  r_nodes.linspace( left, right, N );
  Vector<std::complex<double>> c_nodes( r_nodes.size(), 0.0 );
  // Distribution of nodes in the complex plane
  for ( unsigned i = 0; i < N; ++i )
  {
    c_nodes[i] = r_nodes[i];
    std::complex<double> y( r_nodes[ i ] );
    c_nodes[ i ] -= .2 * std::complex<double>( 0.0, 1.0 ) * y * std::exp( - y );
  }
  //cout << "r_nodes = " << r_nodes << endl;
  //cout << "c_nodes = " << c_nodes << endl;

  // Make a base flow on the complex distribution of nodes
  OneD_node_mesh<std::complex<double>, std::complex<double>> base( c_nodes, 2 );
  for ( unsigned i = 0; i < c_nodes.size(); ++i )
  {
    std::complex<double> y = c_nodes[ i ];
    base( i, 0 ) = Problem::U( y );
    base( i, 1 ) = Problem::Udd( y );
  }

  // Make the Rayleigh eigenvalue problem
  TSL::Rayleigh<std::complex<double>> rayleigh( base, Problem::alpha );
  // Solve the global eigenvalue problem
  rayleigh.global_evp();

  // Output the eigenvalues
  cout << rayleigh.eigenvalues() << endl;


  cout << "FINISHED" << endl;
}
