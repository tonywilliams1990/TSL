/// Solving the Rayleigh equation //TODO - description

#include <cassert>

#include "Core"
#include "Eigenvalue"

// enumerate the variables in the ODE
enum {phi, phid };

namespace TSL
{
  unsigned col( const unsigned& i, const unsigned& k )
  {
      // Return the column number for kth variable at node i
      return 2 * i  + k;
  }
  namespace Problem
  {
    double alpha( 0.5 );              // Rayleigh wavenumber

    // Base flow profile U_B
    double U_B( const double& y)
    {
      return std::sin( y );
    }
    // Curvature of the base flows
    double U_Bdd( const double& y)
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

  // Define the domain
  double left( 0.0 );                 // y1 = 0
  double right(  2 * M_PI );          // y2 = 2*pi
  std::size_t N( 128 );               // Number of nodes
  Vector<double> nodes;
  nodes.linspace( left, right, N );

  cout << "*** Assembling the matrices for the eigenvalue problem." << endl;

  // Create the generalised eigenvalue problem A v = lambda B v
  Matrix<double> A( 2*N, 2*N, 0.0 ); // 2N*2N -> 2nd order system
  Matrix<double> B( 2*N, 2*N, 0.0 );

  unsigned row( 0 );                     // Row counter

  // Boundary condition at y1
  unsigned i( 0 );
  // phi(y1) = 0
  A( row, col( i, phi ) )             =  1.0;
  ++row;

  // Interior points

  double alpha2 = Problem::alpha * Problem::alpha;

  for ( std::size_t i=0; i<N-1; ++i)
  {
    double y = nodes[ i ];
    double yplus1 = nodes[ i + 1 ];
    double delta = yplus1 - y;
    double U_B = 0.5 * ( Problem::U_B( y ) + Problem::U_B( yplus1 ) );
    double U_Bdd = 0.5 * ( Problem::U_Bdd( y ) + Problem::U_Bdd( yplus1 ) );


    // Equation 1
    A( row, col( i + 1, phi ) )     =  1.0 / delta;
    A( row, col( i, phi ) )         = -1.0 / delta;
    A( row, col( i + 1, phid ) )    = -0.5;
    A( row, col( i, phid ) )        = -0.5;
    ++row;

    // Equation 2
    A( row, col( i + 1, phid ) )    =  1.0 / delta;
    A( row, col( i, phid ) )        = -1.0 / delta;
    A( row, col( i + 1, phi ) )     = -0.5 * ( alpha2 + ( U_Bdd / U_B ) );
    A( row, col( i, phi ) )         = -0.5 * ( alpha2 + ( U_Bdd / U_B ) );
    B( row, col( i + 1, phid ) )    =  1.0 / ( delta * U_B );
    B( row, col( i, phid ) )        = -1.0 / ( delta * U_B );
    B( row, col( i + 1, phi ) )     = -0.5 * alpha2 / U_B;
    B( row, col( i, phi ) )         = -0.5 * alpha2 / U_B;
    ++row;
  }
  // Boundary condition at y2
  i = N - 1;
  // phi(y2) = 0
  A( row, col( i, phi ) )             =  1.0;
  ++row;

  cout << "*** Solving the generalised eigenvalue problem A v=lambda B v." << endl;

  // Setup the eigensystem and compute the eigenvalues
  Eigensystem<double> system;
  bool compute_eigenvectors = false;
  Timer timer;
  timer.start();
  system.compute( A, B, compute_eigenvectors );

  Vector< std::complex<double> > evals = system.eigenvalues();

  for (size_t i=0; i < evals.size(); ++i)
  {
      if ( evals[i].real() < 1.0 && evals[i].imag() < 1e-3  )
      {
          cout << evals[i] << endl;
      }
  }

  timer.print();                                      // Output time to screen
  timer.stop();                                       // Stop the timer

  cout << "FINISHED" << endl;
}
