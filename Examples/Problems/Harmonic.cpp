/// Solving the Harmonic equation //TODO - description

#include <cassert>

#include "Core"
#include "ODE_EVP.h"

// enumerate the variables in the ODE
enum {f, fd };

namespace TSL
{
  namespace Problem
  {
    /// Define the harmonic equation by inheriting the Equation base class
    class harmonic_equation : public Equation_2matrix<double>
    {
    public:
      /// The harmonic equation is a 2nd order real ODE
      harmonic_equation() : Equation_2matrix<double>( 2 ) {}

      /// The harmonic equation
      void residual_fn( const Vector<double> &z, Vector<double> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = 0.0;
      }

      /// matrix to multiply the BVP coordinate
      void matrix0( const Vector<double>& z, Matrix<double>& m ) const
      {
        m.fill_diag( 1.0 );
      }

      /// Define the eigenvalue terms by providing the mass matrix
      /// This defines the term lambda * z[ f ] ;
      void matrix1( const Vector<double>& z, Matrix<double>& m ) const
      {
        // eigenvalue multiplies unknown 0 in equation 1
        m( 1, 0 ) = 1.0;
      }

    };

    class harmonic_both_BC : public Residual<double>
    {
    public:
      // 1 reisudal and 2 unknowns
      harmonic_both_BC() : Residual<double> ( 1, 2 ) {}

      void residual_fn( const Vector<double> &z, Vector<double> &B ) const
      {
        B[ 0 ] = z[ f ];
      }
    };

  } // End of namespace Problem
} // End of namespace TSL

using namespace TSL;
using namespace std;

int main()
{
  cout << "Harmonic equation eigenvalue problem (EVP)" << endl;

  Problem::harmonic_equation problem;
  Problem::harmonic_both_BC BC_both;
  // Define the domain
  double left( 0.0 );
  double right( 1.0 );
  // Number of nodes
  std::size_t N( 256 );
  // EV is pi^2 so we'll guess 10
  double guess( 10.0 );
  // Pass/Fail tolerance
  //const double tol( 1e-3 );

  // Setup the problem
  Vector<double> nodes;
  nodes.linspace( left, right, N );
  ODE_EVP<double> ode_evp( &problem, nodes, &BC_both, &BC_both );

  try
  {
    ode_evp.eigensolve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  Vector< std::complex<double> > evals;
  evals = ode_evp.eigenvalues();

  for (size_t i=0; i < evals.size(); ++i)
  {
      if ( evals[i].real() < guess + 1. && evals[i].real() > guess - 1. )
      {
          cout << evals[i] << endl;
      }
  }

  cout << "The exact answer is pi^2 = " << M_PI * M_PI << endl;

  //cout << "eigenvalues = " << ode_evp.eigenvalues() << endl;

  cout << "FINISHED" << endl;
}
