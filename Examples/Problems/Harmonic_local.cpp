/// Solving the Harmonic equation f''(x) + \lambda f(x) = 0 as a local eigenvalue
/// problem with homogeneous boundary conditions. The nonlinear BVP solver is
/// used to refine a guess at an eigenvalue
#include <cassert>

#include "Core"
#include "Eigenvalue"

// enumerate the variables in the ODE
enum {f, fd, lambda };

namespace TSL
{
  namespace Problem
  {
    /// Define the harmonic equation by inheriting the Equation base class
    class harmonic_equation : public Equation<std::complex<double>, double>
    {
    public:
      /// The harmonic equation is a 2nd order real ODE
      harmonic_equation() : Equation<std::complex<double>, double>( 3 ) {}

      /// The harmonic equation
      void residual_fn( const Vector<std::complex<double>> &z, Vector<std::complex<double>> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = -z[ lambda ] * z[ f ];
        g[ lambda ] = 0.0;
      }

    };

    class harmonic_left_BC : public Residual<std::complex<double>>
    {
    public:
      // 1 reisudal and 2 unknowns
      harmonic_left_BC() : Residual<std::complex<double>> ( 2, 3 ) {}

      void residual_fn( const Vector<std::complex<double>> &z, Vector<std::complex<double>> &B ) const
      {
        B[ 0 ] = z[ f ];
        B[ 1 ] = z[ fd ] - 1.0; // arbitrary amplitude
      }
    };

    class harmonic_right_BC : public Residual<std::complex<double>>
    {
    public:
      // 1 reisudal and 2 unknowns
      harmonic_right_BC() : Residual<std::complex<double>> ( 1, 3 ) {}

      void residual_fn( const Vector<std::complex<double>> &z, Vector<std::complex<double>> &B ) const
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
  cout << "*** Harmonic equation local eigenvalue problem (EVP) ***" << endl;

  Problem::harmonic_equation problem;
  Problem::harmonic_left_BC BC_left;
  Problem::harmonic_right_BC BC_right;
  // Define the domain
  double left( 0.0 );
  double right( 1.0 );
  // Guess for the eigenvalue
  std::complex<double> eigenvalue_guess = 9.0;
  // Number of nodes
  cout << "  * The exact answer is pi^2 = " << M_PI * M_PI << endl; // n^2 * pi^2
  cout << "  * -------------------------------------------" << endl;
  for ( std::size_t i = 6; i <=12; ++i )
  {
    // Setup the problem
    std::size_t N( std::pow( 2, i ) );
    Vector<double> nodes;
    nodes.linspace( left, right, N );
    ODE_BVP< std::complex<double> > ode( &problem, nodes, &BC_left, &BC_right );
    bool compute_eigenvectors( false );
    // Set an initial guess
    for ( unsigned i = 0; i < N; ++i )
    {
      double x = ode.solution().coord( i );
      ode.solution()( i, f ) = x * ( 1 - x );
      ode.solution()( i, fd ) = 1.0 - 2 * x;
      ode.solution()( i, lambda ) = eigenvalue_guess;
    }

    try
    {
      ode.solve_bvp();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }


    double abs_error( std::abs( ode.solution()( 1, lambda ) - M_PI * M_PI ) );
    cout << "  * N = " << N << endl;
    cout << "  * lambda = " << ode.solution()( 1, lambda ) << endl;
    cout << "  * |error| = " << abs_error << std::endl;
    cout << "  * -------------------------------------------" << endl;
  }

  cout << "FINISHED" << endl;
}
