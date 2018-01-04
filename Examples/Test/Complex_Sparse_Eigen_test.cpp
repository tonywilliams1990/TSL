/// Solves a 4x4 complex generalised eigenvalue problem
/// \f[ A_{4x4} \,{\underline x}_i = \lambda_i\, B_{4x4}\, {\underline x}_i \f]
/// for the 4 eigenvalues \f$ \lambda_i \f$, \f$i=1,2,3,4.\f$.
/// As a test case we use the SLEPc library.
/// In this case \f$ A_{4x4} \f$ and \f$ B_{4x4} \f$ are such that
/// the eigenvalues are \f$ 3-9i,\, 2-5i,\, 3-i,\, 4-5i \f$. The computation
/// requests eigenvalues that satisfy \f$ \vert\lambda\vert < 4\f$,
/// which in this case is \f$ 3-i \f$.

#include <cassert>
#include <iostream>

#include "Core"
#include "Eigenvalue"

using namespace TSL;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: complex generalised eigenvalue problem  ====\n";
  cout << "\n";

#ifndef PETSC_Z

  cout << " PETSC complex support has not been included\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";

#else

#ifndef SLEPC

  cout << " SCLEPC/PETSC support has not been included\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";

#else

  SlepcInitialize(NULL,NULL,(char*)0,(char*)0);
  SparseMatrix< std::complex<double> > a( 4, 4 );

  a( 0, 0 ) = std::complex<double>( -21.10, -22.50 );
  a( 0, 1 ) = std::complex<double>( 53.50, -50.50 );
  a( 0, 2 ) = std::complex<double>( -34.50, 127.50 );
  a( 0, 3 ) = std::complex<double>( 7.50, 0.50 );

  a( 1, 0 ) = std::complex<double>( -0.46, -7.78 );
  a( 1, 1 ) = std::complex<double>( -3.50, -37.50 );
  a( 1, 2 ) = std::complex<double>( -15.50, 58.50 );
  a( 1, 3 ) = std::complex<double>( -10.50, -1.50 );

  a( 2, 0 ) = std::complex<double>( 4.30, -5.50 );
  a( 2, 1 ) = std::complex<double>( 39.70, -17.10 );
  a( 2, 2 ) = std::complex<double>( -68.50, 12.50 );
  a( 2, 3 ) = std::complex<double>( -7.50, -3.50 );

  a( 3, 0 ) = std::complex<double>( 5.50, 4.40 );
  a( 3, 1 ) = std::complex<double>( 14.40, 43.30 );
  a( 3, 2 ) = std::complex<double>( -32.50, -46.00 );
  a( 3, 3 ) = std::complex<double>( -19.00, -32.50 );

  SparseMatrix< std::complex<double> > b( 4, 4 );

  b( 0, 0 ) = std::complex<double>( 1.00, -5.00 );
  b( 0, 1 ) = std::complex<double>( 1.60, 1.20 );
  b( 0, 2 ) = std::complex<double>( -3.00, 0.00 );
  b( 0, 3 ) = std::complex<double>( 0.00, -1.00 );

  b( 1, 0 ) = std::complex<double>( 0.80, -0.60 );
  b( 1, 1 ) = std::complex<double>( 3.00, -5.00 );
  b( 1, 2 ) = std::complex<double>( -4.00, 3.00 );
  b( 1, 3 ) = std::complex<double>( -2.40, -3.20 );

  b( 2, 0 ) = std::complex<double>( 1.00, 0.00 );
  b( 2, 1 ) = std::complex<double>( 2.40, 1.80 );
  b( 2, 2 ) = std::complex<double>( -4.00, -5.00 );
  b( 2, 3 ) = std::complex<double>( 0.00, -3.00 );

  b( 3, 0 ) = std::complex<double>( 0.00, 1.00 );
  b( 3, 1 ) = std::complex<double>( -1.80, 2.40 );
  b( 3, 2 ) = std::complex<double>( 0.00, -4.00 );
  b( 3, 3 ) = std::complex<double>( 4.00, -5.00 );

  // a vector for the eigenvalues
  Vector< std::complex<double> > lambdas;
  // eigenvalues are: (3,-9), (2,-5), (3,-1), (4,-5)

  SparseEigenSystem< std::complex<double> > system( &a, &b );
  try
  {
    system.eigensolve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  // tag any eigenvalues within a distance of 0.1 of the point 3-i
  //system.set_shift( std::complex<double>( 3.0, -1.0 ) );
  //system.tag_eigenvalues_disc( + 1, 0.1 );
  // get those tagged eigenvalues
  //lambdas = system.get_tagged_eigenvalues();
  lambdas = system.eigenvalues();
  const double tol = 1.e-10;
  if ( std::abs( lambdas[ 0 ] - std::complex<double>( 3.0, -1.0 ) ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout.precision( 12 );
    cout << "    Final error = " << std::abs( lambdas[ 0 ] - std::complex<double>( 3.0, -1.0 ) ) << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

  SlepcFinalize();

  cout << "FINISHED" << endl;

#endif
#endif
}
