// Test the Eigensystem class
#include "Core"
#include "Eigenvalue"

using namespace std;
using namespace TSL;

int main()
{
  cout << "----- TESTING Eigensystem (complex) -----" << endl;

  /* ----- TESTING Eigensystem class ----- */

  // Create eigensystem
  Eigensystem<std::complex<double>> Eig_cmplx;

  // A matrix
  Matrix< std::complex<double> > A( 4, 4, 0.0 );

  A( 0, 0 ) = std::complex<double>( -21.10, -22.50 );
  A( 0, 1 ) = std::complex<double>( 53.50, -50.50 );
  A( 0, 2 ) = std::complex<double>( -34.50, 127.50 );
  A( 0, 3 ) = std::complex<double>( 7.50, 0.50 );

  A( 1, 0 ) = std::complex<double>( -0.46, -7.78 );
  A( 1, 1 ) = std::complex<double>( -3.50, -37.50 );
  A( 1, 2 ) = std::complex<double>( -15.50, 58.50 );
  A( 1, 3 ) = std::complex<double>( -10.50, -1.50 );

  A( 2, 0 ) = std::complex<double>( 4.30, -5.50 );
  A( 2, 1 ) = std::complex<double>( 39.70, -17.10 );
  A( 2, 2 ) = std::complex<double>( -68.50, 12.50 );
  A( 2, 3 ) = std::complex<double>( -7.50, -3.50 );

  A( 3, 0 ) = std::complex<double>( 5.50, 4.40 );
  A( 3, 1 ) = std::complex<double>( 14.40, 43.30 );
  A( 3, 2 ) = std::complex<double>( -32.50, -46.00 );
  A( 3, 3 ) = std::complex<double>( -19.00, -32.50 );

  cout << "A = " << A << endl;

  // B matrix
  Matrix< std::complex<double> > B( 4, 4, 0.0 );

  B( 0, 0 ) = std::complex<double>( 1.00, -5.00 );
  B( 0, 1 ) = std::complex<double>( 1.60, 1.20 );
  B( 0, 2 ) = std::complex<double>( -3.00, 0.00 );
  B( 0, 3 ) = std::complex<double>( 0.00, -1.00 );

  B( 1, 0 ) = std::complex<double>( 0.80, -0.60 );
  B( 1, 1 ) = std::complex<double>( 3.00, -5.00 );
  B( 1, 2 ) = std::complex<double>( -4.00, 3.00 );
  B( 1, 3 ) = std::complex<double>( -2.40, -3.20 );

  B( 2, 0 ) = std::complex<double>( 1.00, 0.00 );
  B( 2, 1 ) = std::complex<double>( 2.40, 1.80 );
  B( 2, 2 ) = std::complex<double>( -4.00, -5.00 );
  B( 2, 3 ) = std::complex<double>( 0.00, -3.00 );

  B( 3, 0 ) = std::complex<double>( 0.00, 1.00 );
  B( 3, 1 ) = std::complex<double>( -1.80, 2.40 );
  B( 3, 2 ) = std::complex<double>( 0.00, -4.00 );
  B( 3, 3 ) = std::complex<double>( 4.00, -5.00 );

  cout << "B = " << B << endl;

  // Solve the complex eigensystem
  bool compute_eigenvectors = true;
  Eig_cmplx.compute( A, B, compute_eigenvectors );

  // Output
  cout << "Eigenvalues = " << endl << Eig_cmplx.eigenvalues() << endl;
  cout << "Alphas = " << endl << Eig_cmplx.alphas() << endl;
  cout << "Betas = " << endl << Eig_cmplx.betas() << endl;
  cout << "Eigenvectors = " << endl << Eig_cmplx.eigenvector_matrix() << endl;

  std::vector< Vector< std::complex<double> > > eigenvectors;
  eigenvectors = Eig_cmplx.eigenvectors();
  cout << "evec[0] = " << endl << eigenvectors[0] << endl;

  // Small test
  /*Matrix< std::complex<double> > A_new( 2, 2, 0.0 );

  A_new( 0, 0 ) = std::complex<double>( 1, 0 );
  A_new( 0, 1 ) = std::complex<double>( 0, 0 );
  A_new( 1, 0 ) = std::complex<double>( 1, 0 );
  A_new( 1, 1 ) = std::complex<double>( 3, 0 );

  Matrix< std::complex<double> > B_new( 2, 2, 0.0 );

  B_new( 0, 0 ) = std::complex<double>( 1, 0 );
  B_new( 0, 1 ) = std::complex<double>( 0, 0 );
  B_new( 1, 0 ) = std::complex<double>( 0, 0 );
  B_new( 1, 1 ) = std::complex<double>( 1, 0 );

  Eigensystem<std::complex<double>> Eig_cmplx_2;
  Eig_cmplx_2.compute( A_new, B_new, true );
  cout << "Eigenvalues = " << endl << Eig_cmplx_2.eigenvalues() << endl;
  cout << "Eigenvectors = " << endl << Eig_cmplx_2.eigenvector_matrix() << endl;*/

	cout << "FINISHED" << endl;

}
