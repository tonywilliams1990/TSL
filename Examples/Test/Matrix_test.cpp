// Test the Matrix class
#include "Core"

using namespace std;

int main()
{
  cout << "----- TESTING Matrix -----" << endl;

  /* ----- TESTING Tmatrix class ----- */	

	TSL::Matrix<double> A;                              // Test basic constructor
  TSL::Matrix<double> B( A );                         // Test copy constructor
	TSL::Matrix<double> C(3,3,0.5);                     // Test constructor with elements
  
  C(0,0) = 2.22; C(2,2) = 1.2;                        // Test indexing operator

  cout << "C = " << endl << C << endl;                // Test output operator

  A = C;                                              // Test assignment 
  cout << "A = " << endl << A << endl;
  A(2,0) = 0.7;   
  cout << "+A = " << endl << +A << endl;              // Test unary +
  cout << "-A = " << endl << -A << endl;              // Test unary -
  cout << "A + C = " << endl << A + C << endl;        // Test binary +  
  cout << "A - C = " << endl << A - C << endl;        // Test binary -

  cout << "A * 3.0 = " << endl << A * 3.0 << endl;    // Test scalar multiplication 
  cout << "2.0 * A = " << endl << 2.0 * A << endl;    // Test scalar multiplication
  cout << "A / 3.0 = " << endl << A / 3.0 << endl;    // Test scalar division
  
  cout << "A * C = " << endl << A * C << endl;        // Test matrix multiplication 
  
  A += C;                                             // Test addition assignment +=
  cout << "A += C: A = " << endl << A << endl;

  A -= C;                                             // Test subtraction assignment -=
  cout << "A -= C: A = " << endl << A << endl;

  cout << "A.rows() = " << A.rows() << endl;          // Test rows method
  cout << "A.cols() = " << A.cols() << endl;          // Test cols method
  cout << "A.numel() = " << A.numel() << endl;        // Test numel method

  cout << "A.transpose() = " << endl << A.transpose() << endl;  // Test transpose method
  cout << "A.adjoint() = " << endl << A.adjoint() << endl;      // Test adjoint method
  cout << "A.conjugate() = " << endl << A.conjugate() << endl;  // Test conjugate method

  A.resize(2,2);                                      // Test resize method
  cout << "A = " << endl << A << endl;

  A(0,0) = 0.3; A(1,0) = 2.1;
  cout << "A = " << endl << A << endl;
  A.transpose_in_place();                             // Test transpose in place
  cout << "A = A^T: A = " << endl << A << endl;

  A.fill( 1.1 );                                      // Test fill method
  cout << "A = " << endl << A << endl;

  C.fill_diag( 1.3 );                                 // Test fill_diag method
  cout << "C = " << endl << C << endl;

  C.fill_band( 1, 2.7 );                              // Test fill_band method
  cout << "C = " << endl << C << endl;

  C.fill_tridiag( 0.5, 0.2, 0.3 );                    // Test fill_tridiag method

  C(0,0) = 4.0; C(0,1) = 1.0; C(0,2) = 2.0;
  C(1,0) = 2.0; C(1,1) = -1.; C(1,2) = 3.0;
  C(2,0) = 1.0; C(2,1) = 2.0; C(2,2) = 7.0;
  cout << "C = " << endl << C << endl;

  cout << "C.norm_1() = " << C.norm_1() << endl;      // Test norm_1 method
  cout << "C.norm_inf() = " << C.norm_inf() << endl;  // Test norm_inf method
  cout << "C.norm_frob() = " << C.norm_frob() << endl;// Test norm_frob method
  cout << "C.norm_max() = " << C.norm_max() << endl;  // Test norm_max method
  cout << "C.norm_p(3.0) = " << C.norm_p(3.0) << endl;// Test norm_p method

  cout << "C.determinant() = " << C.determinant() << endl; // Test determinant method
  
	cout << "FINISHED" << endl;

}
