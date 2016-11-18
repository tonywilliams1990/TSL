// Test the sparse matrix class
#include "Sparse"

using namespace std;

int main()
{
  cout << "----- TESTING SparseMatrix -----" << endl;

  TSL::SparseMatrix<double> A(100,100);               // Create a SparseMatrix object

  cout << "A.rows() = " << A.rows() << endl;          // Test rows method
  cout << "A.cols() = " << A.cols() << endl;          // Test cols method
  cout << "A.size() = " << A.size() << endl;          // Test size method
  A.eye();                                            // Test eye method
  A.set_elem( 30, 56, 0.2 );                          // Test set_elem
  cout << "A(30,56) = " << A.get_elem(30,56) << endl; // Test get_elem
  A.value_ref(41,27) = 0.77;                          // Test value_ref method
  cout << "A(41,27) = " << A(41,27) << endl;          // Test indexing (read only)
  A(4,3) = 12.1;                                      // Test indexing (read/write)
  cout << "A(4,3) = " << A(4,3) << endl;
  cout << "A.numel() = " << A.numel() << endl;        // Test numel method

  TSL::SparseMatrix<double> B;
  B = A;                                              // Test assignment
  B = A + B;                                          // Test binary +
  cout << "B(1,1) = " << B(1,1) << endl;
  B = A - B;                                          // Test binary -
  cout << "B(1,1) = " << B(1,1) << endl;
  B.scale( 3.0 );
  cout << "B(1,1) = " << B(1,1) << endl;

  B = B * 2.2;
  cout << "B(1,1) = " << B(1,1) << endl;
  A = 3.3 * A;
  cout << "A(1,1) = " << A(1,1) << endl;

  /* ----- TESTING sparse linear system solver ----- */

  size_t N = 1000000;
  TSL::SparseMatrix<double> A_matrix(N,N);
  for (size_t i=0; i<N; ++i)
  {
     A_matrix(i,i) = 0.435 * (i + 1.0);
  }
  
  TSL::Vector<double> B_vector(N,0.435);
  TSL::Vector<double> X_vector(N,0.0);
  
	X_vector = A_matrix.solve( B_vector );
  //cout << "B_vector = " << endl << B_vector << endl;
  //cout << "X_vector = " << endl << X_vector << endl;
  cout << "X[N-1] = " << X_vector[ N - 1 ] << endl;
  cout << "1 / N  = " << 1.0 / N << endl;

	cout << "FINISHED" << endl;

}
