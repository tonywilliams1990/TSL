// Test the sparse matrix class
#include "Core"
#include "Sparse"

using namespace std;

int main()
{
  cout << "----- TESTING SparseMatrix -----" << endl;

  /* ----- TESTING sparse linear system solver ----- */

  size_t N = 1000000;

  TSL::SparseMatrix<double> A_matrix(N,N);
  for (size_t i=0; i<N; ++i)
  {
     A_matrix(i,i) = 0.435 * (i + 1.0);
  }

  cout << "A_matrix.rows() = " << A_matrix.rows() << endl;
  cout << "A_matrix.cols() = " << A_matrix.cols() << endl;
  cout << "A_matrix.size() = " << A_matrix.size() << endl;
  cout << "A_matrix.numel() = " << A_matrix.numel() << endl;
  
  TSL::Vector<double> B_vector(N,0.435);
  TSL::Vector<double> X_vector(N,0.0);

  // Test SparseLU solve method
  TSL::Timer timer;
  timer.start();
	X_vector = A_matrix.solve( B_vector );              // SparseLU solve
  //cout << "B_vector = " << endl << B_vector << endl;
  //cout << "X_vector = " << endl << X_vector << endl;
  cout << "X[N-1] = " << X_vector[ N - 1 ] << endl;
  cout << "1 / N  = " << 1.0 / N << endl;

  timer.print();
  timer.stop();

  A_matrix.clear();
  A_matrix.eye();
  A_matrix.scale( 2.0 );
  //A_matrix.print();

  TSL::SparseMatrix<double> B_matrix(N,N);
  B_matrix = A_matrix;
  B_matrix = A_matrix * 2.0;
  B_matrix = 2.0 * A_matrix;
  //B_matrix.print();

	cout << "FINISHED" << endl;

}
