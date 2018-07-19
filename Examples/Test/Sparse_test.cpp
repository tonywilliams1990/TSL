// Test the sparse matrix class
#include "Core"
#include "Sparse"

using namespace std;

int main()
{
  cout << "----- TESTING SparseMatrix -----" << endl;

  /* ----- TESTING sparse linear system solver ----- */

  size_t N = 1000;

  TSL::SparseMatrix<double> A_matrix(N,N);
  for (size_t i=0; i<N; ++i)
  {
     A_matrix(i,i) = 0.435 * (i + 1.0);
  }

  cout << "A_matrix.rows() = " << A_matrix.rows() << endl;
  cout << "A_matrix.cols() = " << A_matrix.cols() << endl;
  cout << "A_matrix.size() = " << A_matrix.size() << endl;
  cout << "A_matrix.numel() = " << A_matrix.numel() << endl;

  //A_matrix(33, 100) = 2.0;

  cout << "A_matrix.numel_row( 33 ) = " << A_matrix.numel_row( 33 ) << endl;

  // Test copy constructor
  TSL::SparseMatrix<double> A_copy( A_matrix );
  cout << "A_copy.numel() = " << A_copy.numel() << endl;

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

  cout << "A_matrix( 10, 10 ) = " << A_matrix( 10, 10 ) << endl;
  A_matrix.clear();
  cout << "A_matrix( 10, 10 ) = " << A_matrix( 10, 10 ) << endl;
  A_matrix.eye();
  cout << "A_matrix( 10, 10 ) = " << A_matrix( 10, 10 ) << endl;
  A_matrix.scale( 2.0 );
  cout << "A_matrix( 10, 10 ) = " << A_matrix( 10, 10 ) << endl;
  //A_matrix.print();

  TSL::SparseMatrix<double> B_matrix(N,N);
  B_matrix = A_matrix;
  B_matrix = A_matrix * 2.0;
  B_matrix = 2.0 * A_matrix;
  //B_matrix.print();

  TSL::SparseMatrix<double> C_matrix(3,3);
  C_matrix(0,0) = 1.0; C_matrix(0,1) = -2.0; //C_matrix(0,2) =  3.0;
  C_matrix(1,0) = 5.0; C_matrix(1,1) =  8.0; C_matrix(1,2) = -1.0;
  //C_matrix(2,0) = 2.0;
  C_matrix(2,1) =  1.0; C_matrix(2,2) =  1.0;
  C_matrix.print();
  C_matrix.output( "./C_matrix" );

  //B_matrix = C_matrix;
  //cout << "B_matrix = ";
  //B_matrix.print();

  TSL::SparseMatrix<double> D_matrix(N,N);
  Eigen::SparseMatrix<double, Eigen::ColMajor, long long> C_matrix_Eigen(3,3);
  C_matrix_Eigen = C_matrix.convert_to_Eigen();

  for (int k = 0; k < C_matrix_Eigen.outerSize(); ++k){
    for (Eigen::SparseMatrix<double, Eigen::ColMajor, long long>::InnerIterator it(C_matrix_Eigen, k); it; ++it){
        std::cout << it.row() <<"\t";
        std::cout << it.col() << "\t";
        std::cout << it.value() << std::endl;
    }
  }

  Eigen::SparseLU< Eigen::SparseMatrix<double, Eigen::ColMajor, long long> > solver;
  solver.compute( C_matrix_Eigen );

  // Test complex sparse matrices
  TSL::SparseMatrix<std::complex<double>> E_matrix( 2, 2 );
  TSL::Vector<std::complex<double>> E_vector( 2, 0.0 );
  E_matrix( 0, 0 ) = std::complex<double>(1.0,1.0);
  E_matrix( 0, 1 ) = std::complex<double>(-1.0,0.0);
  E_matrix( 1, 0 ) = std::complex<double>(1.0,-1.0);
  E_matrix( 1, 1 ) = std::complex<double>(1.0,1.0);
  E_vector[ 0 ] = std::complex<double>(0.0,1.0);
  E_vector[ 1 ] = std::complex<double>(1.0,0.0);

  //E_matrix.print();
  TSL::Vector<std::complex<double>> E_sol( 2, 0.0 );
  E_sol = E_matrix.solve( E_vector );
  cout << " E_sol = " << endl << E_sol << endl;

	cout << "FINISHED" << endl;

}
