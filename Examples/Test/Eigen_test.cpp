// Test the Eigensystem class
#include "Eigenvalue"

using namespace std;

int main()
{
  cout << "----- TESTING Eigensystem -----" << endl;

  /* ----- TESTING Eigensystem class ----- */	

  TSL::Eigensystem<double> Eig;           // Test constructor

  // Test eigenvectors_computed method
  if ( Eig.eigenvectors_computed() ) { cout << "Computed" << endl; }
  else { cout << "Not computed" << endl; }

  double N = 4.0;
  double d2 = (1.0 / (N - 1.0)) * (1.0 / (N - 1.0));

  TSL::Matrix<double> A_mat(4,4,0.0);
  TSL::Matrix<double> B_mat(4,4,0.0);

  A_mat(0,0) = 1.0;  
  A_mat(1,0) = 1.0/d2; A_mat(1,1) = -2.0/d2; A_mat(1,2) = 1.0/d2; A_mat(1,3) = 0.0; 
  A_mat(2,1) = 1.0/d2; A_mat(2,2) = -2.0/d2; A_mat(2,3) = 1.0/d2;
  A_mat(3,3) = 1.0;

  B_mat(0,0) = -1.0; B_mat(1,1) = -1.0; B_mat(2,2) = -1.0; B_mat(3,3) = -1.0;

  cout << "A = " << endl << A_mat << endl;
  cout << "B = " << endl << B_mat << endl;

  // Test compute method
  bool compute_eigenvectors = true;
  Eig.compute( A_mat, B_mat, compute_eigenvectors );

  cout << "Eigenvalues = " << endl << Eig.eigenvalues() << endl;
  cout << "Alphas = " << endl << Eig.alphas() << endl;
  cout << "Betas = " << endl << Eig.betas() << endl;
  cout << "Eigenvectors = " << endl << Eig.eigenvector_matrix() << endl;
  std::vector< TSL::Vector< std::complex<double> > > evecs;
  evecs = Eig.eigenvectors();
  for (std::size_t i=0; i<evecs.size(); ++i)
  {
    cout << "Eigenvector[" << i << "] = " << endl << evecs[i] << endl;
  }

	cout << "FINISHED" << endl;

}
