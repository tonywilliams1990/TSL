// Test the Vector class
#include "Core"

using namespace std;

int main()
{
  cout << "----- TESTING linear system solver -----" << endl;

  TSL::Vector<double> b(3,1.0);
  cout << "b = " << endl << b << endl;

  TSL::Matrix<double> A(3,3,1.0);
  A(0,1) = 2.0; A(0,2) = 3.0;
  A(1,0) = 4.0; A(1,1) = 5.0;
  A(2,1) = 8.0; A(2,2) = 9.0;
  cout << "A = " << endl << A << endl;

  TSL::Vector<double> x(3,0.0);
  TSL::Timer timer;                                        // Timer
  cout << "Solving using Eigen and LU factorisation" << endl;
  timer.start();
  x = A.solve( b, "LU" );
  cout << "x = " << endl << x << endl;
  timer.print();
  timer.stop();

  TSL::Matrix<double> B(3,2,1.0);
  cout << "B = " << endl << B << endl;
  TSL::Matrix<double> X(3,2,0.0);
  X = A.solve( B, "QR" );

  cout << "X = " << endl << X << endl;

	cout << "FINISHED" << endl;

}
