// Test the Dual class
#include "Core"

using namespace std;

template <class T>
T square( T& x )
{
  return x*x;
}


int main()
{
  cout << "----- TESTING Dual -----" << endl;

  TSL::Dual<double> dual;                                    // Test default constructor
  dual.real() = 5.0;
  TSL::Vector<double> eps(3,1.1);
  eps[1] = 2.22; eps[2] = 3.333;
  dual.epsilons() = eps;
  cout << "dual.real() = " << dual.real() << endl;           // Test real method
  cout << "dual.epsilons() = " << dual.epsilons() << endl;   // Test epsilons method

  TSL::Vector<double> eps_2(3,2.0);
  eps_2[1] = -1.0; eps_2[2] = 1.0;
  TSL::Dual<double> dual_2( 2.0, eps_2 );                    // Test initialised constructor

  cout << "dual_2.real() = " << dual_2.real() << endl;
  cout << "dual_2.epsilons() = " << dual_2.epsilons() << endl;

  TSL::Dual<double> dual_3;
  dual_3 = dual + dual_2;                                     // Test operators
  cout << "dual_3 = " << dual_3 << endl;

  dual_3 = dual_2 - dual;
  cout << "dual_3 = " << dual_3 << endl;

  dual_3 = dual * dual_2;
  cout << "dual_3 = " << dual_3 << endl;

  TSL::Vector<double> eps_deriv(1,1.0);
  TSL::Dual<double> dual_deriv( 3.0, eps_deriv );

  dual_3 = square( dual_deriv );
  cout << "dual_3 = " << dual_3 << endl;

  cout << "dual = " << dual << endl;

  cout << "dual[0] = " << dual[0] << endl;                    // Test indexing
  dual[0] = 3.2;
  cout << "dual[0] = " << dual[0] << endl;

  dual_3 = TSL::cos( dual );
  cout << "dual_3 = " << dual_3 << endl;
  cout << "cos( dual ) = " << TSL::cos( dual ) << endl;

  double p = 3.2;
  cout << "cos(p) = " << TSL::cos( p ) << endl;
  cout << "cos(3.2) = " << TSL::cos( 3.2 ) << endl;

  cout << "sin( dual ) = " << TSL::sin( dual ) << endl;

  cout << "dual_deriv = " << dual_deriv << endl;
  cout << "sin( dual_deriv ) = " << TSL::sin( dual_deriv ) << endl;
  cout << "exp( dual_deriv ) = " << TSL::exp( dual_deriv ) << endl;
  cout << "exp( 3.0 ) = " << TSL::exp( 3.0 ) << endl;
  cout << "log( dual_deriv ) = " << TSL::log( dual_deriv ) << endl;
  cout << "log( 3.0 ) = " << TSL::log( 3.0 ) << endl;
  cout << "tan( dual_deriv ) = " << TSL::tan( dual_deriv ) << endl;
  cout << "tan( 3.0 ) = " << TSL::tan( 3.0 ) << endl;
  cout << "pow( dual_deriv ) = " << TSL::pow( dual_deriv, 2.0 ) << endl;
  cout << "pow( 3.0, 2.0 ) = " << TSL::pow( 3.0, 2.0 ) << endl;
  cout << "sqrt( dual_deriv ) = " << TSL::sqrt( dual_deriv ) << endl;
  cout << "sqrt( 3.0 ) = " << TSL::sqrt( 3.0 ) << endl;
  cout << "abs( dual_deriv ) = " << TSL::abs( dual_deriv ) << endl;
  cout << "abs( 3.0 ) = " << TSL::abs( 3.0 ) << endl;
  cout << "cosh( dual_deriv ) = " << TSL::cosh( dual_deriv ) << endl;
  cout << "cosh( 3.0 ) = " << TSL::cosh( 3.0 ) << endl;
  cout << "sinh( dual_deriv ) = " << TSL::sinh( dual_deriv ) << endl;
  cout << "sinh( 3.0 ) = " << TSL::sinh( 3.0 ) << endl;
  cout << "tanh( dual_deriv ) = " << TSL::tanh( dual_deriv ) << endl;
  cout << "tanh( 3.0 ) = " << TSL::tanh( 3.0 ) << endl;

  // Test scalar operators
  cout << "dual = " << dual << endl;
  cout << "dual + p = " << dual + p << endl;
  cout << "p + dual = " << p + dual << endl;
  cout << "dual - p = " << dual - p << endl;
  cout << "dual * p = " << dual * p << endl;
  cout << "p * dual = " << p * dual << endl;
  cout << "dual / p = " << dual / p << endl;

	cout << "FINISHED" << endl;

}
