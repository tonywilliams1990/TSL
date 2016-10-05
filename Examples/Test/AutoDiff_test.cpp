// Test the AutoDiff class
#include "Core"

using namespace std;

template <class F>
F test_function( const F& x )
{
  return TSL::pow( x, 3.0 );
}

template <class F>
F test_function_2( const F& x )
{
  return x*x + 3.0*x;
}

template <class F>
F test_grad_func( const TSL::Vector<F>& x )
{
  return x[0] * x[1] * x[2];
}


int main()
{
  cout << "----- TESTING AutoDiff -----" << endl;

  double (*function_ptr)(const double&) = test_function_2;
  cout << "5*5 + 3*5 = " << function_ptr( 5.0 ) << endl;

  // Test constructor
  TSL::AutoDiff<double> diff;

  // Test derivative function (single input)
  double f_dash;

  f_dash = diff.derivative( &test_function, 1.3 );
  cout << "f' = " << f_dash << endl;

  f_dash = diff.derivative( &test_function_2, 2.4 );
  cout << "f' = " << f_dash << endl;

  // Test derivative function (vector input)
  TSL::Vector<double> x_vals;
  x_vals.linspace( 0.0, 1.0, 11 );
  //cout << "x_vals = " << endl << x_vals << endl;
  
  TSL::Vector<double> f_deriv;
  f_deriv = diff.derivative( &test_function, x_vals );
  //cout << "f_deriv = " << endl << f_deriv << endl;
  
  // Test the gradient function (single input)
  TSL::Vector<double> x_vec( 3, 1.0 );
  x_vec[ 1 ] = -2.0; x_vec[ 2 ] = 3.0;
  cout << "x_vec = " << endl << x_vec << endl;
  TSL::Vector<double> grad_f( 3, 0.0 );
  grad_f = diff.gradient( &test_grad_func, x_vec );
  cout << "grad_f = " << endl << grad_f << endl;

	cout << "FINISHED" << endl;

}
