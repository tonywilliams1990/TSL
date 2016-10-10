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
  if ( x.size() != 3 ){ throw TSL::Error( "test_grad_func input size error" );}
  return x[0] * x[1] * x[2];
}

template <class F>
TSL::Vector<F> test_jac_func( const TSL::Vector<F>& x )
{
  if ( x.size() != 2 ){ throw TSL::Error( "test_jac_func input size error" );}
  TSL::Vector<F> temp( 3 );
  temp[ 0 ] = x[ 0 ] * x[ 0 ] * x[ 1 ];
  temp[ 1 ] = 5.0 * x[ 0 ] + 2.0 * x[ 1 ];
  temp[ 2 ] = 3.0 * x[ 0 ];
  // f(x,y) = [ x^2 * y, 5*x + 2*y, 3*x ]^T
  return temp;
}

namespace TSL
{
  template <class T>

  class test_residual : public Residual<T>
  {
    public:
			// The test equation is 2nd order
			test_residual() : Residual<T> ( 2 ) {}

      // Define the equation
			void residual_fn( const Vector<T>& x_k, Vector<T>& f ) const
			{
				f[ 0 ] = TSL::pow( x_k[ 0 ], 3.0 ) + x_k[ 1 ] - 1;
				f[ 1 ] = TSL::pow( x_k[ 1 ], 3.0 ) - x_k[ 0 ] + 1; 
				/*
					x^3 + y - 1 = 0,
					y^3 - x + 1 = 0,
					(x,y) = (1,0) is the only solution 
				*/
			}		
  }; // End of class test_residual
} // End of namespace TSL

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

  // Test the gradient function (vector of inputs)
  TSL::Vector< TSL::Vector<double> > x_vecs( 3, x_vec );
  TSL::Vector< TSL::Vector<double> > grads( 3 );
  grads = diff.gradient(  &test_grad_func, x_vecs );

  for ( std::size_t i=0; i<x_vecs.size(); ++i )
  {
    cout << "grad[" << i << "] = " << endl << grads[ i ] << endl;
  }

  // Test the Jacobian function (single input)
  TSL::Vector<double> x_jac( 2 );
  x_jac[ 0 ] = 2.0; x_jac[ 1 ] = 3.0;
  cout << "x_jac = " << endl << x_jac << endl;
  cout << "f( x_jac ) = " << endl << test_jac_func( x_jac ) << endl;
  TSL::Matrix<double> jac( 3, 2, 0.0);
  jac = diff.jacobian( &test_jac_func, x_jac );
  cout << "J = " << endl << jac << endl;

  // Test the Jacobian function (vector of inputs)
  
  

	cout << "FINISHED" << endl;

}
