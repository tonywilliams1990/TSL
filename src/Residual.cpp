/* Residual - A residual class for use with Newton for solving systems of 
			   non-linear equations of the form F(x) = 0 where x is a vector 	
			   and F is vector valued function. Given a current approximation x_k
			   of x we can calculate the residual and approximate the Jacobian.
*/

#include "Residual.h"

namespace TSL
{

  template <typename _Type>
  void Residual<_Type>::jacobian( const Vector<_Type>& state, Matrix<_Type>& jac ) const
  {
    Vector<_Type> new_state( state );
    // evaluation of the function
    Vector<_Type> f_at_new_state( ORDER_OF_SYSTEM, 0.0 );
    // default is to FD the Jacobian
    for ( std::size_t i = 0; i < NUMBER_OF_VARS; ++i )
    {
      new_state[ i ] += DELTA;
      residual_fn( new_state, f_at_new_state );
      new_state[ i ] -= DELTA;
      jac.set_col( i, ( f_at_new_state - FN_AT_LAST_STATE ) / DELTA );
    }
  }

  // Templated versions
  template class Residual<double>;
	template class Residual< std::complex<double> >;

}
