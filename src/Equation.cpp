/* Equation - Implementation for Equation class
*/

#include <Equation.h>
#include <Residual_with_coords.h>

namespace TSL
{

  template <class T, class X>
  Equation<T, X>::Equation( const unsigned& order ) :
      Residual_with_coords<T, X>( order, 1 )
  {
  }

  template <class T, class X>
  Equation<T, X>::~Equation()
  {
  }

  // the required templated versions are:
  template class Equation<double>
  ;
  template class Equation<std::complex<double> >
  ;
  template class Equation<std::complex<double>, std::complex<double> >
  ;

} // End of namespace TSL
