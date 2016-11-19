/* Matrix class - Encapsulates the Eigen matrix container
*/
#include "Matrix.h"

namespace TSL
{

	

	// Destructor
	//template <class T>
  //Matrix<T>::~Matrix( )
  //{}


  // Templated versions
  template class Matrix<double>;
	template class Matrix< std::complex<double> >;

} // End of namespace TSL
