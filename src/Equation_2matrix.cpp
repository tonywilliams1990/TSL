/// Implementation for an equations class that can be inherited from
/// to allow instantiation of ODE objects using the resulting class.

#include <Equation_2matrix.h>
#include <Residual_with_coords.h>
//#include <Utility.h>

namespace TSL
{

  template <class T, class X>
  Equation_2matrix<T, X >::Equation_2matrix( const unsigned& order ) :
      Equation_1matrix<T, X>( order )
  {
    // initialise the container for the extra matrix
    MATRIX1_AT_LAST_STATE = Matrix<T>( order, order, 0.0 );
    // add an extra coordinate to the vector stored in the residual_with_coords baseclass 
    Residual_with_coords<T,X>::coords.resize( 2, 0.0 );
  }

  template <class T, class X>
  Equation_2matrix<T, X >::~Equation_2matrix()
  {
    // timer reporting is done via the Equation (base) class
  }

  template <class T, class X>
  void Equation_2matrix<T, X >::update( const Vector<T> &state )
  {
    // call the base class's update method
    Equation_1matrix<T, X>::update( state );
    // now deal with the additional matrix separately
    matrix1( state, MATRIX1_AT_LAST_STATE );
  }

  template <class T, class X>
  void Equation_2matrix<T, X>::get_jacobian_of_matrix1_mult_vector( 
        const Vector<T> &state, const Vector<T> &vec, Matrix<T> &h ) const
  {
    // we dont need state in the default implementation as its already been set by
    // the update method. You do need it for the user to overload this method with an      
    // explicit analytical version however.
    //
    // copy some items for FD computation of Jacobian of mass matrix
    Vector<T> copy_of_state( this -> LAST_STATE );
    Matrix<T> copy_of_matrix( MATRIX1_AT_LAST_STATE );
    std::vector< Matrix<T> > jacmatrix;
    // update the Jacobian of the mass matrix
    for ( std::size_t i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      copy_of_state[ i ] += this -> DELTA;
      matrix1( copy_of_state, copy_of_matrix );
      copy_of_state[ i ] -= this -> DELTA;
      copy_of_matrix -= MATRIX1_AT_LAST_STATE;
      copy_of_matrix = copy_of_matrix * ( 1. / this -> DELTA );
      // the 3D object that represents the Jacobian of the mass matrix
      jacmatrix.push_back( copy_of_matrix );
    }

    // evaluate the jacabian of mass contribution
    for ( unsigned i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      for ( unsigned j = 0; j < this -> ORDER_OF_SYSTEM; ++j )
      {
        h( i, j ) = jacmatrix[ j ][ i ].dot( vec );
      }
    }
  }

  // the required templated versions are:
  template class Equation_2matrix<double>
  ;
  template class Equation_2matrix<std::complex<double> >
  ;
  template class Equation_2matrix<std::complex<double>, std::complex<double> >
  ;

} // end namespace
