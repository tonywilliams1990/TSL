/// A templated class for equations that can be inherited from
/// to allow instantiation of PDE_double_IBVP objects (amongst others).

#ifndef EQUATION_2MATRIX_H
#define EQUATION_2MATRIX_H

#include <Equation_1matrix.h>

namespace TSL
{

  /// An equation object base class used in the PDE_double_IBVP class.
  /// An equation object is essentially a ('square') residual object (although
  /// it doesn't currently inherit) with some independent variable
  /// data members and access methods. By 'square' we mean that it defines
  /// N residuals and N state variables. In this case the equation also defines
  /// 2 matrices (amongst other data). This inherits from the Equation_1matrix
  /// and adds the functionality for the additional matrix.

  template < class T, class X = double >
  class Equation_2matrix : public Equation_1matrix<T, X>
  {    
  public:

    /// Constructor for equation class.
    /// \param order The order of the system
    explicit Equation_2matrix( const unsigned &order );

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation_2matrix();

    /// Update the Equation object for the current set of state variables
    void update( const Vector<T> &state );

    /// Return a handle to the matrix member data
    const Matrix<T>& matrix1() const;

    /// Return the product of the Jacobian-of-the-matrix and a vector 'vec'
    /// when the equation has a given 'state'. The user should overload this
    /// if concerned about performance of the solver. If not overloaded, the
    /// default is to finite difference the Jacobian-of-the-matrix.
    virtual void get_jacobian_of_matrix1_mult_vector( const Vector<T> &state, 
                              const Vector<T> &vec, Matrix<T> &h ) const;

  protected:

    /// Define the matrix in terms of the current state vector.
    virtual void matrix1( const Vector<T> &state, Matrix<T> &m ) const
    {
      std::string problem;
      problem = "The equation::matrix1 method has not been implemented.\n";
      problem += "You have to implement this method to define the equation.\n";
      throw Error( problem );
    }

  private:

    /// Matrix for the last state vector
    Matrix<T> MATRIX1_AT_LAST_STATE;

  }
  ; // End of class Equation_2matrix
 
  template <class T, class X>
  inline const Matrix<T>& Equation_2matrix<T, X>::matrix1() const
  {
    return MATRIX1_AT_LAST_STATE;
  }

} // End of namespace TSL

#endif
