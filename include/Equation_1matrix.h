/// A templated class for equations that can be inherited from
/// to allow instantiation of PDE_IBVP objects.

#ifndef EQUATION_1MATRIX_H
#define EQUATION_1MATRIX_H

#include <Residual_with_coords.h>

namespace TSL
{

  /// An equation object base class used in the IBVP classes (and others).
  /// An equation object is essentially a ('square') residual object with an
  /// independent variable data member and access methods. By 'square' we mean that
  /// it defines N residuals and N state variables. The equation is defined
  /// using an NxN matrix that multiplies the derivative of unknowns and a residual RHS.
  template < class T, class X = double >
  class Equation_1matrix : public Residual_with_coords<T, X>
  {
  public:
    /// Constructor for equation class.
    /// \param order The order of the system
    explicit Equation_1matrix( const unsigned &order );

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation_1matrix();

    /// Update the Equation object for the current set of state variables
    /// \param state The state vector at which to set the equation object
    void update( const Vector<T> &state );

    /// Return a handle to the matrix
    const Matrix<T>& matrix0() const;

    /// Return the product of the Jacobian-of-the-matrix and a vector 'vec'
    /// when the equation has a given 'state'. The user should overload this
    /// if concerned about performance of the solver. If not overloaded, the
    /// default is to finite difference the Jacobian-of-the-matrix.
    virtual void get_jacobian_of_matrix0_mult_vector( const Vector<T> &state, 
                                    const Vector<T> &vec, Matrix<T> &h ) const;


  protected:

    /// Define the matrix in terms of the current state vector.
    virtual void matrix0( const Vector<T> &x, Matrix<T> &m ) const
    {
      std::string problem;
      problem = "The equation::matrix0 method has not been implemented!\n";
      problem += "You have to implement this method to define the equation.\n";
      throw Error( problem );
    }

  private:
    /// Matrix0 evaluated for the last state vector
    Matrix<T> MATRIX0_AT_LAST_STATE;

  }
  ; // End of class Equation_1matrix

  template <class T, class X>
  inline const Matrix<T>& Equation_1matrix<T, X>::matrix0() const
  {
    return MATRIX0_AT_LAST_STATE;
  }

} // End namespace TSL

#endif
