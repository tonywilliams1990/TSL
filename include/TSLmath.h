/* TSLmath - Specialise the cmath fuctions for use with dual number types.

*/

#ifndef TSLMATH_H
#define TSLMATH_H

#include <cmath>

#include "Error.h"
#include "Dual.h"

namespace TSL
{
	/* ----- Trigonometric functions ----- */
  
  /// Cos
  template <class T>
  Dual<T> cos( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::cos( dual.real() ), -std::sin( dual.real() ) * eps );
    return result;
  }

  template<class T>
  T cos( const T& x ) { return std::cos( x ); }

  /// Sin
  template <class T>
  Dual<T> sin( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::sin( dual.real() ), std::cos( dual.real() ) * eps );
    return result;
  }

  template<class T>
  T sin( const T& x ) { return std::sin( x ); }

  /// Tan
  template <class T>
  Dual<T> tan( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::tan( dual.real() ), ( 1.0 / ( std::cos( dual.real() ) * 
                                                       std::cos( dual.real() ) ) ) * eps );
    return result;
  }
  
  template<class T>
  T tan( const T& x ) { return std::tan( x ); }


  /* ----- Hyperbolic functions ----- */

  /// Cosh
  template <class T>
  Dual<T> cosh( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::cosh( dual.real() ), std::sinh( dual.real() ) * eps );
    return result;
  }

  template<class T>
  T cosh( const T& x ) { return std::cosh( x ); }
	
  /// Sinh
  template <class T>
  Dual<T> sinh( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::sinh( dual.real() ), std::cosh( dual.real() ) * eps );
    return result;
  }

  template<class T>
  T sinh( const T& x ) { return std::sinh( x ); }

  /// Tanh
  template <class T>
  Dual<T> tanh( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::tanh( dual.real() ), ( 1.0 / ( std::cosh( dual.real() ) * 
                                                     std::cosh( dual.real() ) ) ) * eps );
    return result;
  }
  
  template<class T>
  T tanh( const T& x ) { return std::tanh( x ); }

  /* ----- Exponential and logarithmic functions ----- */

  /// Exp
  template <class T>
  Dual<T> exp( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::exp( dual.real() ), std::exp( dual.real() ) * eps );
    return result;
  }

  template<class T>
  T exp( const T& x ) { return std::exp( x ); }

  /// Log
  template <class T>
  Dual<T> log( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::log( dual.real() ), ( 1.0 / dual.real() ) * eps );
    return result;
  }

  template<class T>
  T log( const T& x ) { return std::log( x ); }

  /* ----- Power functions ----- */

  /// Pow
  template <class T>
  Dual<T> pow( const Dual<T>& dual, const T& p )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::pow( dual.real(), p ), p* std::pow( dual.real(), p - 1.0 )* eps );
    return result;
  }

  template<class T>
  T pow( const T& x, const T& p ) { return std::pow( x, p ); }

  /// Sqrt
  template <class T>
  Dual<T> sqrt( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::sqrt( dual.real() ), ( 1.0 / ( 2.0 * std::sqrt( dual.real() ) ) )
                                                                                 * eps );
    return result;
  }

  template<class T>
  T sqrt( const T& x ) { return std::sqrt( x ); }

  /* ----- Other functions ----- */

  /// Abs
  template <class T>
  Dual<T> abs( const Dual<T>& dual )
  {
    Vector<T> eps;
    eps = dual.epsilons();
    Dual<T> result( std::abs( dual.real() ), ( dual.real() / std::abs( dual.real() ) )
                                                                                  * eps );
    return result;
  }

  template<class T>
  T abs( const T& x ) { return std::abs( x ); }


} // End of namespace TSL

#endif
