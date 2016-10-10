/* AutoDiff - A class for performing forward mode automatic differentiation
*/

#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <iostream>

#include "Vector.h"
#include "Matrix.h"
#include "Residual.h"
#include "Error.h"
#include "Dual.h"

namespace TSL
{

  template<class T, class F = Dual<T> >
	
	class AutoDiff 
	{
		
      protected:
        

      public:
        
        /// Default constructor
        AutoDiff() { }
        
        /// Destructor
        ~AutoDiff() { }

        /* ----- Methods ----- */
  
        /// Calculate the derivative of a scalar function at a point (x_val) 
        T derivative( F (*func_ptr)(const F&), const T& x_val )
        {
          Vector<T> eps( 1, 1.0 );
          Dual<T> dual( x_val, eps );
          Dual<T> result;
          result = func_ptr( dual );
          return result.epsilons()[ 0 ];
        } 

        /// Calculate the derivative of a scalar function at points stored in a Vector 
        Vector<T> derivative( F (*func_ptr)(const F&), const Vector<T>& x_vals )
        {
          Vector<T> eps( 1, 1.0 );
          Dual<T> dual( x_vals[ 0 ], eps );
          Dual<T> result;
          Vector<T> ans( x_vals.size() );
          for ( std::size_t i=0; i<x_vals.size(); ++i )
          {
            dual[ 0 ] = x_vals[ i ];
            result = func_ptr( dual );
            ans[ i ] = result.epsilons()[ 0 ];
          }
          return ans;
        }

        /// Calculate the gradient of a scalar field at a point
        Vector<T> gradient( F (*func_ptr)(const Vector<F>&), const Vector<T>& x_val )
        {
          const std::size_t n( x_val.size() );
          Vector<T> eps( n, 0.0 );
          Dual<T> result;
          Dual<T> dual( x_val[ 0 ], eps );
          Vector< Dual<T> > dual_vec( n, dual );
          // Fill the vector of dual numbers
          for ( std::size_t i=0; i < n; ++i )
          {
            dual_vec[ i ][ 0 ] = x_val[ i ];
            dual_vec[ i ].epsilons()[ i ] = 1.0; 
          }
          result = func_ptr( dual_vec );
          return result.epsilons();
        }

        /// Calculate the gradient of a scalar field at points stored in a vector
        Vector< Vector<T> > gradient( F (*func_ptr)(const Vector<F>&), 
                                      const Vector< Vector<T> >& x_vals )
        {
          Vector< Vector<T> > ans( x_vals.size() );
          for ( std::size_t i=0; i<x_vals.size(); ++i )
          {
            ans[ i ] = gradient( func_ptr, x_vals[ i ] ); 
          }
          return ans;
        }

        /// Calculate the Jacobian of a vector valued function at a point
        Matrix<T> jacobian( Vector<F> (*func_ptr)(const Vector<F>&), 
                            const Vector<T>& x_val )
        {
          const std::size_t n( x_val.size() );
          Vector<T> eps( n, 0.0 );
          Dual<T> dual( x_val[ 0 ], eps );
          Vector< Dual<T> > dual_vec( n, dual );
          // Fill the vector of dual numbers
          for (std::size_t i=0; i<n; ++i )
          {
            dual_vec[ i ][ 0 ] = x_val[ i ];
            dual_vec[ i ].epsilons()[ i ] = 1.0;
          }
          Vector< Dual<T> > result;
          result = func_ptr( dual_vec );
          const std::size_t m( result.size() );
          // Fill the matrix
          Matrix<T> temp( m, n, x_val[0] );
          for ( std::size_t i=0; i < m; ++i )
          {
            for ( std::size_t j=0; j< n; ++j )
            {
              temp( i, j ) = result[ i ].epsilons()[ j ];
            }
          }
          return temp;
        }
        
        ///TODO Calculate the Jacobian at points stored in a vector
        

        //TODO Use dual numbers in Residual class to calculate the Jacobian etc?

	}; // End of class AutoDiff

} // End of namespace TSL

#endif
