/* AutoDiff - A class for performing forward mode automatic differentiation
*/

#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <iostream>

#include "Vector.h"
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
          for (std::size_t i=0; i<x_vals.size(); ++i)
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
          Vector<T> eps( x_val.size(), 0.0 );
          Dual<T> result;
          Dual<T> dual( x_val[ 0 ], eps );
          Vector< Dual<T> > dual_vec( x_val.size(), dual );
          // Fill the vector of dual numbers
          for ( std::size_t i=0; i < x_val.size(); ++i )
          {
            dual_vec[ i ][ 0 ] = x_val[ i ];
            dual_vec[ i ].epsilons()[ i ] = 1.0; 
          }
          result = func_ptr( dual_vec );
          return result.epsilons();
        }
   

	}; // End of class AutoDiff

} // End of namespace TSL

#endif
