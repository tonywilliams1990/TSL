/* Newton - A class for using Newton's method for solving systems of nonlinear
			 equations of the form F(x) = 0 where x is a vector and F is vector
			 valued function. 
*/

#ifndef NEWTON_H
#define NEWTON_H

#include <iostream>

#include "Matrix.h"
#include "Vector.h"
#include "Error.h"
#include "Timer.h"
#include "Residual.h"

namespace TSL
{
  

	/// A templated Newton iteration class
	template <class T>
	
	class Newton 
	{
    protected:
      Residual<T>* ptr_RESIDUAL;                // Pointer to the residual object
      std::size_t MAX_ITER;                     // Maximum number of iterations
      double TOL;                               // Convergence tolerance
      double DELTA;                             // Derivative step
      std::size_t ORDER;                        // Order of the system

    public:
      /// Constructor
      Newton( Residual<T>* ptr_residual, std::size_t max_iter = 20, 
              double tolerance = 1.0e-8, double delta_step = 1.0e-8 ) :
      TOL( tolerance ),
      MAX_ITER( max_iter ),
      ptr_RESIDUAL( ptr_residual )
      {
        ptr_RESIDUAL -> delta() = delta_step;
        DELTA = delta_step;
        ORDER = ptr_RESIDUAL -> get_order();
      }

      /// Destructor
      ~Newton() {}

      /* ----- Methods ----- */

      /// Newton iterate using intial guess x_0 (stores the solution in x_0)
      void iterate( Vector<T>& x_0 )
      {
        if ( x_0.size() != ORDER ) { throw Error( "Newton error: size does not agree " ) ; }
        Vector<T> F( ORDER, 0.0 );        // Residual function evaluation
        Matrix<T> J( ORDER, ORDER, 0.0 ); // Jacobian matrix
        Vector<T> dx( ORDER, 0.0 );       // Increment vector
        std::size_t iter( 0 );            // Iteration counter
        double max_residual;              // Maximum residual

        do
        {
          ++iter;
          ptr_RESIDUAL -> update( x_0 );  // Update the residual object
          F = ptr_RESIDUAL -> residual(); // Get the residual
          max_residual = F.norm_inf();    // Find the maximum residual
          J = ptr_RESIDUAL -> jacobian(); // Get the Jacobian matrix
          
          dx = J.solve( -F );             // Solve for the increment
          x_0 += dx;                      // Increment the vector 
        }while( iter < MAX_ITER && max_residual > TOL );

        if ( iter == MAX_ITER )
        {
          std::cout << "NEWTON WARNING! MAXIMUM NUMBER OF ITERATIONS REACHED." << std::endl;
        } 
      }

      /// Solve the system using Newton iteration (points to iterate)
      void solve( Vector<T>& x_0 )
      {
        iterate( x_0 );
      }
		
      

	}; // End of class Newton


} // End of namespace TSL

#endif
