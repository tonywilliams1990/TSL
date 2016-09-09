/* Residual - A residual class for use with Newton for solving systems of 
			   non-linear equations of the form F(x) = 0 where x is a vector 	
			   and F is vector valued function. Given a current approximation x_k
			   of x we can calculate the residual and approximate the Jacobian.
*/

#ifndef RESIDUAL_H
#define RESIDUAL_H

#include <iostream>

#include "Matrix.h"
#include "Vector.h"
#include "Error.h"
#include "Timer.h"

namespace TSL
{
  

	/// A templated base class to be inherited by objects that define residuals
	template <class T>
	
	class Residual 
	{
		
      protected:
        /// Because the residual evaluation at the current state is assumed
        /// to have already been done by the 'update' method, this routine is
        /// protected. This default uses a finite-differenced Jacobian.
        /// You can overload this to provide an analytic Jacobian if you wish
        virtual void jacobian( const Vector<T>& state, Matrix<T>& jac ) const;

        /// Jacobian for the last state vector
        Matrix<T> JAC_AT_LAST_STATE;
        /// Residual for the last state vector
        Vector<T> FN_AT_LAST_STATE;
        /// The last state vector
        Vector<T> LAST_STATE;
        /// A default step for FD computation of the Jacobian
        T DELTA;
        /// The order of the system of equations
        std::size_t ORDER_OF_SYSTEM;
        /// The number of elements in the state vector
        std::size_t NUMBER_OF_VARS;
      #ifdef TIME
        Timer T_UPDATER;
      #endif

      public:
      /// Constructor for a 'square' residual object (N residuals for N unknowns)
      Residual( const std::size_t& order ) : DELTA( 1.e-8 ) 
      {
        ORDER_OF_SYSTEM = order;
        NUMBER_OF_VARS = order;
        LAST_STATE = Vector<T>( NUMBER_OF_VARS, 0.0 );
        FN_AT_LAST_STATE = Vector<T>( ORDER_OF_SYSTEM, 0.0 );
        JAC_AT_LAST_STATE = Matrix<T>( ORDER_OF_SYSTEM, NUMBER_OF_VARS, 0.0 );
      #ifdef TIME
        T_UPDATER = Timer( "Updating of the residual object:" );
      #endif
      }

      /// Constructor for a 'non-square' residual object 
      Residual( const std::size_t& order, const std::size_t& nvars ) : DELTA( 1.e-8 )
      {
        ORDER_OF_SYSTEM = order;
        NUMBER_OF_VARS = nvars;
        LAST_STATE = Vector<T>( NUMBER_OF_VARS, 0.0 );
        FN_AT_LAST_STATE = Vector<T>( ORDER_OF_SYSTEM, 0.0 );
        JAC_AT_LAST_STATE = Matrix<T>( ORDER_OF_SYSTEM, NUMBER_OF_VARS, 0.0 );
      #ifdef TIME
        T_UPDATER = Timer( "Updating of the residual object:" );
      #endif
      }

      /// Destructor ( virtual since we have virtual functions )
      virtual ~Residual() 
      { 
        #ifdef TIME
        std::cout << std::endl;
        T_UPDATER.stop();
        T_UPDATER.print();
      #endif
      }

      
      /* ----- Methods ----- */

      /// Get the order of the residual vector
      std::size_t get_order() const { return ORDER_OF_SYSTEM; } 

      /// Get the number of variables
      std::size_t get_number_of_vars() const { return NUMBER_OF_VARS; }

      /// Update the Residual object for the current set of state variables
      void update( const Vector<T>& state )
      {
      #ifdef TIME
        T_UPDATER.start();
      #endif
        LAST_STATE = state;
        residual_fn( LAST_STATE, FN_AT_LAST_STATE );
        jacobian( LAST_STATE, JAC_AT_LAST_STATE );
      #ifdef TIME
        T_UPDATER.stop();
      #endif
      }

      /// Return a handle to the residuals corresponding to the last update state
      const Vector<T>& residual() const { return FN_AT_LAST_STATE; }

      /// Return a handle to the Jacobian of the residual
      const Matrix<T>& jacobian() const { return JAC_AT_LAST_STATE; }

      /// Return a handle to the state vector
      const Vector<T>& state() const { return LAST_STATE; }

      /// Return a handle to the step size used when finite-differencing
      T& delta() { return DELTA; }

      /// Return a handle to the step size used when finite-differencing
      const T& delta() const { return DELTA; }

      virtual void residual_fn( const Vector<T>& state, Vector<T>& f ) const
      {
        std::string problem;
        problem = "The Residual::residual_fn method has not been implemented.\n";
        problem += "You have to implement this method to define the residual.\n";
        throw Error( problem );
      }

	}; // End of class Residual


} // End of namespace TSL

#endif
