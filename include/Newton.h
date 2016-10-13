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
#include "Arclength.h"

namespace TSL
{
  

	/// A templated Newton iteration class
	template <class T>
	
	class Newton : public Arclength< T >
	{
    protected:
      Residual<T>* ptr_RESIDUAL;                // Pointer to the residual object
      std::size_t MAX_ITER;                     // Maximum number of iterations
      double TOL;                               // Convergence tolerance
      double DELTA;                             // Derivative step
      std::size_t ORDER;                        // Order of the system
      int LAST_DET_SIGN;                        // Last sign of the determinant of Jacobian
      bool MONITOR_DET;                         // True if determinant is to be monitored

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
        //TODO determinant monitoring etc
        if ( x_0.size() != ORDER ) { throw Error( "Newton error: size does not agree." ); }
        Vector<T> F( ORDER, 0.0 );                        // Residual function evaluation
        Matrix<T> J( ORDER, ORDER, 0.0 );                 // Jacobian matrix
        Vector<T> dx( ORDER, 0.0 );                       // Increment vector
        std::size_t iter( 0 );                            // Iteration counter
        double max_residual;                              // Maximum residual

        do
        {
          ++iter;
          ptr_RESIDUAL -> update( x_0 );                  // Update the residual object
          F = ptr_RESIDUAL -> residual();                 // Get the residual
          max_residual = F.norm_inf();                    // Find the maximum residual
          J = ptr_RESIDUAL -> jacobian();                 // Get the Jacobian matrix
          
          dx = J.solve( -F );                             // Solve for the increment
          x_0 += dx;                                      // Increment the vector 
        }while( iter < MAX_ITER && max_residual > TOL );

        if ( iter == MAX_ITER )
        {
          std::cout << "NEWTON WARNING! MAXIMUM NUMBER OF ITERATIONS REACHED." << std::endl;
        } 
      }

      /// If set the system will monitor the sign of the determinant
      void set_monitor_det( bool flag ) { MONITOR_DET = flag; }

      /// Solve the system using Newton iteration (points to iterate)
      void solve( Vector<T>& x_0 )
      {
        iterate( x_0 );
      }

      /// Arc-length solve the system (init_arc must be called first)
      void arclength_solve( Vector<T>& x )
      {
        // Throw an error if init_arc has not been called
        if ( !this->INITIALISED ) { throw Error( "Newton arc error: not initialised." ); }
        // Backup the state / parameter in case we fail
        Vector<T> backup_state( x );
        T backup_parameter( *(this->ptr_PARAM) );
        int det_sign( 1 );                                // Sign of the determinant
        bool step_succeeded( false );                     // Check for success
        std::size_t itn( 0 );                             // Iteration counter
        // Guess the next solution
        *( this->ptr_PARAM ) = this->LAST_PARAM + this->PARAM_DERIV_S * this->DS;
        x = this->LAST_X + this->X_DERIV_S * this->DS;

        Matrix<T> J( ORDER, ORDER, 0.0 );                 // Jacobian matrix
        Vector<T> R1( ORDER, 0.0 );
        Vector<T> R2( ORDER, 0.0 );
        Vector<T> dR_dp( ORDER, 0.0 );
        Vector<T> y( ORDER, 0.0 );
        Vector<T> z( ORDER, 0.0 );
        Vector<T> JacE; //TODO set order?
        
        do
        {       
          this->ptr_RESIDUAL -> update( x );              // Update the residual
          J = this->ptr_RESIDUAL -> jacobian();           // Get the Jacobian
          R1 = this->ptr_RESIDUAL -> residual();          // Get the residual
          double E1 = this -> arclength_residual( x );    // Arclength residual
          // Compute derivatives wrt the parameter
          *(this->ptr_PARAM) += this->DELTA;              
          this->ptr_RESIDUAL -> residual_fn( x, R2 );
          double E2 = this -> arclength_residual( x ); 
          *(this->ptr_PARAM) -= this->DELTA; 
          dR_dp = ( R2 - R1 ) / this->DELTA;
          T dE_dp = ( E2 - E1 ) / this->DELTA;
          // Bordering algorithm
          y = J.solve( -R1 );
          z = J.solve( -dR_dp );
          JacE = this-> Jac_arclength_residual( x );
          T delta_p = -( E1 + JacE.dot( y ) ) / ( dE_dp + JacE.dot( z ) );
          Vector<T> delta_x = y + z * delta_p;
          double max_correction = std::max( delta_x.norm_inf(), std::abs( delta_p ) );
          if ( max_correction < this->TOL )
          {
            step_succeeded = true;
            break;
          }
          // add the corrections to the state variables
          x += delta_x;
          *(this->ptr_PARAM) += delta_p;
          ++itn;
          if ( itn > MAX_ITER )
          {
            step_succeeded = false;
            break;
          }
        }while( true );
        
        // If this isn't a successful step
        if ( !step_succeeded )
        {
          x = backup_state;                                 // Restore state
          *(this->ptr_PARAM) = backup_parameter;            // Restore the parameter
          this->ptr_RESIDUAL -> update( x );                // Restore the residual
          this -> DS /= this -> ARCSTEP_MULTIPLIER;         // Reduce our step length

        }
        else // If this is a successful step
        { 
          this -> update( x );
          // TODO Bifurcation detection using sign of determinant (Newton.cpp ln 233 ...)
          if ( itn >= 7 )                                   // Converging to slowly
          {
            this -> DS /= this -> ARCSTEP_MULTIPLIER;
          }
          if ( itn <= 2 )                                   // Converging to quickly
          {
            if ( std::abs( this -> DS * this -> ARCSTEP_MULTIPLIER ) < this -> MAX_DS )
            {
              this -> DS *= this -> ARCSTEP_MULTIPLIER;
            } 
          }
        }

      
		  }
	}; // End of class Newton


} // End of namespace TSL

#endif
