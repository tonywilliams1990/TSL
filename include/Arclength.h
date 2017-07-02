/* Arclength - A base class for arc length solvers.		  
*/

#ifndef ARCLENGTH_H
#define ARCLENGTH_H

#include <iostream>

#include "Vector.h"
#include "Error.h"


namespace TSL
{
  

	/// A templated arc length solver class
	template <class T>
	
	class Arclength 
	{
    protected:
      T *ptr_PARAM;                         // Pointer to the arc-length parameter
      Vector<T> LAST_X;                     // State variable at the last computed solution
      Vector<T> X_DERIV_S;                  // Deriv of the state variable wrt arc-length
      T LAST_PARAM;                         // Parameter value at the last computed solution
      T PARAM_DERIV_S;                      // Deriv of the parameter wrt arc-length
      double DS;                            // Size of the arc-length step
      double MAX_DS;                        // Maximum arc length step to be taken
      double ARCSTEP_MULTIPLIER;            // Step change multiplier
      bool INITIALISED;                     // True if the arc-length solver is initialised
      double THETA;                         // Theta weighting paramter (usually = 1/2)

      /* ----- Methods ----- */

      /// A method which stores the current converged state and parameter + computes derivs
      void update( const Vector<T>& x );

      /// Extra constraint that is to be used to replace the unknown arc-length
      double arclength_residual( const Vector<T>& x ) const;

      /// Return the derivative of the arclength_residual function wrt each state variable
      Vector<T> Jac_arclength_residual( Vector<T>& x ) const
      {
        Vector<T> Jx( x - LAST_X );
        Jx = Jx * ( THETA / ( x.size() * ( x - LAST_X ).norm_2() ) );
        return Jx;
      }

      /// Automatically update theta value only if RESCALE_THETA = true
      void update_theta( const Vector<T>& x );

    private:
      double DESIRED_ARC_PROPORTION;        // Desired proportion of the arc-length
      bool RESCALE_THETA;                   // If true theta is rescaled at each step

    public:

      /// Constructor
      Arclength() : ARCSTEP_MULTIPLIER( 2.0 ),
                    INITIALISED( false ),
                    THETA( 0.5 ),
                    DESIRED_ARC_PROPORTION( 0.5 ),
                    RESCALE_THETA( false )
      {}

      /// Destructor
      virtual ~Arclength() {}

      /* ----- Methods ----- */

      /// Initialise the arc-length continuation class
      void init_arc( Vector<T> x, T* ptr_param, const double& ds, const double& max_ds ); 

      /// Solve the system for a given initial guess
      virtual void solve( Vector<T>& x ) = 0;

      /// Return a handle to the arclength step
      double& ds() { return DS; }

      /// Return a handle to the arc step multiplier
      double& arcstep_multiplier() { return ARCSTEP_MULTIPLIER; }

      /// Return a handle to the RESCALE_THETA flag 
      bool& rescale_theta() { return RESCALE_THETA; }
      
      /// Return a handle to the theta parameter
      double& theta() { return THETA; }

      /// Return a handle to the desired proportion of the parameter to be used
      double& desired_arc_proportion() { return DESIRED_ARC_PROPORTION; } 

      
      

	}; // End of class Arclength


} // End of namespace TSL

#endif
