/* Arclength - A base class for arc length solvers.		  
*/

#include "Arclength.h"

namespace TSL
{

  /* ----- Protected methods ----- */

  /// A method which stores the current converged state and parameter + computes derivs
  template <class T>
  void Arclength<T>::update( const Vector<T>& x )
  {
    if ( RESCALE_THETA ) { update_theta( x ); }
    X_DERIV_S = ( x - LAST_X ) / DS;
    PARAM_DERIV_S = ( *ptr_PARAM - LAST_PARAM ) / DS;
    LAST_X = x;
    LAST_PARAM = *ptr_PARAM;
  }

  /// Extra constraint that is to be used to replace the unknown arc-length
  template <class T>
  double Arclength<T>::arclength_residual( const Vector<T>& x ) const
  {
    return THETA * ( x - LAST_X ).norm_2() / x.size()
            + ( 1.0 - THETA ) * std::pow( std::abs( *ptr_PARAM - LAST_PARAM ), 2 )
            - DS * DS;
  }

  /// Automatically update theta value only if RESCALE_THETA = true
  template <class T>  
  void Arclength<T>::update_theta( const Vector<T>& x )
  {
    if ( RESCALE_THETA )
    {
      double Delta_p2 = std::pow( std::abs( *ptr_PARAM - LAST_PARAM ), 2);
      double Delta_x2 = ( x - LAST_X ).norm_2() / x.size();
      THETA = Delta_p2 * ( DESIRED_ARC_PROPORTION - 1.0 ) /
             ( Delta_p2 * ( DESIRED_ARC_PROPORTION - 1.0 ) 
              - DESIRED_ARC_PROPORTION * Delta_x2 );
    }
  }

  /* ----- Public methods ----- */

  /// Initialise the arc-length continuation class
  template <class T>
  void Arclength<T>::init_arc( Vector<T> x, T* ptr_param, 
                               const double& ds, const double& max_ds )
  {
    ptr_PARAM = ptr_param;
    // Compute the solution at this parameter value in the usual way
    solve( x ); 
    // Store the converged solution at this point                        
    LAST_X = x;
    LAST_PARAM = *ptr_PARAM;
    // We now have one solution which we can use to arc-length continue
    DS = ds;
    *ptr_PARAM += DS;
    // Recompute the state at this new parameter value
    solve( x );
    // Update the derivatives of state and parameter variables
    update( x );
    INITIALISED = true;
    MAX_DS = max_ds;
  } 


  template <class T>
  void Arclength<T>::solve( Vector<T>& x ) 
  {}

  
  

  // Templated versions
  template class Arclength<double>;
	template class Arclength< std::complex<double> >;

}
