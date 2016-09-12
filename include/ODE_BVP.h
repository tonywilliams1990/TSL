/* ODE_BVP - Here we define the ODE_BVP class used for solving ODE
			       boundary value problems. 
*/

#ifndef ODE_BVP_H
#define ODE_BVP_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Vector.h"
#include "Error.h"
#include "Matrix.h"
#include "SparseMatrix.h"
#include "Equation.h"
#include "OneD_node_mesh.h"


namespace TSL
{
  

	/// A templated Newton iteration class
	template <class T, class X = double>
	
	class ODE_BVP
	{
    private:	
			Equation<T, X> *ptr_EQUATION; 	  // Pointer to ODE equation
      Residual<T> *ptr_LEFT_RESIDUAL;   // Pointer to the left residual
      Residual<T> *ptr_RIGHT_RESIDUAL;  // Pointer to the right residual
      OneD_node_mesh<T, X> SOLUTION;    // Solution mesh
      std::size_t MAX_ITER;			        // Maximum number of iterations
      double TOL;						            // Tolerance for convergence
      Vector<X> NODES;                  // Vector of nodes ( defines the domain ) 
      double DELTA;                     // Perturbation step for computation of the Jacobian
    
		public:	
			
      /* ----- Constructors and destructors ----- */

			// Constructor
      ODE_BVP( Equation<T, X>* ptr_equation, const Vector<X>& nodes, Residual<T>* ptr_left_residual, Residual<T>* ptr_right_residual );

      // Copy constructor
			ODE_BVP(  const ODE_BVP& source );

			// Destructor
      virtual ~ODE_BVP();
			
      /* ----- Methods ----- */

      // Return a pointer to DELTA
      double& delta();

      // Return a pointer to MAX_ITER
      std::size_t& max_iterations();

      // Return a pointer to TOL
      double& tolerance();

      // Return the current solution mesh
			OneD_node_mesh<T, X>& solution();
            
      // Solve the BVP - do iterate while residual > tol and iter < max_iter
			void solve();
    

	}; // End of class ODE_BVP


} // End of namespace TSL

#endif
