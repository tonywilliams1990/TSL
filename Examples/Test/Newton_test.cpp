// Test the Residual class
#include "Core"

using namespace std;

namespace TSL
{
  class test_residual : public Residual<double>
  {
    public:
			// The test equation is 2nd order
			test_residual() : Residual<double> ( 2 ) {}

      // Define the equation
			void residual_fn( const Vector<double>& x_k, Vector<double>& F ) const
			{
				F[ 0 ] = std::pow( x_k[ 0 ], 3.0 ) + x_k[ 1 ] - 1;
				F[ 1 ] = std::pow( x_k[ 1 ], 3.0 ) - x_k[ 0 ] + 1; 
				/*
					x^3 + y - 1 = 0,
					y^3 - x + 1 = 0,
					(x,y) = (1,0) is the only solution 
				*/
			}		
  }; // End of class test_residual

  // Define the residual for arc-length continuation of a circle
  class Arc_problem : public Residual<double>
  {
    public:
      double p;

      Arc_problem() : Residual<double>( 1 ) {}

      void residual_fn( const Vector<double>& x, Vector<double>& F ) const
      {
        F[ 0 ] = x[ 0 ] * x[ 0 ] + p * p - 2.0; 
      }       

  }; // End of class Arc_problem

} // End of namespace TSL

int main()
{
  cout << "----- TESTING Newton -----" << endl;

  /* ----- TESTING Newton class ----- */	

  TSL::test_residual Res;                                   // Create residual

  TSL::Vector<double> x_0(2, 0.5);                          // Initial guess
  x_0[1] = 0.25;	
  cout << "x_0 = " << endl << x_0 << endl;

  TSL::Newton<double> newton( &Res );                       // Test constructor

  newton.solve( x_0 );                                      // Test solve method 
  //newton.iterate( x_0 );                                    // Test iterate method
  
  cout << "Solution = " << endl << x_0 << endl; 

  // Arc-length continuation of a circle
  cout << endl << "----- Arc continuation of x^2 + p^2 = 2.0 -----" << endl;

  // Instantiate the problem
  TSL::Arc_problem arc_prob;
  arc_prob.p = 1.0;

  const double tol = 1.e-10;
  // Scalar newton iteration problem
  TSL::Newton<double> newton_arc( &arc_prob, 8, tol );
  newton_arc.set_monitor_det( true );                   // Not implemented yet

  // Initial guess
  TSL::Vector<double> state( 1, 1.0 );
  newton_arc.init_arc( state, &arc_prob.p, 0.001, 0.1 );

  do
  {
    newton_arc.arclength_solve( state );
    cout << "x = " << state[ 0 ] << ", p = " << arc_prob.p << endl;
  }while( state[ 0 ] < 1.0 );

  //TODO should try plotting some of the points 

	cout << "FINISHED" << endl;

}
