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
				F[ 0 ] = pow( x_k[ 0 ], 3.0 ) + x_k[ 1 ] - 1;
				F[ 1 ] = pow( x_k[ 1 ], 3.0 ) - x_k[ 0 ] + 1; 
				/*
					x^3 + y - 1 = 0,
					y^3 - x + 1 = 0,
					(x,y) = (1,0) is the only solution 
				*/
			}		
  }; // End of class test_residual
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

	cout << "FINISHED" << endl;

}
