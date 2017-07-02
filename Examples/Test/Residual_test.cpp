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
} // End of namespace TSL

int main()
{
  cout << "----- TESTING Residual -----" << endl;

  /* ----- TESTING Residual class ----- */	

  TSL::test_residual Res;                                   // Test constructor
  cout << "Order = " << Res.get_order() << endl;            // Test get_order
  cout << "Vars = " << Res.get_number_of_vars() << endl;    // Test get_number_of_vars

  TSL::Vector<double> x_0(2, 0.5);                          // Initial guess
  x_0[1] = 0.25;	
  cout << "x_0 = " << endl << x_0 << endl;

  // Test function evaluation
	TSL::Vector<double> F( 2, 0.0 );	
	Res.residual_fn( x_0, F );
	cout << "F(x_0) = " << endl << F << endl;

  Res.update( x_0 );                                        // Test update method

  TSL::Vector<double> residual = Res.residual();            // Test residual method
  cout << "residual = " << endl << residual << endl;

  TSL::Matrix<double> jacobian = Res.jacobian();            // Test jacobian method
  cout << "jacobian = " << endl << jacobian << endl;

  TSL::Vector<double> state = Res.state();                  // Test state method
  cout << "state = " << endl << state << endl;

	cout << "FINISHED" << endl;

}
