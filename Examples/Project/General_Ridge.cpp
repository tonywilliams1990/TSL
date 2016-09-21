#include <cassert>
#include <cmath>

#include "Core"

// Enumerations
enum{ f, fd, fdd, g, gd, gdd };                         // ODE enumeration
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };           // Base flow enumeration

// Either Base_2D or Base_3D for 2D or 3D base flows
#define Base_2D

namespace TSL
{
    namespace Param
    {
      double hzeta_right( 20.0 );       // Size of the domain in the zeta_hat direction
      double eta_top( 30.0 );           // Size of the domain in the eta direction
      std::size_t N( 400 );             // Number of intervals in the zeta_hat direction
      std::size_t M( 400 );             // Number of intervals in the eta direction
      std::size_t Nvar( 4 );            // Number of variables
      double n( 0.1 );                  // Pressure gradient parameter
      double beta = (2.0*n)/(n+1.0);    // Hartree parameter

    } // End of namespace Param

    namespace Base_Flow
    {
      double K( 0.0 );                  // Base transpiration parameter (+ve = blowing)   
      double beta = Param::beta;        // Hartree parameter

#ifdef Base_2D
      class equation : public Equation<double>
      {
        public:
          // Falkner-Skan equation is 3rd order
          equation() : Equation<double> ( 3 ) {}
          // Define the equation
          void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
          {
            F[ f ]   = u[ fd ];
            F[ fd ]  = u[ fdd ];
            F[ fdd ] = - u[ f ] * u[ fdd ] - beta * ( 1.0 - u[ f ] * u[ f ] ); 
          }
      }; // End Falkner-Skan equation class

      class plate_BC : public Residual<double>
      {
        public:
          plate_BC() : Residual<double> ( 2, 3 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + K;
            B[ 1 ] = z[ fd ];
          }
      }; // End Falkner-Skan plate_BC class

      class far_BC : public Residual<double>
      {
        public:
          far_BC() : Residual<double> ( 1, 3 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ fd ] - 1.0;
          }
      }; // End Falkner-Skan far_BC class
#endif

#ifdef Base_3D
      class equation : public Equation<double>
      {
        public:
          // The 3D alternative equation is 6th order
          equation() : Equation<double> ( 6 ) {}
          // Define the equation
          void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
          {
            //TODO
          }
      }; // End 3D alternative equation class

      class plate_BC : public Residual<double>
      {
        public:
          plate_BC() : Residual<double> ( 4, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            //TODO
          }
      }; // End 3D alternative plate_BC class

      class far_BC : public Residual<double>
      {
        public:
          far_BC() : Residual<double> ( 2, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            //TODO  
          }
      }; // End 3D alternative far_BC class
#endif                  
    } // End of namespace Base_Flow
   
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{ 
  cout << "*** ---------- General Ridge Code ---------- ***" << endl;

  /* ----- Solve the base flow ODE ----- */

  //TODO we should be using a non-uniform mesh (for base ODE too)

  // Vector of nodes
  Vector<double> nodes;
  nodes.linspace( 0.0, Param::eta_top, Param::M );

  // Setup the base flow ODE problem
  Base_Flow::equation equation;
  Base_Flow::plate_BC plate_BC;
  Base_Flow::far_BC far_BC;
  ODE_BVP<double> base( &equation, nodes, &plate_BC, &far_BC );

  // Set the initial guess
#ifdef Base_2D
  for (std::size_t j=0; j < Param::N; ++j )
	{
    //TODO maybe think about a better guess involving n
		double eta = nodes[ j ];				                      // eta value at node j
		base.solution()( j, f )  		= eta + exp( -eta );
    base.solution()( j, fd ) 		= 1.0 - exp( -eta ); 
		base.solution()( j, fdd )  	= exp( -eta );
	}
#endif

#ifdef Base_3D

#endif

  // Solve the ODE system and store in a mesh
  base.solve();
  OneD_node_mesh<double> Base_soln( nodes, 6 );
#ifdef Base_2D
  for (std::size_t j=0; j < Param::N; ++j )
	{
		Base_soln( j, UB )      =   base.solution()( j, fd );
    Base_soln( j, UBd )     =   base.solution()( j, fdd );
    Base_soln( j, PhiB )    =   base.solution()( j, f );
    Base_soln( j, ThetaB )  =   ( 1.0 - Param::beta ) * base.solution()( j, fdd );
    Base_soln( j, ThetaBd ) =   ( 1.0 - Param::beta ) * ( - base.solution()( j, f ) *
                                base.solution()( j, fdd ) - Param::beta * ( 1.0 -   
                                base.solution()( j, fd ) * base.solution()( j, fd ) ) );
    Base_soln( j, PsiB )    =   ( 1.0 - Param::beta ) * base.solution()( j, fd );
	}
#endif

#ifdef Base_3D
  
#endif

  Base_soln.output( "./DATA/Base_soln.dat" );             // Output the solution

  cout << "FINISHED" << endl;
}
