#include <cassert>
#include <cmath>

#include "Newton.h"
#include "Residual.h"
#include "ODE_BVP.h"
#include "Vector.h"
#include "Eigensystem.h"
#include "Timer.h"


// ODE enumeration
enum{ f, fd, fdd, g, gd, gdd };

// Base flow enumeration
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };

// Eigenvalue enumeration
enum{ u, ud, phi, theta, thetad, psi };

// Either Base_2D or Base_3D for 2D or 3D base flows
#define Base_3D

namespace TSL
{
    unsigned col( const unsigned& i, const unsigned& k )
    {
        // Return the column number for kth variable at node i
        return 6 * i  + k;
    }

    namespace Base_Flow
    {
        double K( 0.0 );   // Transpiration parameter (+ve is blowing)

#ifdef Base_2D
	    class equation : public Equation<double>
	    {
		    public:
			    // The Blasius equation is 3rd order
			    equation() : Equation<double> ( 3 ) {}

			    // Define the equation
			    void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			    {
				    F[ f ]   = u[ fd ];
				    F[ fd ]  = u[ fdd ];
				    F[ fdd ] = - u[ f ] * u[ fdd ];
			    }
	    };

        class plate_BC : public Residual<double>
        {
            public:
                plate_BC() : Residual<double> ( 2, 3 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ f ] + K;
                B[ 1 ] = z[ fd ];
            }
        };

        class far_BC : public Residual<double>
        {
            public:
                far_BC() : Residual<double> ( 1, 3 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ fd ] - 1.0;
            }
        };
#endif

#ifdef Base_3D
        class equation : public Equation<double>
	    {
		    public:
			    // The equation is 6th order
			    equation() : Equation<double> ( 6 ) {}

			    // Define the equation
			    void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			    {
				    F[ f ]   = u[ fd ];
				    F[ fd ]  = u[ fdd ];
				    F[ fdd ] = - ( u[ f ] +  u[ g ] ) * u[ fdd ];
            F[ g ]   = u[ gd ];
            F[ gd ]  = u[ gdd ];
            F[ gdd ] = - ( u[ f ] + u[ g ] ) * u[ gdd ] - 2.0 * u[ gd ] * u[ fd ] + u[ gd ] * u[ gd ];
			    }
	    };

        class plate_BC : public Residual<double>
        {
            public:
                plate_BC() : Residual<double> ( 4, 6 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ f ] + K;
                B[ 1 ] = z[ fd ];
                B[ 2 ] = z[ g ];
                B[ 3 ] = z[ gd ];
            }
        };

        class far_BC : public Residual<double>
        {
            public:
                far_BC() : Residual<double> ( 2, 6 ) {}

            void residual_fn( const Vector<double> &z, Vector<double> &B ) const
            {
                B[ 0 ] = z[ fd ] - 1.0;
                B[ 1 ] = z[ gd ];
            }
        };
#endif

    } // End of namespace Base_Flow
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{
    cout << "----- Eigenvalue problem -----" << endl;

	// Define the domain
	double Inf( 30.0 );											           // Infinite boundary
	size_t N_nodes( 200 );                             // Number of nodes
	Vector<double> nodes;								               // Declare vector of nodes (uniform)
	nodes.linspace(0,Inf,N_nodes);
  const double delta = Inf / ( N_nodes - 1 );        // Mesh grid spacing

  cout << "*** Solving the base flow ODE using " << N_nodes << " points." << endl;

	// Create instances of the equation and BCs
  Base_Flow::equation equation;
  Base_Flow::plate_BC left_BC;
  Base_Flow::far_BC right_BC;

  ODE_BVP<double> bvp( &equation, nodes, &left_BC, &right_BC ); // Create boundary value problem

    // Set the initial guess
#ifdef Base_2D
	for (std::size_t j=0; j < N_nodes; ++j )
	{
		double eta = nodes[ j ];				                      // eta value at node j
		bvp.solution()( j, f )  		= eta + exp( -eta );
    bvp.solution()( j, fd ) 		= 1.0 - exp( -eta );
		bvp.solution()( j, fdd )  	= exp( -eta );
	}
#endif
#ifdef Base_3D
    for (std::size_t j=0; j < N_nodes; ++j )
	{
		double eta = nodes[ j ];					                   // eta value at node j
		bvp.solution()( j, f )  		= eta + exp( -eta );
    bvp.solution()( j, fd ) 		= 1.0 - exp( -eta );
		bvp.solution()( j, fdd )  	= exp( -eta );
    bvp.solution()( j, g )  		= 0.35 * (1.0 - exp( -eta ));
    bvp.solution()( j, gd ) 		= 1 - exp( -eta ) - exp( -1 / (eta * eta) );
		bvp.solution()( j, gdd )  	= exp( -eta ) - 0.5 * tanh( eta ) + 0.5 * tanh( eta - 2.0 );
	}
#endif

    // Solve with no transpiration for a good initial guess
    Base_Flow::K = 0.0;
    bvp.solve_bvp();

    // Set the the transpiration value and solve again
    Base_Flow::K = 0.0;
    bvp.solve_bvp();

    cout << "*** K = " << Base_Flow::K << endl;

    // Create a mesh for storing the base flow solution
    OneD_node_mesh<double> Base_soln( nodes, 6 );

    // Populate the mesh with the base flow solution
#ifdef Base_2D
    for (std::size_t j=0; j < N_nodes; ++j )
	{
		Base_soln( j, UB )      =   bvp.solution()( j, fd );
    Base_soln( j, UBd )     =   bvp.solution()( j, fdd );
    Base_soln( j, PhiB )    =   bvp.solution()( j, f );
    Base_soln( j, ThetaB )  =   bvp.solution()( j, fdd );
    Base_soln( j, ThetaBd ) = - bvp.solution()( j, f ) * bvp.solution()( j, fdd );
    Base_soln( j, PsiB )    =   bvp.solution()( j, fd );
	}
#endif
#ifdef Base_3D
    for (std::size_t j=0; j < N_nodes; ++j )
	{
		Base_soln( j, UB )      =   bvp.solution()( j, fd );
    Base_soln( j, UBd )     =   bvp.solution()( j, fdd );
    Base_soln( j, PhiB )    =   bvp.solution()( j, f ) + bvp.solution()( j, g );
    Base_soln( j, ThetaB )  =   bvp.solution()( j, fdd ) - bvp.solution()( j, gdd );
    Base_soln( j, ThetaBd ) =   ( bvp.solution()( j, f ) + bvp.solution()( j, g ) ) * ( bvp.solution()( j, gdd ) -
                                    bvp.solution()( j, fdd ) ) + 2.0 * bvp.solution()( j, gd ) * bvp.solution()( j, fd )
                                    - bvp.solution()( j, gd ) * bvp.solution()( j, gd ) ;
    Base_soln( j, PsiB )    =   bvp.solution()( j, fd ) - bvp.solution()( j, gd );
	}
#endif

    Base_soln.output( "./DATA/Base_soln.dat" );              // Output the solution

    cout << "*** Assembling the matrices for the eigenvalue problem." << endl;

    // Create the generalised eigenvalue problem A v = lambda B v
    Matrix<double> A( 6*N_nodes, 6*N_nodes, 0.0 ); // 6N*6N -> 6th order system
    Matrix<double> B( 6*N_nodes, 6*N_nodes, 0.0 );

    unsigned row( 0 );                     // Row counter

    // Plate BCs
    unsigned i( 0 );
    // u(0) = 0
    A( row, col( i, u ) )               =  1.0;
    ++row;
    // phi(0) = 0
    A( row, col( i, phi ) )             =  1.0;
    ++row;
    // psi(0) = 0
    A( row, col( i, psi ) )             =  1.0;
    ++row;

    // Interior points
    for ( std::size_t i=0; i<N_nodes-1; ++i)
    {
        // Base solution at the mid-node location i + 1/2
        double U_B      = 0.5 * ( Base_soln( i, UB      ) + Base_soln( i + 1, UB      )  );
        double U_Bd     = 0.5 * ( Base_soln( i, UBd     ) + Base_soln( i + 1, UBd     )  );
        double Phi_B    = 0.5 * ( Base_soln( i, PhiB    ) + Base_soln( i + 1, PhiB    )  );
        double Theta_B  = 0.5 * ( Base_soln( i, ThetaB  ) + Base_soln( i + 1, ThetaB  )  );
        double Theta_Bd = 0.5 * ( Base_soln( i, ThetaBd ) + Base_soln( i + 1, ThetaBd )  );
        double Psi_B    = 0.5 * ( Base_soln( i, PsiB    ) + Base_soln( i + 1, PsiB    )  );
        // Equation 1
        A( row, col( i, u ) )           =  1.0;
        A( row, col( i + 1, u ) )       =  1.0;
        A( row, col( i, phi ) )         =  1.0/delta;
        A( row, col( i + 1, phi ) )     = -1.0/delta;
        A( row, col( i, psi ) )         = -0.5;
        A( row, col( i + 1, psi ) )     = -0.5;
        B( row, col( i, psi ) )         =  0.5;
        B( row, col( i + 1, psi ) )     =  0.5;
        ++row;
        // Equation 2
        A( row, col( i, theta ) )       =  0.5;
        A( row, col( i + 1, theta ) )   =  0.5;
        A( row, col( i, psi ) )         =  1.0/delta;
        A( row, col( i + 1, psi ) )     = -1.0/delta;
        ++row;
        // Equation 3
        A( row, col( i, ud ) )          =  0.5;
        A( row, col( i + 1, ud ) )      =  0.5;
        A( row, col( i, u ) )           =  1.0/delta;
        A( row, col( i + 1, u ) )       = -1.0/delta;
        ++row;
        // Equation 4
        A( row, col( i, ud ) )          =  0.5 * Phi_B - 1.0/delta;
        A( row, col( i + 1, ud ) )      =  0.5 * Phi_B + 1.0/delta;
        A( row, col( i, phi ) )         =  0.5 * U_Bd;
        A( row, col( i + 1, phi ) )     =  0.5 * U_Bd;
        B( row, col( i, u ) )           = -0.5 * Psi_B;
        B( row, col( i + 1, u ) )       = -0.5 * Psi_B;
        ++row;
        // Equation 5
        A( row, col( i, theta ) )       =  1.0/delta;
        A( row, col( i + 1, theta ) )   = -1.0/delta;
        A( row, col( i, thetad ) )      =  0.5;
        A( row, col( i + 1, thetad ) )  =  0.5;
        ++row;
        // Equation 6
        A( row, col( i, u ) )           =  Theta_B - U_Bd;
        A( row, col( i + 1, u ) )       =  Theta_B - U_Bd;
        A( row, col( i, ud ) )          = -U_B;
        A( row, col( i + 1, ud ) )      = -U_B;
        A( row, col( i, phi ) )         =  0.5 * Theta_Bd;
        A( row, col( i + 1, phi ) )     =  0.5 * Theta_Bd;
        A( row, col( i, theta ) )       =  0.5 * Psi_B + U_B;
        A( row, col( i + 1, theta ) )   =  0.5 * Psi_B + U_B;
        A( row, col( i, thetad ) )      =  0.5 * Phi_B - 1.0/delta;
        A( row, col( i + 1, thetad ) )  =  0.5 * Phi_B + 1.0/delta;
        A( row, col( i, psi ) )         =  0.5 * Theta_B;
        A( row, col( i + 1, psi ) )     =  0.5 * Theta_B;
        B( row, col( i, theta ) )       = -0.5 * Psi_B;
        B( row, col( i + 1, theta ) )   = -0.5 * Psi_B;
        ++row;
    }

    // Far BCs
    i = N_nodes - 1;
    // u(inf) = 0
    A( row, col( i, u ) )               =  1.0;
    ++row;
    // theta(inf) = 0
    A( row, col( i, theta ) )           =  1.0;
    ++row;
    // psi(inf) = 0
    A( row, col( i, psi ) )             =  1.0;
    //++row;

    cout << "*** Solving the generalised eigenvalue problem A v=lambda B v." << endl;

    // Compute the eigenvalues
    Eigensystem<double> system;                        // Construct Teigensystem object
    bool compute_eigenvectors = false;
    TSL::Timer timer;								// Create a Timer object
	  timer.start();								  // Start the timer
    system.compute( A, B, compute_eigenvectors );

    TSL::Vector< std::complex<double> > evals = system.eigenvalues();

    for (size_t i=0; i < evals.size(); ++i)
    {
        if ( evals[i].real() < 2.0 && evals[i].real() > -1.0 && abs( evals[i] ) < 50.0 )
        {
            cout << evals[i] << endl;
        }
    }

    timer.print();                                      // Output time to screen
	  timer.stop();                                       // Stop the timer

    cout << "FINISHED" << endl;
}
