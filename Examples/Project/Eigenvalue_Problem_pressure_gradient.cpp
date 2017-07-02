#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <sstream>

#include "Core"
#include "Eigenvalue"


// ODE enumeration
enum{ f, fd, fdd, g, gd, gdd };

// Base flow enumeration
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };

// Eigenvalue enumeration
enum{ u, ud, phi, theta, thetad, psi };

// Either Base_2D or Base_3D for 2D or 3D base flows
#define Base_2D
#define UNIFORM

//TODO K iteration is temperamental -> just use single values of K

namespace TSL
{
    namespace Param
    {
      double eta_top( 96.0 );             // Size of the domain in the eta direction
      const std::size_t M( 203 );         // Number of intervals in the eta direction
      double beta( 0.5 );                 // Hartree parameter
      double K( 9.0 );                    // Transpiration parameter ( +ve = blowing )
      double K_max( 3.0 );                // Maximum value of K
      double K_min( -2.0 );               // Minimum value of K
      double dK( 0.1 );                   // K increment
      bool compute_eigenvectors( true );  // Compute the eigenvectors or not

    } // End of namespace Param

    namespace Example
    {
      std::string output_path;          // Output path
    } // End of namespace Example

    unsigned col( const unsigned& i, const unsigned& k )
    {
        // Return the column number for kth variable at node i
        return 6 * i  + k;
    }

    namespace Mesh
    {
#ifdef UNIFORM
      double Y( const double& eta )
      {
        return eta;
      }
      double Yd( const double& eta )
      {
        return 1;
      }
      double Ydd( const double& eta )
      {
        return 0;
      }
#endif
#ifdef NONUNIFORM
      const double b1( 0.5 );
      const double b2( 0.5 );   // Y = (eta_hat + b1)^b2

      double Y( const double& eta )
      {
        return std::pow(eta + b1, b2) - pow(b1,b2);
      }
      double Yd( const double& eta )
      {
        return b2 * std::pow(eta + b1, b2 - 1);
      }
      double Ydd( const double& eta )
      {
        return b2 * (b2 - 1) * std::pow(eta + b1, b2 - 2);
      }
#endif

      class invert_eta : public Residual<double>
      {
        // USED TO INVERT THE NON-UNIFORM MESH MAPPING
        public:
        double Y0;

        invert_eta() : Residual<double>( 1 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &f ) const
        {
          f[ 0 ] = Y( z[0] ) - Y0;
        }
      };

    } // End of namespace Mesh

    namespace Base_Flow
    {
#ifdef Base_2D
      class equation : public Equation<double>
      {
        public:
          double beta;                            // Hartree parameter
          // Falkner-Skan equation is 3rd order
          equation() : Equation<double> ( 3 ) {}
          // Define the equation
          void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
          {
            F[ f ]   = u[ fd ];
            F[ fd ]  = u[ fdd ];
            F[ fdd ] = - u[ f ] * u[ fdd ] - beta * ( 1.0 - u[ fd ] * u[ fd ] );
          }
      }; // End Falkner-Skan equation class

      class plate_BC : public Residual<double>
      {
        public:
          double K;                              // Transpiration parameter

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
          double beta;                     // Hartree parameter
          // The 3D alternative equation is 6th order
          equation() : Equation<double> ( 6 ) {}
          // Define the equation
          void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
          {
            F[ f ]    =  u[ fd ];
            F[ fd ]   =  u[ fdd ];
            F[ fdd ]  = -( u[ f ] + ( 2.0 - beta ) * u[ g ] ) * u[ fdd ]
                        - beta * ( 1.0 - u[ fd ] * u[ fd ] );
            F[ g ]    =  u[ gd ];
            F[ gd ]   =  u[ gdd ];
            F[ gdd ]  = -( u[ f ] + ( 2.0 - beta ) * u[ g ] ) * u[ gdd ]
                        -( 2.0 * ( 1.0 - beta ) * u[ fd ]
                        - ( 2.0 - beta) * u[ gd ] ) * u[ gd ];
          }
      }; // End 3D alternative equation class

      class plate_BC : public Residual<double>
      {
        public:
          double K;                        // Transpiration parameter

          plate_BC() : Residual<double> ( 4, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + K;
            B[ 1 ] = z[ fd ];
            B[ 2 ] = z[ g ];
            B[ 3 ] = z[ gd ];
          }
      }; // End 3D alternative plate_BC class

      class far_BC : public Residual<double>
      {
        public:
          far_BC() : Residual<double> ( 2, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ fd ] - 1.0;
            B[ 1 ] = z[ gd ];
          }
      }; // End 3D alternative far_BC class
#endif
    } // End of namespace Base_Flow
} // End of namespace TSL

using namespace std;
using namespace TSL;

int main()
{
  cout << "----- Eigenvalue problem -----" << endl;

  /* ----- Make the output directory ----- */
  std::ostringstream ss;
  ss << "./DATA/Eigenvalue_problem_beta_" << Param::beta << "/";
  Example::output_path = ss.str();
  int status = mkdir( Example::output_path.c_str(), S_IRWXU );
  if ( status == 0 ) {
  cout << "  * Output directory " + Example::output_path +
          " has been made successfully." << endl;
  }

  /* ----- Setup the mesh ----- */

  // define the remapped (non-uniform mesh) domain
  double bottom = Mesh::Y(0.0);
  double top    = Mesh::Y( Param::eta_top );

  // number of points to solve for
  std::size_t N_eta = Param::M + 1;
  std::size_t N_Y( N_eta );

  // nodal positions in the remapped domain (spanned by X,Y)
  Vector<double> Y_nodes;
  Y_nodes.linspace( bottom, top, N_Y );

  // Vectors for original coordinates for writing data on the original zeta-eta domain
  Vector<double> eta_nodes;
  eta_nodes.linspace( 0.0, Param::eta_top, N_eta );

  // to find eta=eta(Y) we will use Newton iteration
  Mesh::invert_eta find_eta;
  Newton<double> newton_eta( &find_eta );
  for ( unsigned j = 0; j < N_Y; ++j )
  {
    unsigned kmin(0); double min(99);
    for ( unsigned k = 0; k < N_Y; ++k )
    {
      if ( std::abs( Mesh::Y( eta_nodes[k] ) - Y_nodes[j] ) < min )
      {
        min = std::abs( Mesh::Y( eta_nodes[k] ) - Y_nodes[j] );
        kmin = k;
      }
    }
    find_eta.Y0 = Y_nodes[ j ];
    Vector<double> guess( 1, 1.0 );
    guess[ 0 ] = eta_nodes[ kmin ];
    newton_eta.iterate( guess );
    eta_nodes[j] = guess[ 0 ];
  }

  // step size in the remapped domain
  const double dY( Y_nodes[ 1 ] - Y_nodes[ 0 ] );
  //const double delta = Param::eta_top / ( N_eta - 1 ); // Mesh grid spacing

  cout << "*** Solving the base flow ODE using " << N_eta << " nodes." << endl;

	// Create instances of the equation and BCs
  Base_Flow::equation equation;
  Base_Flow::plate_BC plate_BC;
  Base_Flow::far_BC far_BC;

  ODE_BVP<double> bvp( &equation, eta_nodes, &plate_BC, &far_BC );

  // Set the initial guess
#ifdef Base_2D
	for (std::size_t j=0; j < N_eta; ++j )
	{
		double eta = eta_nodes[ j ];				                      // eta value at node j
		bvp.solution()( j, f )  		= eta + exp( -eta );
    bvp.solution()( j, fd ) 		= 1.0 - exp( -eta );
		bvp.solution()( j, fdd )  	= exp( -eta );
	}
#endif
#ifdef Base_3D
    for (std::size_t j=0; j < N_eta; ++j )
	{
		double eta = eta_nodes[ j ];					                   // eta value at node j
		bvp.solution()( j, f )  		= eta + exp( -eta );
    bvp.solution()( j, fd ) 		= 1.0 - exp( -eta );
		bvp.solution()( j, fdd )  	= exp( -eta );
    bvp.solution()( j, g )  		= 0.35 * (1.0 - exp( -eta ));
    bvp.solution()( j, gd ) 		= 1 - exp( -eta ) - exp( -1 / (eta * eta) );
		bvp.solution()( j, gdd )  	= exp( -eta ) - 0.5 * tanh( eta ) + 0.5 * tanh( eta - 2.0 );
	}
#endif
// Loop over values of K

do{
    // Set the parameters
    equation.beta = Param::beta;
    plate_BC.K = Param::K;
    // Solve the BVP
    bvp.solve_bvp();

    cout << "*** K = " << Param::K << endl;

    // Store the solution in a mesh
    OneD_node_mesh<double> Base_soln( eta_nodes, 6 );
#ifdef Base_2D
  for (std::size_t j=0; j < N_eta; ++j )
	{
		Base_soln( j, UB )      =   bvp.solution()( j, fd );
    Base_soln( j, UBd )     =   bvp.solution()( j, fdd );
    Base_soln( j, PhiB )    =   bvp.solution()( j, f );
    Base_soln( j, ThetaB )  =   ( 1.0 - Param::beta ) * bvp.solution()( j, fdd );
    Base_soln( j, ThetaBd ) =   ( 1.0 - Param::beta ) * ( - bvp.solution()( j, f ) *
                                bvp.solution()( j, fdd ) - Param::beta * ( 1.0 -
                                bvp.solution()( j, fd ) * bvp.solution()( j, fd ) ) );
    Base_soln( j, PsiB )    =   ( 1.0 - Param::beta ) * bvp.solution()( j, fd );
	}
#endif
#ifdef Base_3D
  for (std::size_t j=0; j < N_eta; ++j )
	{
		Base_soln( j, UB )      =   bvp.solution()( j, fd );
    Base_soln( j, UBd )     =   bvp.solution()( j, fdd );
    Base_soln( j, PhiB )    =   bvp.solution()( j, f )
                              + ( 2.0 - Param::beta ) * bvp.solution()( j, g );
    Base_soln( j, ThetaB )  =   ( 1.0 - Param::beta ) * bvp.solution()( j, fdd )
                              - ( 2.0 - Param::beta ) * bvp.solution()( j, gdd );
    Base_soln( j, ThetaBd ) =   ( 1.0 - Param::beta ) * ( -(bvp.solution()( j, f ) +
                                (2.0 - Param::beta) * bvp.solution()( j, g )) *
                                bvp.solution()( j, fdd ) - Param::beta * ( 1.0 -
                                bvp.solution()( j, fd ) * bvp.solution()( j, fd ) ) )
                              - ( 2.0 - Param::beta ) * ( -(bvp.solution()( j, f ) +
                                (2.0 - Param::beta) * bvp.solution()( j, g )) *
                                bvp.solution()( j, gdd ) - Param::beta * ( 1.0 -
                                bvp.solution()( j, gd ) * bvp.solution()( j, gd ) ) -
                                2.0 * (1.0 - Param::beta ) * (bvp.solution()( j, fd ) -
                                bvp.solution()( j, gd )) * bvp.solution()( j, gd ) );
    Base_soln( j, PsiB )    =   ( 1.0 - Param::beta ) * bvp.solution()( j, fd )
                              - ( 2.0 - Param::beta ) * bvp.solution()( j, gd );
	}
#endif

    Base_soln.output( Example::output_path + "Base_soln_K_"
                    + Utility::stringify( Param::K, 1, "fixed" )+".dat" );  // Output the solution

    cout << "*** Assembling the matrices for the eigenvalue problem." << endl;

    // Create the generalised eigenvalue problem A v = lambda B v
    Matrix<double> A( 6*N_eta, 6*N_eta, 0.0 ); // 6N*6N -> 6th order system
    Matrix<double> B( 6*N_eta, 6*N_eta, 0.0 );

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
    for ( std::size_t i=0; i<N_eta-1; ++i)
    {
        // eta at mid-node location i + 1/2
        double eta( ( eta_nodes[ i + 1 ] + eta_nodes[ i ] ) / 2 );
        double Yd( Mesh::Yd( eta ) );

        // Base solution at the mid-node location i + 1/2
        double U_B      =   0.5 * ( Base_soln( i, UB ) + Base_soln( i + 1, UB )  );
        double U_Bd     =   0.5 * ( Base_soln( i, UBd ) + Base_soln( i + 1, UBd )  );
        double Phi_B    =   0.5 * ( Base_soln( i, PhiB ) + Base_soln( i + 1, PhiB )  );
        double Theta_B  =   0.5 * ( Base_soln( i, ThetaB ) + Base_soln( i + 1, ThetaB )  );
        double Theta_Bd =   0.5 * ( Base_soln( i, ThetaBd ) + Base_soln( i + 1, ThetaBd )  );
        double Psi_B    =   0.5 * ( Base_soln( i, PsiB ) + Base_soln( i + 1, PsiB )  );

        // Equation 1
        A( row, col( i, u ) )           =  ( 2.0 - Param::beta ) / 2.0;
        A( row, col( i + 1, u ) )       =  ( 2.0 - Param::beta ) / 2.0;
        A( row, col( i, phi ) )         =  Yd / dY;
        A( row, col( i + 1, phi ) )     = -Yd / dY;
        A( row, col( i, psi ) )         = -0.5;
        A( row, col( i + 1, psi ) )     = -0.5;
        B( row, col( i, psi ) )         =  0.5;
        B( row, col( i + 1, psi ) )     =  0.5;
        ++row;
        // Equation 2
        A( row, col( i, theta ) )       =  0.5;
        A( row, col( i + 1, theta ) )   =  0.5;
        A( row, col( i, psi ) )         =  Yd / dY;
        A( row, col( i + 1, psi ) )     = -Yd / dY;
        ++row;
        // Equation 3
        A( row, col( i, ud ) )          =  0.5;
        A( row, col( i + 1, ud ) )      =  0.5;
        A( row, col( i, u ) )           =  Yd / dY;
        A( row, col( i + 1, u ) )       = -Yd / dY;
        ++row;
        // Equation 4
        A( row, col( i, ud ) )          =  0.5 * Phi_B - Yd / dY;
        A( row, col( i + 1, ud ) )      =  0.5 * Phi_B + Yd / dY;
        A( row, col( i, u ) )           = -2.0 * Param::beta * U_B / 2.0;
        A( row, col( i + 1, u ) )       = -2.0 * Param::beta * U_B / 2.0;
        A( row, col( i, phi ) )         =  0.5 * U_Bd;
        A( row, col( i + 1, phi ) )     =  0.5 * U_Bd;
        B( row, col( i, u ) )           = -0.5 * Psi_B;
        B( row, col( i + 1, u ) )       = -0.5 * Psi_B;
        ++row;
        // Equation 5
        A( row, col( i, theta ) )       = -Yd / dY;
        A( row, col( i + 1, theta ) )   =  Yd / dY;
        A( row, col( i, thetad ) )      = -0.5;
        A( row, col( i + 1, thetad ) )  = -0.5;
        ++row;
        // Equation 6
        A( row, col( i, u ) )           =  ( ( 2.0 - Param::beta ) * Theta_B / 2.0 )
                                          -( 2.0 * ( 1.0 - Param::beta ) * U_Bd / 2.0 );
        A( row, col( i + 1, u ) )       =  ( ( 2.0 - Param::beta ) * Theta_B / 2.0 )
                                          -( 2.0 * ( 1.0 - Param::beta ) * U_Bd / 2.0 );
        A( row, col( i, ud ) )          = -( 2.0 * ( 1.0 - Param::beta ) * U_B / 2.0 );
        A( row, col( i + 1, ud ) )      = -( 2.0 * ( 1.0 - Param::beta ) * U_B / 2.0 );
        A( row, col( i, phi ) )         =  0.5 * Theta_Bd;
        A( row, col( i + 1, phi ) )     =  0.5 * Theta_Bd;
        A( row, col( i, theta ) )       =  0.5 * Psi_B + ( 2.0 - Param::beta ) * U_B * 0.5;
        A( row, col( i + 1, theta ) )   =  0.5 * Psi_B + ( 2.0 - Param::beta ) * U_B * 0.5;
        A( row, col( i, thetad ) )      =  0.5 * Phi_B - Yd / dY;
        A( row, col( i + 1, thetad ) )  =  0.5 * Phi_B + Yd / dY;
        A( row, col( i, psi ) )         =  0.5 * Theta_B;
        A( row, col( i + 1, psi ) )     =  0.5 * Theta_B;
        B( row, col( i, theta ) )       = -0.5 * Psi_B;
        B( row, col( i + 1, theta ) )   = -0.5 * Psi_B;
        ++row;
    }

    // Far BCs
    i = N_eta - 1;
    // u(inf) = 0
    A( row, col( i, u ) )               =  1.0;
    ++row;
    // theta(inf) = 0
    A( row, col( i, theta ) )           =  1.0;
    ++row;
    // psi(0) = 0
    A( row, col( i, psi ) )             =  1.0;
    //++row;

    cout << "*** Solving the generalised eigenvalue problem A v=lambda B v." << endl;

    // Compute the eigenvalues ( and possibly the eigenvectors )
    Eigensystem<double> system;
    if ( Param::compute_eigenvectors ){ cout << "*** Computing the eigenvectors." << endl;}
    Timer timer;
	  timer.start();
    system.compute( A, B, Param::compute_eigenvectors );

    if ( system.eigenvectors_computed() )
    {
      cout << "*** The Eigenvectors have been computed. " << endl;
    }
    // Get eigenvalues and eigenvectors
    Vector< std::complex<double> > evals = system.eigenvalues();
    std::vector< Vector< std::complex<double> > > evecs = system.eigenvectors();

    // Find the real parts of the eigenvalues
    Vector<double> real_evals; // Vector to store the real part of relevant eigenvalues
    for (size_t i=0; i < evals.size(); ++i)
    {
        if ( evals[i].real() < 3.0 && evals[i].real() > -1.0 && abs( evals[i] ) < 50.0 )
        {
            cout << evals[i] << endl;
            real_evals.push_back( evals[i].real() );
        }
    }

    if ( Param::compute_eigenvectors )
    {

      // Find the eigenvector that corresponds to eigenvalue with negative real part
      Vector< std::complex<double> > negative_evec;
      for (size_t i=0; i < evals.size(); ++i)
      {
          if ( evals[i].real() < 0.0 && evals[i].real() > -1.0 && abs( evals[i] ) < 50.0 )
          {
              negative_evec = evecs[i];
          }
      }

      // Separate the eigenvectors into [u,ud,phi,theta,thetad,psi] and put into a mesh
      OneD_node_mesh< double > evec_mesh( eta_nodes, 6 );
      for ( std::size_t i=0; i<eta_nodes.size(); ++i )
      {
        for ( std::size_t j=0; j<6; ++j )
        {
          // We only want the real part
          evec_mesh( i, j ) = negative_evec[ i * 6 + j ].real();
        }
      }
      evec_mesh.output( Example::output_path + "Eigenvectors_K_"
                      + Utility::stringify( Param::K, 1, "fixed" )+ ".dat", 10);
    }

    timer.print();
	  timer.stop();

    real_evals.output( Example::output_path + "Eigenvalues_K_"
                     + Utility::stringify( Param::K, 1, "fixed" )+ ".dat", 10);

    Param::K += Param::dK;

}while( Param::K <= Param::K_max + 0.0001 && Param::K > Param::K_min );

    cout << "FINISHED" << endl;
}
