/// \file Falkner_Skan_marching.cpp
/// Solving the Falkner-Skan equation with downstream variation
/// \f[ f''' + f f'' + \beta (1 - f'^2 ) = 2xf'f'_x \f]
/// subject to \f$ f(x,0) = s(x) \f$, \f$ f'(x,0) = 0 \f$
/// and \f$ f'(x,\infty) = 1 \f$. Here ' denotes the derivative with
/// respect to eta and s(0)=0.
/// The solution is computed over a range of \f$ x \f$.

#include <cassert>

#include "Core"

// Enumerate the variables
enum { f, fd, fdd };
// Either UNIFORM or NONUNIFORM for uniform of non-uniform mesh
#define NONUNIFORM

namespace TSL
{
  namespace Example
  {
    // Global parameters
    double x_max( 6.28 );               // Maximum downstreaml location
    const std::size_t N_x( 100 );       // Number of steps in the x direction
    double eta_top( 96.0 );             // Size of the domain in the eta direction
    const std::size_t M( 400 );         // Number of intervals in the eta direction
    double beta( 0.0 );                 // Hartree parameter
    double injection( const double& x ) // Injection function
    {
      return sin(x);
    }

    class FS_equation : public Equation_2matrix<double>
    {
    public:

      /// The problem is 3rd order and real
      FS_equation() : Equation_2matrix<double> ( 3 ) {}

      /// Define the equation
      void residual_fn( const Vector<double>& z, Vector<double>& g ) const
      {
        g[ f ]   = z[ fd ];
        g[ fd ]  = z[ fdd ];
        g[ fdd ] = z[ f ] * z[ fdd ] + Example::beta * ( 1.0 - z[ fd ] * z[ fd ]);
      }

      /// Define the (BVP) deriv by providing the identity matrix
      void matrix0( const Vector<double>& z, Matrix<double>& m ) const
      {
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = 1.0;
        m( 2, 2 ) = -1.0;
      }

      /// To speed things up we'll overload this to say the mass matrix is constant
      void get_jacobian_of_matrix0_mult_vector( const Vector<double> &state,
                                    const Vector<double> &vec, Matrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }

      /// Define the unsteady terms by providing the mass matrix
      void matrix1( const Vector<double>& z, Matrix<double>& m ) const
      {
        // x = coord(1)
        m( 2, 1 ) = 2.0 * coord(1) * z[ fd ];
      }

      /// To speed things up we'll overload this to say the mass matrix is constant
      /*void get_jacobian_of_matrix1_mult_vector( const Vector<double> &state,
                                    const Vector<double> &vec, Matrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }*/

    };

    class FS_left_BC : public Residual_with_coords<double>
    {
    public:
      // 2 BCs and 3 unknowns + 1 time coordinate
      FS_left_BC() : Residual_with_coords<double> ( 2, 3, 1 ) {}

      void residual_fn( const Vector<double>& z, Vector<double>& b ) const
      {
        // x = coord(0)
        b[ 0 ] = z[ f ] - Example::injection( coord(0) );
        b[ 1 ] = z[ fd ];
      }
    };

    class FS_right_BC : public Residual_with_coords<double>
    {
    public:
      // 1 BC and 3 unknowns + 1 time coordinate
      FS_right_BC() : Residual_with_coords<double> ( 1, 3, 1 ) {}

      void residual_fn( const Vector<double>& z, Vector<double>& b ) const
      {
        b[ 0 ] = z[ fd ] - 1.0;
      }
    };

  } // end Example namespace

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
    const double b1( 0.3 );
    const double b2( 0.3 );   // Y = (eta_hat + b1)^b2

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
          F[ fdd ] = - u[ f ] * u[ fdd ] - Example::beta * ( 1.0 - u[ fd ] * u[ fd ] );
        }
    }; // End Falkner-Skan equation class

    class plate_BC : public Residual<double>
    {
      public:
        plate_BC() : Residual<double> ( 2, 3 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
          B[ 0 ] = z[ f ];
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

  } // End of namespace Base_Flow

} // end TSL namespace

using namespace TSL;
using namespace std;

int main()
{
  cout << "*** ---------- Falkner-Skan with marching ---------- ***" << endl;

  /* ----- Setup the mesh ----- */

  // define the remapped (non-uniform mesh) domain
  double bottom = Mesh::Y(0.0);
  double top    = Mesh::Y( Example::eta_top );

  // number of points to solve for
  std::size_t N_eta = Example::M + 1;
  std::size_t N_Y( N_eta );

  // nodal positions in the remapped domain (spanned by X,Y)
  Vector<double> Y_nodes;
  Y_nodes.linspace( bottom, top, N_Y );

  // Vectors for original coordinates for writing data on the original zeta-eta domain
  Vector<double> eta_nodes;
  eta_nodes.linspace( 0.0, Example::eta_top, N_eta );

  // to find eta=eta(Y) and zeta=zeta(X) we will use Newton iteration
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
  //

  /* ------ Solve the ODE to get the initial condition at x=0 ------ */
  cout << "*** Solving the ODE to get the initial condition at x = 0 ***" << endl;

  // Setup the base flow ODE problem
  Base_Flow::equation equation;
  Base_Flow::plate_BC plate_BC;
  Base_Flow::far_BC far_BC;
  ODE_BVP<double> base( &equation, eta_nodes, &plate_BC, &far_BC );
  // Set the initial guess
  for (std::size_t j=0; j < N_eta; ++j )
  {
    double eta = eta_nodes[ j ];				                      // eta value at node j
    base.solution()( j, f )  	= eta + exp( -eta );
    base.solution()( j, fd ) 	= 1.0 - exp( -eta );
    base.solution()( j, fdd ) = exp( -eta );
  }
  // Solve the equation
  base.solve_bvp();


  /* ------ Solve the PDE stepping in x ------ */
  cout << "*** Solving the PDE ***" << endl;
  // Falkner-Skan equation
  Example::FS_equation problem;
  // boundary conditions
  Example::FS_left_BC BC_left;
  Example::FS_right_BC BC_right;

  double x( 0.0 );
  // vector of x nodes
  Vector<double> x_nodes;
  x_nodes.linspace( 0, Example::x_max, Example::N_x );
  // x step
  double dx = x_nodes[1] - x_nodes[0];
  // test tolerance
  double tol = 1.e-5;

  // construct our IBVP
  PDE_IBVP<double> Falkner( &problem, eta_nodes, &BC_left, &BC_right );

  // Setup the initial condition
  for ( unsigned i = 0; i < N_eta; ++i )
  {
    //double eta = Falkner.solution().coord( i );
    Falkner.solution()( i, f )   = base.solution()( i, f );
    Falkner.solution()( i, fd )  = base.solution()( i, fd );
    Falkner.solution()( i, fdd ) = base.solution()( i, fdd );
  }
  cout << "x = " << x << ", f''(x,0) = " << Falkner.solution()( 0, fdd ) << endl;

  // maximum difference between numerical and series solution at centre point
  double max_diff( 0.0 );
  // time step
  for ( unsigned i = 1; i < Example::N_x; ++i )
  {
    // take a time step
    try
    {
      Falkner.step2( dx );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    x += dx;
    cout << "x = " << x << ", f''(x,0) = " << Falkner.solution()( 0, fdd ) << endl;
  }

  bool failed = true;
  if ( abs( max_diff ) < tol )
  {
    failed = false;
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
