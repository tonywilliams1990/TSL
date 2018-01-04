/* Rayleigh - Here we define the Rayleigh class which is useful for solving
              the Rayleigh equation.
*/

#ifndef RAYLEIGH_H
#define RAYLEIGH_H


#include "Vector.h"
#include "Matrix.h"
#include "OneD_node_mesh.h"
#include "ODE_BVP.h"
#include "Equation_1matrix.h"
#include "Residual.h"
#include "Eigensystem.h"


namespace TSL
{

  template <typename T>
  class Rayleigh
  {
    /// Class to define the Rayleigh equation
    class Rayleigh_equation : public Equation<std::complex<double>, T>
    {
      public:
        // this is 3rd order for local refinement
        Rayleigh_equation( ) : Equation<std::complex<double> , T>( 3 )
        {}

        void residual_fn( const Vector<std::complex<double> > &z, Vector<std::complex<double> > &g ) const
        {
          T y_pos( this -> coord(0) );
          T U( p_BASEFLOW -> get_interpolated_vars( y_pos )[ 0 ] );
          T Udd( p_BASEFLOW -> get_interpolated_vars( y_pos )[ 1 ] );
          g[ 0 ] = z[ 1 ];
          g[ 1 ] = *p_ALPHA * *p_ALPHA * z[ 0 ] + Udd * z[ 0 ] / ( U - z[ 2 ] );
          g[ 2 ] = 0.0;
        }

        /// The matrix for the BVP coord -- in this case its identity
        void matrix0( const Vector< std::complex<double> >& z, Matrix< std::complex<double> >& m ) const
        {
          m.fill_diag( 1.0 );
        }

        double* p_ALPHA;
        OneD_node_mesh<T, T>* p_BASEFLOW;

    }; // End of class Rayleigh_equation

    /// Class for the left boundary condition
    class Rayleigh_left_BC : public Residual<std::complex<double> >
    {
      public:
        // 2 boundary conditions and 3 unknowns
        Rayleigh_left_BC() : Residual<std::complex<double> > ( 2, 3 )
        { }

        void residual_fn( const Vector<std::complex<double> > &z, Vector<std::complex<double> > &B ) const
        {
          B[ 0 ] = z[ 0 ];
          B[ 1 ] = z[ 1 ] - AMPLITUDE;
        }

        std::complex<double> AMPLITUDE;

    }; // End of class Rayleigh_left_BC

    /// Class for the right boundary condition
    class Rayleigh_right_BC : public Residual<std::complex<double> >
    {
      public:
        // 1 boundary condition and 3 unknowns
        Rayleigh_right_BC() : Residual<std::complex<double> > ( 1, 3 ) {}

        void residual_fn( const Vector<std::complex<double> > &z, Vector<std::complex<double> > &B ) const
        {
          B[ 0 ] = z[ 0 ];
        }
    }; // End of class Rayleigh_right_BC

    private:

      double ALPHA;
      OneD_node_mesh<T, T> BASEFLOW;
      OneD_node_mesh< std::complex<double>, T> EIGENVECTORS;
      Vector<std::complex<double> > EIGENVALUES;

    public:

      /// Constructor
      Rayleigh( OneD_node_mesh<T, T> &base_flow, double &alpha )
      {
        ALPHA = alpha;
        BASEFLOW = base_flow;
        EIGENVALUES = Vector< std::complex<double> >( BASEFLOW.get_nnodes() , 0.0 );
        EIGENVECTORS = OneD_node_mesh< std::complex<double>, T>( BASEFLOW.nodes(), 1);
      }

      /* ----- Methods ----- */

      /// Return a pointer to the eigenvectors meshes
      OneD_node_mesh< std::complex<double>, T>& eigenvectors()
      {
        return EIGENVECTORS;
      }

      /// Return a pointer to the vector of eigenvalues
      Vector< std::complex<double> >& eigenvalues()
      {
        return EIGENVALUES;
      }

      /// Return a pointer to the wavenumber
      double& alpha()
      {
        return ALPHA;
      }

      /// Solve the global eigenvalue problem using 2nd order finite-differences
      void global_evp();

      /// Solve the EVP locally as a nonlinear BVP
      void local_evp( std::size_t i_ev );

      /// Iterate on the wavenumber to drive a selected eigenvalue to neutral
      void iterate_to_neutral( std::size_t i_ev );

  }; // End of class Rayleigh

  template <>
  void Rayleigh<std::complex<double> >::global_evp( )
  {
    unsigned N( BASEFLOW.get_nnodes() );
    Matrix<std::complex<double> > A( N, N, 0.0 );
    Matrix<std::complex<double> > B( N, N, 0.0 );

    // Left BC
    A( 0, 0 ) = 1.0;
    B( 0, 0 ) = 0.0;

    // Interior nodes
    for ( std::size_t i = 1; i < N - 1; ++i )
    {
      std::complex<double>  h( BASEFLOW.coord( i ) - BASEFLOW.coord( i - 1 ) );
      std::complex<double>  k( BASEFLOW.coord( i + 1 ) - BASEFLOW.coord( i ) );
      std::complex<double>  sigma( k / h ); // sigma = 1 => uniform mesh
      std::complex<double>  y = BASEFLOW.coord( i );
      std::complex<double>  U = BASEFLOW.get_interpolated_vars( y )[ 0 ];
      std::complex<double>  Udd = BASEFLOW.get_interpolated_vars( y )[ 1 ];
      std::complex<double>  h2 = 0.5 * h * h * sigma * ( sigma + 1.0 );
      A( i, i ) = U * ( -( sigma + 1.0 ) / h2 - ALPHA * ALPHA ) - Udd;
      A( i, i - 1 ) = sigma * U / h2;
      A( i, i + 1 ) = U / h2;
      //
      B( i, i ) = - ( sigma + 1.0 ) / h2 - ALPHA * ALPHA;
      B( i, i - 1 ) = sigma / h2;
      B( i, i + 1 ) = 1. / h2;
    }

    // Right BC
    A( N - 1, N - 1 ) = 1.0;
    B( N - 1, N - 1 ) = 0.0;

    // Max and min real part of the velocity in base flow
    double U_max( BASEFLOW( 0, 0 ).real() );
    double U_min( BASEFLOW( 0, 0 ).real() );
    for ( unsigned i = 1; i < N; ++i )
    {
      U_max = std::max( U_max, BASEFLOW( i, 0 ).real() );
      U_min = std::min( U_min, BASEFLOW( i, 0 ).real() );
    }

    // Setup and solve the generalised eigenvalue problem
    Eigensystem<std::complex<double> > rayleigh_evp;
    bool compute_eigenvectors = true;
    rayleigh_evp.compute( A, B, compute_eigenvectors );

    EIGENVALUES.resize(0);

    // Return eigenvalues based on Howard's semi-circle theorem
    for ( unsigned i=0; i < N; ++i)
    {
      double c_r, c_i;
      c_r = rayleigh_evp.eigenvalues()[i].real();
      c_i = rayleigh_evp.eigenvalues()[i].imag();
      if ( c_r * c_r + c_i * c_i - ( U_max + U_min ) * c_r + U_max * U_min <= 0 )
      {
        EIGENVALUES.push_back( rayleigh_evp.eigenvalues()[i] );
      }
    }
    // Eigenvectors
    Matrix< std::complex<double> > evec_mat = rayleigh_evp.eigenvector_matrix();
    EIGENVECTORS = OneD_node_mesh<std::complex<double>, std::complex<double> >( BASEFLOW.nodes(), N );
    for ( unsigned evec = 0; evec < EIGENVALUES.size(); ++evec )
    {
      for ( unsigned node = 0; node < N; ++node )
      {
        EIGENVECTORS( node, evec ) = evec_mat( evec, node );
      }
    }

  }

  template <typename T>
  void Rayleigh<T>::local_evp( std::size_t i_ev  )
  {
    // number of nodes in the mesh
    std::size_t N = BASEFLOW.get_nnodes();
    // formulate the Rayleigh equation as a BVP
    Rayleigh_equation Rayleigh_problem;
    // boundary conditions
    Rayleigh_left_BC Rayleigh_left;
    Rayleigh_right_BC Rayleigh_right;
    // set the private member data in the objects
    Rayleigh_problem.p_BASEFLOW = &BASEFLOW;
    Rayleigh_problem.p_ALPHA = &ALPHA;

    // pointer to the equation
    ODE_BVP<std::complex<double>, T>* p_Rayleigh;
    p_Rayleigh = new ODE_BVP<std::complex<double>, T>( &Rayleigh_problem, BASEFLOW.nodes(), &Rayleigh_left, &Rayleigh_right );

    p_Rayleigh -> max_iterations() = 30;

    // set the initial guess using the global_evp solve data
    p_Rayleigh -> solution()( 0, 0 ) = EIGENVECTORS( 0, i_ev );
    p_Rayleigh -> solution()( 0, 1 ) = ( EIGENVECTORS( 1, i_ev ) - EIGENVECTORS( 0, i_ev ) ) / ( BASEFLOW.coord( 1 ) - BASEFLOW.coord( 0 ) );
    p_Rayleigh -> solution()( 0, 2 ) = EIGENVALUES[ i_ev ];
    // set the (arbitrary) amplitude from the global_evp solve data
    Rayleigh_left.AMPLITUDE = p_Rayleigh -> solution()( 0, 1 );
    for ( unsigned i = 1; i < N - 1; ++i )
    {
      p_Rayleigh -> solution()( i, 0 ) = EIGENVECTORS( i, i_ev );
      p_Rayleigh -> solution()( i, 1 ) = ( EIGENVECTORS( i + 1, i_ev ) - EIGENVECTORS( i - 1, i_ev ) ) / ( BASEFLOW.coord( i + 1 ) - BASEFLOW.coord( i - 1 ) );
      p_Rayleigh -> solution()( i, 2 ) = EIGENVALUES[ i_ev ];
    }
    p_Rayleigh -> solution()( N - 1, 0 ) = EIGENVECTORS( N - 1, i_ev );
    p_Rayleigh -> solution()( N - 1, 1 ) = ( EIGENVECTORS( N - 1, i_ev ) - EIGENVECTORS( N - 2, i_ev ) ) / ( BASEFLOW.coord( N - 1 ) - BASEFLOW.coord( N - 2 ) );
    p_Rayleigh -> solution()( N - 1, 2 ) = EIGENVALUES[ i_ev ];

    // do a local solve
    p_Rayleigh -> solve_bvp();
    // write the eigenvalue and eigenvector back to private member data store
    EIGENVALUES[ i_ev ] = p_Rayleigh -> solution()( 0, 2 );
    for ( unsigned i = 0; i < N; ++i )
    {
      EIGENVECTORS( i, i_ev ) = p_Rayleigh -> solution()( i, 0 );
    }
    // delete the equation object
    delete p_Rayleigh;
  }

  template<>
  void Rayleigh<double>::iterate_to_neutral( std::size_t i_ev ) {}

  template<>
  void Rayleigh<std::complex<double> >::iterate_to_neutral( std::size_t i_ev )
  {
    double delta( 1.e-8 );
    do
    {
      std::cout << "ALPHA = " << ALPHA << "\n";
      local_evp( i_ev );
      std::complex<double> copy_of_ev( EIGENVALUES[ i_ev ] );
      ALPHA += delta;
      local_evp( i_ev );
      ALPHA -= delta;
      double d_ev = ( std::imag( EIGENVALUES[ i_ev ] ) - std::imag( copy_of_ev ) ) / delta;
      ALPHA -= std::imag( copy_of_ev ) / d_ev;
      std::cout << "ITERATING: " << ALPHA << " " << EIGENVALUES[ i_ev ] << " " << d_ev << "\n";
    }
    while ( std::abs( std::imag( EIGENVALUES[ i_ev ] ) ) > 1.e-6 );
  }

} // End of namespace TSL

#endif
