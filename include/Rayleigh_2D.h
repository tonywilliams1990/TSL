/* Rayleigh_2D - Here we define the Rayleigh_2D class which is useful for solving
              the Rayleigh equation.
*/

#ifndef RAYLEIGH_2D_H
#define RAYLEIGH_2D_H

#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <sstream>
#include <string>

#include "Vector.h"
#include "SparseMatrix.h"
#include "SelfSimInjection.h"
#include "SparseEigenSystem.h"
#include "TwoD_node_mesh.h"
#include "Error.h"

namespace TSL
{
  class Rayleigh_2D
  {
    protected:
      SelfSimInjection SSI;                                   // Self similar injection object
      Vector< std::complex<double> > EIGENVALUES;             // Eigenvalues
      TwoD_node_mesh< std::complex<double> > EIGENVECTORS;    // Eigenvectors
      double ALPHA;                                           // Wavenumber
      std::size_t EIGS_REQUESTED;                             // Number of eigenvalues requested
      double LEFT;                                            // Region in which to look for
      double RIGHT;                                           // eigenvalues
      double BOTTOM;
      double TOP;
      std::complex<double> TARGET;                            // Target eigenvalue
      std::string ORDER;                                      // Ordering of returned eigenvalues
      bool CALC_EIGENVECTORS;                                 // Calculate the eigenvectors?
      std::string OUTPUT_PATH;

    public:

      /// Constructor
      Rayleigh_2D( SelfSimInjection& ssi, double& alpha, std::size_t& nev )
      {
        SSI = ssi;
        ALPHA = alpha;
        EIGS_REQUESTED = nev;
        //TODO give default values to other members
        CALC_EIGENVECTORS = true;
        OUTPUT_PATH = SSI.output_path();
      }

      /// Destructor
	   	~Rayleigh_2D()
      {
        SlepcFinalize();
      }

      /* ----- Methods ----- */

      /// Return a pointer to the eigenvectors
      TwoD_node_mesh< std::complex<double> >& eigenvectors()
      {
        return EIGENVECTORS;
      }

      /// Return a pointer to the vector of eigenvalues
      Vector< std::complex<double> >& eigenvalues()
      {
        return EIGENVALUES;
      }

      /// Return a pointer to the wavenumber
      double& wavenumber()
      {
        return ALPHA;
      }

      /// Set the region of the complex plane in which to look for eigenvalues
      void set_region( const double& left, const double& right, const double& bottom, const double& top )
      {
        LEFT   = left;
        RIGHT  = right;
        BOTTOM = bottom;
        TOP    = top;
      }

      /// Set target eigenvalue
      void set_target( const std::complex<double>& target )
      {
        TARGET = target;
      }

      /// Set order
      void set_order( const std::string& order )
      {
        //TODO probably should do some checking to make sure it is a valid string
        ORDER = order;
      }

      /// Calculate eigenvectors (return a handle to CALC_EIGENVECTORS)
      bool& calc_eigenvectors()
      {
        return CALC_EIGENVECTORS;
      }

      /// Solve the sparse eigenvalue problem using SLEPc
      void solve_evp();

      /// Solve the eigenvalue problem for whilst stepping in ALPHA
      void step_in_alpha( const double& step, const double& max );

      /// Solve the eigenvalue problem for whilst stepping in ALPHA while c_i > tol
      void step_alpha_tol( const double& step, const double& tol );

      /// Solve the eigenvalue problem for whilst stepping (backwards) in ALPHA while c_i > tol
      void step_back_alpha_tol( const double& step, const double& tol );

      /// Iterate on the wavenumber to drive first eigenvalue to neutral
      void iterate_to_neutral( const double& epsilon = 1e-3 );

      /// Output the eigenvalues (and eigenvectors) to a directory
      void output()
      {
        make_output_directory();
        output_eigenvalues();
        if ( CALC_EIGENVECTORS ){
          output_eigenvectors();
        }
      }

      /// Output the eigenvalues to a directory
      void output_eigenvalues()
      {
        std::string eval_output_path( OUTPUT_PATH + "eigenvalues/" );
        EIGENVALUES.output( eval_output_path + "alpha_"
                      + Utility::stringify( ALPHA, 3 ) + "_evals.dat", 15 );
      }

      /// Output the eigenvectors to a directory
      void output_eigenvectors()
      {
        std::string evec_output_path( OUTPUT_PATH + "eigenvectors/" );

        EIGENVECTORS.dump_gnu( evec_output_path + "alpha_"
                      + Utility::stringify( ALPHA, 3 ) + "_evecs.dat" );
      }

      /// Make output directory
      void make_output_directory()
      {
        // Make eigenvalue directory
        std::string eval_output_path( OUTPUT_PATH + "eigenvalues/" );
        int status = mkdir( eval_output_path.c_str(), S_IRWXU );
        if ( status == 0 ) {
        std::cout << "  * Eigenvalue directory created successfully" << std::endl;
        }
        // Make eigenvector directory
        if ( CALC_EIGENVECTORS ){
          std::string evec_output_path( OUTPUT_PATH + "eigenvectors/" );
          int status = mkdir( evec_output_path.c_str(), S_IRWXU );
          if ( status == 0 ) {
          std::cout << "  * Eigenvector directory created successfully" << std::endl;
          }
        }
      }


  }; // End of class Rayleigh_2D

  void Rayleigh_2D::solve_evp()
  {
    //std::cout << "*** Setting up the generalised eigenvalue problem ***" << std::endl;
    //std::cout << "--- K = " << SSI.injection() << ", alpha = " << ALPHA << std::endl;
    SlepcInitialize(NULL,NULL,(char*)0,(char*)0);

    // Create the sparse matrices
    std::size_t N( SSI.hzeta_intervals() ); //TODO make N and M member variables?
    std::size_t M( SSI.eta_intervals() );
    std::size_t N_hzeta( N + 1 );
    std::size_t N_eta( M + 1 );
    Vector<double> ETA_NODES      = SSI.eta_nodes();
    Vector<double> HZETA_NODES    = SSI.hzeta_nodes();
    Vector<double> X_NODES        = SSI.x_nodes();
    Vector<double> Y_NODES        = SSI.y_nodes();
    Vector<double> BASE_ETA_NODES = SSI.base_eta_nodes();
    TwoD_node_mesh<double> sol = SSI.solution();
    OneD_node_mesh<double> base = SSI.base_flow_solution();

    SparseMatrix< std::complex<double> > A( N_eta * N_hzeta, N_eta * N_hzeta );
    SparseMatrix< std::complex<double> > B( N_eta * N_hzeta, N_eta * N_hzeta );

    // Create the SparseEigenSystem
    SparseEigenSystem< std::complex<double> > system( &A, &B );
    system.set_nev( EIGS_REQUESTED );
    system.set_region( LEFT, RIGHT, BOTTOM, TOP );
    system.set_target( TARGET );
    system.set_order( ORDER );
    system.calc_eigenvectors() = CALC_EIGENVECTORS;

    // Fill the sparse matrices
    std::size_t row( 0 );                               // Initialise row counter
    const double dY( SSI.y_nodes()[ 1 ] - SSI.y_nodes()[ 0 ] );
    const double dX( SSI.x_nodes()[ 1 ] - SSI.x_nodes()[ 0 ] );
    const double zeta0( SSI.injection_width() );

    // hzeta = 0 boundary ( left boundary )
    std::size_t i( 0 );

    for ( std::size_t j = 0; j < M + 1 ; ++j )
    {
      double hzeta( HZETA_NODES[ 0 ] );
      double Xd( SSI.mesh_Xd( hzeta ) );
      //double eta( ETA_NODES[ j ] );

      // P_hzeta = 0
      A( row, i * N_eta + j )           = -3 * Xd / ( 2 * dX );
      A( row, ( i + 1 ) * N_eta + j )   =  4 * Xd / ( 2 * dX );
      A( row, ( i + 2 ) * N_eta + j )   = -1 * Xd / ( 2 * dX );

      ++row;

    } // End of for loop over LHS eta nodes

    // Interior points between the hzeta boundaries
    for ( std::size_t i = 1; i < N; ++i )
    {
      // hzeta location
      double hzeta( HZETA_NODES[ i ] );
      double Xd( SSI.mesh_Xd( hzeta ) );
      double Xdd( SSI.mesh_Xdd( hzeta ) );

      // eta = 0 boundary ( bottom boundary )
      std::size_t j( 0 );
      double eta( ETA_NODES[ j ] );
      double Yd( SSI.mesh_Yd( eta ) );

      // P_eta = 0
      A( row, i * N_eta + j )        = -3 * Yd / ( 2 * dY );
      A( row, i * N_eta + j + 1 )    =  4 * Yd / ( 2 * dY );
      A( row, i * N_eta + j + 2 )    = -1 * Yd / ( 2 * dY );

      ++row;

      // Main interior grid points
      for ( std::size_t j = 1; j < M; ++j )
      {
        // eta location
        double eta( ETA_NODES[ j ] );
        double Yd( SSI.mesh_Yd( eta ) );
        double Ydd( SSI.mesh_Ydd( eta ) );

        // Self similar solution
        double U = sol( i, j, 6 ); // UB + U_pert
        double U_hzeta = Xd * ( sol( i + 1, j, 2 ) - sol( i - 1, j, 2 ) ) / ( 2 * dX ); // U_pert_hzeta
        double U_eta = Yd * ( sol( i, j + 1, 2 ) - sol( i, j - 1, 2 ) ) / ( 2 * dY ); // U_pert_eta
        U_eta += base.get_interpolated_vars( eta )[ 1 ]; // + UBd

        // Laplacian coefficients for finite-differencing

        // X(i,j-1)
        double laplace_1 =  ( Yd*Yd / (dY*dY) - Ydd / (2.*dY) ) ;
        // X(i-1,j)
        double laplace_3 = ( Xd*Xd / (dX*dX) - Xdd / (2.*dX) ) / ( zeta0 * zeta0 );
        // X(i,j)
        double laplace_4 = -2.*( Yd*Yd / (dY*dY)
                           + Xd*Xd / ( zeta0 * zeta0 * dX * dX ) );
        // X(i+1,j)
        double laplace_5 = ( Xdd / (2.*dX) + Xd*Xd / (dX*dX) ) / ( zeta0 * zeta0 );
        // X(i,j+1)
        double laplace_7 = ( Yd*Yd / (dY*dY) + Ydd / (2.*dY) );

        //////////////////////////
        // 2D Rayleigh equation //
        //////////////////////////

        // A matrix

        // Laplacian of P * U
        A( row, i * N_eta + j - 1 )    = laplace_1 * U;
        A( row, (i - 1) * N_eta + j )  = laplace_3 * U;
        A( row, i * N_eta + j )        = laplace_4 * U;
        A( row, (i + 1) * N_eta + j )  = laplace_5 * U;
        A( row, i * N_eta + j + 1 )    = laplace_7 * U;

        // - alpha^2 * U * P
        A( row, i * N_eta + j )       += - ALPHA * ALPHA * U;

        // - 2 * U_eta * P_eta
        A( row, i * N_eta + j + 1 )   += -2 * U_eta * ( Yd / ( 2 * dY ) );
        A( row, i * N_eta + j - 1 )   +=  2 * U_eta * ( Yd / ( 2 * dY ) );

        // - 2 * U_hzeta * P_hzeta / zeta0^2
        A( row, (i + 1) * N_eta + j ) += -2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );
        A( row, (i - 1) * N_eta + j ) +=  2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );

        // B matrix

        // Laplacian of P
        B( row, i * N_eta + j - 1 )    = laplace_1;
        B( row, (i - 1) * N_eta + j )  = laplace_3;
        B( row, i * N_eta + j )        = laplace_4;
        B( row, (i + 1) * N_eta + j )  = laplace_5;
        B( row, i * N_eta + j + 1 )    = laplace_7;

        // - alpha^2 * P
        B( row, i * N_eta + j )       += - ALPHA * ALPHA;

        ++row;

      }

      // eta = eta_inf boundary ( top boundary )
      j = M ;
      eta = ETA_NODES[ j ];
      Yd = SSI.mesh_Yd( eta );

      // P = 0 (Dirichlet)
      //A( row, i * N_eta + j ) = 1;

      // Robin condition (P_eta + alpha*P = 0)
      // P_eta
      A( row, i * N_eta + j )        =  3 * Yd / ( 2 * dY );
      A( row, i * N_eta + j - 1 )    = -4 * Yd / ( 2 * dY );
      A( row, i * N_eta + j - 2 )    =  1 * Yd / ( 2 * dY );

      // + alpha * P
      A( row, i * N_eta + j )       += ALPHA;

      ++row;

    } // End of for loop over interior nodes

    // hzeta = hzeta_inf boundary ( right boundary )
    for ( std::size_t j = 0; j < M + 1; ++j )
    {
      std::size_t i( N );
      // P = 0
      A( row, i * N_eta + j ) = 1;
      ++row;

    } // End of loop over RHS eta nodes

    // Solve the sparse eigenvalue problem
    try
    {
      system.eigensolve();
    }
    catch ( std::runtime_error )
    {
      throw Error( "Rayleigh_2D solve_evp() failed through exception being raised." );
    }

    EIGENVALUES = system.eigenvalues();
    // Convert the eigenvector matrix to a 2D mesh
    if ( CALC_EIGENVECTORS )
    {
      Matrix< std::complex<double> > evecs;
      evecs = system.eigenvectors();
      std::size_t nev( system.get_nconv() ); // Number of converged eigenvalues
      TwoD_node_mesh< std::complex<double> > output( HZETA_NODES, ETA_NODES, nev );
      EIGENVECTORS = output;
      // Push the data back into the unmapped domain
      for ( std::size_t n = 0; n < nev; ++n )
      {
        for ( std::size_t i = 0; i < N_hzeta; ++i )
        {
          for ( std::size_t j = 0; j < N_eta; ++j )
          {
            EIGENVECTORS( i, j, n ) = evecs( n, i * N_eta + j );
          }
        }
      }
    }

    //SlepcFinalize();
  }

  void Rayleigh_2D::step_in_alpha( const double& step, const double& max )
  {
    std::cout << "*** Iterating on the wavenumber ***" << std::endl;
    make_output_directory();

    TrackerFile metric( OUTPUT_PATH + "First_eval.dat" );
    std::complex<double> first_eval;
    metric.push_ptr( &ALPHA, "Wavenumber (alpha)");
    metric.push_ptr( &first_eval, "1st eigenvalue" );
    metric.header();

    Timer timer;
    timer.start();

    do {
      std::cout << "*** alpha = " << ALPHA << std::endl;
      solve_evp();
      first_eval = EIGENVALUES[ 0 ];
      metric.update();
      output_eigenvalues();
      if ( CALC_EIGENVECTORS ){ output_eigenvectors(); }
      TARGET = first_eval;
      ALPHA += step;

      timer.print();
      timer.stop();
      timer.reset();
      timer.start();

    }while ( ALPHA <= max );
  }

  void Rayleigh_2D::step_alpha_tol( const double& step, const double& tol )
  {
    std::cout << "*** Iterating on the wavenumber ***" << std::endl;
    make_output_directory();

    TrackerFile metric( OUTPUT_PATH + "First_eval.dat" );
    std::complex<double> first_eval;
    metric.push_ptr( &ALPHA, "Wavenumber (alpha)");
    metric.push_ptr( &first_eval, "1st eigenvalue" );
    metric.header();

    Timer timer;
    timer.start();

    do {
      std::cout << "*** alpha = " << ALPHA << std::endl;
      solve_evp();
      first_eval = EIGENVALUES[ 0 ];
      metric.update();
      output_eigenvalues();
      if ( CALC_EIGENVECTORS ){ output_eigenvectors(); }
      TARGET = first_eval;
      ALPHA += step;

      timer.print();
      timer.stop();
      timer.reset();
      timer.start();

    }while ( first_eval.imag() > tol );
  }

  void Rayleigh_2D::step_back_alpha_tol( const double& step, const double& tol )
  {
    std::cout << "*** Iterating on the wavenumber ***" << std::endl;
    make_output_directory();

    TrackerFile metric( OUTPUT_PATH + "First_eval_back.dat" );
    std::complex<double> first_eval;
    metric.push_ptr( &ALPHA, "Wavenumber (alpha)");
    metric.push_ptr( &first_eval, "1st eigenvalue" );
    metric.header();

    Timer timer;
    timer.start();

    do {
      std::cout << "*** alpha = " << ALPHA << std::endl;
      solve_evp();
      first_eval = EIGENVALUES[ 0 ];
      metric.update();
      output_eigenvalues();
      if ( CALC_EIGENVECTORS ){ output_eigenvectors(); }
      TARGET = first_eval;
      ALPHA -= step;

      timer.print();
      timer.stop();
      timer.reset();
      timer.start();

    }while ( first_eval.imag() > tol && ALPHA >= 0.0 );
  }

  void Rayleigh_2D::iterate_to_neutral( const double& epsilon )
  {
    std::cout << "*** Finding the critical wavenumber for K = " << SSI.injection() << std::endl;
    double delta( 1.e-8 );
    do
    {
      std::cout << "*** alpha = " << ALPHA << std::endl;
      solve_evp();
      std::complex<double> copy_of_ev( EIGENVALUES[ 0 ] );
      ALPHA += delta;
      solve_evp();
      ALPHA -= delta;
      double d_ev = ( std::imag( EIGENVALUES[ 0 ] ) - std::imag( copy_of_ev ) ) / delta;
      ALPHA -= std::imag( copy_of_ev ) / d_ev;
      std::cout << "ITERATING: " << ALPHA << " " << EIGENVALUES[ 0 ] << " " << d_ev << "\n";
      TARGET = EIGENVALUES[ 0 ];
    }
    while ( std::abs( std::imag( EIGENVALUES[ 0 ] ) ) > epsilon );
  }

} // End of namespace TSL

#endif
