/* OrrSommerfeld_2D - Here we define the OrrSommerfeld_2D class which is useful for solving
              the OrrSommerfeld equation.
*/

#ifndef ORRSOMMERFELD_2D_H
#define ORRSOMMERFELD_2D_H

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

#include "Eigenvalue"

enum{ v, w, q, s };

namespace TSL
{
  class OrrSommerfeld_2D
  {
    protected:
      SelfSimInjection SSI;                                   // Self similar injection object
      Vector< std::complex<double> > EIGENVALUES;             // Eigenvalues
      TwoD_node_mesh< std::complex<double> > EIGENVECTORS;    // Eigenvectors
      double ALPHA;                                           // Wavenumber
      double RX;                                              // Local Reynolds number
      std::size_t EIGS_REQUESTED;                             // Number of eigenvalues requested
      double LEFT;                                            // Region in which to look for
      double RIGHT;                                           // eigenvalues
      double BOTTOM;
      double TOP;
      std::complex<double> TARGET;                            // Target eigenvalue
      std::string ORDER;                                      // Ordering of returned eigenvalues
      bool CALC_EIGENVECTORS;                                 // Calculate the eigenvectors?
      std::string OUTPUT_PATH;

      std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
      {
        // Return the column number for the kth variable at node (i,j)
        return 4 * ( i * ( SSI.eta_intervals() + 1 ) + j ) + k;
        //TODO might not be the most efficient way of getting M / eta_intervals - maybe have a member M
      }

    public:

      /// Constructor
      OrrSommerfeld_2D( SelfSimInjection& ssi, double& alpha, double& rx, std::size_t& nev )
      {
        SSI = ssi;
        ALPHA = alpha;
        RX = rx;
        EIGS_REQUESTED = nev;
        //TODO give default values to other members
        CALC_EIGENVECTORS = true;
        OUTPUT_PATH = SSI.output_path();
      }

      /// Destructor
	   	~OrrSommerfeld_2D()
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
      double& alpha()
      {
        return ALPHA;
      }

      /// Return a pointer to the local Reynolds number
      double& Reynolds()
      {
        return RX;
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

  }; // End of class OrrSommerfeld_2D

  void OrrSommerfeld_2D::solve_evp()
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
    double beta = SSI.hartree();

    SparseMatrix< std::complex<double> > A( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta );
    SparseMatrix< std::complex<double> > B( 4 * N_eta * N_hzeta, 4 * N_eta * N_hzeta );

    // Create the SparseEigenSystem
    SparseEigenSystem< std::complex<double> > system( &A, &B );
    system.set_nev( EIGS_REQUESTED );
    system.set_region( LEFT, RIGHT, BOTTOM, TOP );
    system.set_target( TARGET );
    system.set_order( ORDER );
    system.calc_eigenvectors() = CALC_EIGENVECTORS;

    //Eigensystem< std::complex<double> > system;

    // Fill the sparse matrices
    std::size_t row( 0 );                               // Initialise row counter
    const double dY( SSI.y_nodes()[ 1 ] - SSI.y_nodes()[ 0 ] );
    const double dX( SSI.x_nodes()[ 1 ] - SSI.x_nodes()[ 0 ] );
    const double zeta0( SSI.injection_width() );

    std::complex<double> iaR ( 0.0, 1.0 / ( ALPHA * sqrt( RX ) ) );
    std::complex<double> Ralpha ( 1.0 / sqrt( RX ), ALPHA );

    // hzeta = 0 boundary ( left boundary )
    std::size_t i( 0 );

    for ( std::size_t j = 0; j < M + 1 ; ++j )
    {
      double hzeta( HZETA_NODES[ 0 ] );
      double Xd( SSI.mesh_Xd( hzeta ) );
      double Xdd( SSI.mesh_Xdd( hzeta ) );
      // eta location
      double eta( ETA_NODES[ j ] );
      double Yd( SSI.mesh_Yd( eta ) );
      double Ydd( SSI.mesh_Ydd( eta ) );

      // Symmetry boundary conditions

      // w = 0
      A( row, col( i, j, w ) ) = 1;

      ++row;

      // v_hzeta = 0
      A( row, col( i, j, v ) )       = - 3 * Xd / ( 2 * dX );
      A( row, col( i + 1, j, v ) )   =   4 * Xd / ( 2 * dX );
      A( row, col( i + 2, j, v ) )   = - 1 * Xd / ( 2 * dX );

      ++row;

      // s = 0
      A( row, col( 0, j, s ) ) = 1;

      ++row;

      // q_{hzeta} = 0
      A( row, col( i, j, q ) )      = - 3 * Xd / ( 2 * dX );
      A( row, col( i + 1, j, q ) )  =   4 * Xd / ( 2 * dX );
      A( row, col( i + 2, j, q ) )  = - 1 * Xd / ( 2 * dX );

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
      double Ydd( SSI.mesh_Ydd( eta ) );

      // v = 0
      A( row, col( i, 0, v ) ) = 1;

      ++row;

      // w = 0
      A( row, col( i, 0, w ) ) = 1;

      ++row;

      //s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
      A( row, col( i, 0, s ) )    =   1;
      A( row, col( i, 1, w ) )    = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) );
      A( row, col( i, 2, w ) )    =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) );
      A( row, col( i, 3, w ) )    =   iaR * Yd * Yd / ( dY * dY );

      ++row;

      //q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
      //A( row, col( i, 0, q ) )    =   1;
      //A( row, col( i, 1, v ) )    = - 2. * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, 0, q ) )    =   1;
      A( row, col( i, 1, v ) )    = - 6. * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, 2, v ) )    =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, 3, v ) )    = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

      ++row;

      // Main interior grid points
      for ( std::size_t j = 1; j < M; ++j )
      {
        // eta location
        double eta( ETA_NODES[ j ] );
        double Yd( SSI.mesh_Yd( eta ) );
        double Ydd( SSI.mesh_Ydd( eta ) );

        // Laplacian coefficients for finite-differencing

        // X(i,j-1)
        double laplace_1 =  ( Yd*Yd / (dY*dY) - Ydd / (2.*dY) ) ;
        // X(i-1,j)
        double laplace_3 = ( Xd*Xd / (dX*dX) - Xdd / (2.*dX) ) / ( zeta0 * zeta0 );
        // X(i,j)
        double laplace_4 = - 2 *( Yd*Yd / (dY*dY) + Xd*Xd / ( zeta0 * zeta0 * dX * dX ) )
                           - ALPHA * ALPHA;
        // X(i+1,j)
        double laplace_5 = ( Xdd / (2.*dX) + Xd*Xd / (dX*dX) ) / ( zeta0 * zeta0 );
        // X(i,j+1)
        double laplace_7 = ( Yd*Yd / (dY*dY) + Ydd / (2.*dY) );

        // Self similar solution
        double U = sol( i, j, 6 ); // UB + U_pert
        double U_hzeta = Xd * ( sol( i + 1, j, 2 ) - sol( i - 1, j, 2 ) ) / ( 2 * dX ); // U_pert_hzeta
        double U_eta = Yd * ( sol( i, j + 1, 2 ) - sol( i, j - 1, 2 ) ) / ( 2 * dY ); // U_pert_eta
        double UBd = base.get_interpolated_vars( eta )[ 1 ];
        U_eta += UBd;

        double UB = base.get_interpolated_vars( eta )[ 0 ];
        double PhiB = base.get_interpolated_vars( eta )[ 2 ];
        double U_eta_eta = laplace_1 * sol( i, j - 1, 2 ) + laplace_7 * sol( i, j + 1, 2 );
        U_eta_eta += - 2 * Yd * Yd * sol( i, j, 2 ) / ( dY * dY );
        U_eta_eta += beta * ( UB * UB - 1 ) - PhiB * UBd; // + UBdd = beta * [UB^2 - 1 ]- PhiB * UBd

        double U_hzeta_hzeta = laplace_3 * sol( i - 1, j, 2 ) + laplace_5 * sol( i + 1, j, 2 );
        U_hzeta_hzeta += - 2 * Xd * Xd * sol( i, j, 2 ) / ( dX * dX );

        double U_eta_hzeta = Yd * Xd * ( sol( i + 1, j + 1, 2 ) + sol( i - 1, j - 1, 2 )
                                       - sol( i + 1, j - 1, 2 ) - sol( i - 1, j + 1, 2 ) )
                                       / ( 4 * dY * dX );

        ///////////////////////////////
        // 2D OrrSommerfeld equation //
        ///////////////////////////////

        ////////////////
        // v equation //
        ////////////////

        // A matrix

        // ( i / alpha * Rx^(1/2) ) * Laplacian of v
        A( row, col( i, j - 1, v ) )     = iaR * laplace_1;
        A( row, col( i - 1, j, v ) )     = iaR * laplace_3;
        A( row, col( i, j, v ) )         = iaR * laplace_4;
        A( row, col( i + 1, j, v ) )     = iaR * laplace_5;
        A( row, col( i, j + 1, v ) )     = iaR * laplace_7;

        // + (1 - i / (alpha*Rx^(1/2)) ) * U * v
        //A( row, col( i, j, v ) )        +=  U;
        A( row, col( i, j, v ) )        += ( 1.0 - iaR ) * U;

        // - ((1 + beta) / Rx) * v
        A( row, col( i, j, v ) )        += - ( 1.0 + beta ) / RX;

        // + (beta - 1) * ( i / alpha * Rx^(3/2) ) * v
        A( row, col( i, j, v ) )        += ( beta - 1.0 ) * iaR / RX;

        // - q
        A( row, col( i, j, q ) )         = - 1.;

        // B matrix

        B( row, col( i, j, v ) ) =  1;

        ++row;

        ////////////////
        // w equation //
        ////////////////

        // A matrix

        // ( i / alpha * Rx^(1/2) ) * Laplacian of w
        A( row, col( i, j - 1, w ) )     = iaR * laplace_1;
        A( row, col( i - 1, j, w ) )     = iaR * laplace_3;
        A( row, col( i, j, w ) )         = iaR * laplace_4;
        A( row, col( i + 1, j, w ) )     = iaR * laplace_5;
        A( row, col( i, j + 1, w ) )     = iaR * laplace_7;

        // + (1 - i / (alpha*Rx^(1/2)) ) * U * w
        //A( row, col( i, j, w ) )        +=  U;
        A( row, col( i, j, w ) )        += ( 1.0 - iaR ) * U;

        // - ((1 + beta) / Rx) * w
        A( row, col( i, j, w ) )        += - ( 1.0 + beta ) / RX;

        // + (beta - 1) * ( i / alpha * Rx^(3/2) ) * w
        A( row, col( i, j, w ) )        += ( beta - 1.0 ) * iaR / RX;

        // - s
        A( row, col( i, j, s ) )         = - 1.;

        // B matrix

        B( row, col( i, j, w ) ) =  1;

        ++row;

        ////////////////
        // q equation //
        ////////////////

        // A matrix

        // q_{eta eta}
        A( row, col( i, j - 1, q ) )      =   Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ;
        A( row, col( i, j, q ) )          = - 2 * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 1, q ) )      =   Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ;

        // + ( i * alpha + Rx^(-1/2) )^2 * q
        //A( row, col( i, j, q ) )         += - ALPHA * ALPHA;
        A( row, col( i, j, q ) )         += Ralpha * Ralpha;

        // + ( 1 / zeta0 ) * s_{hzeta eta}
        A( row, col( i + 1, j + 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j - 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i + 1, j - 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j + 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );

        // + (1 - i / (alpha*Rx^(1/2)) ) * ( 1 / zeta0 ) * U_{eta} * w_{hzeta}
        A( row, col( i + 1, j, w ) )     =   ( 1.0 - iaR ) * U_eta * Xd / ( 2 * dX * zeta0 );
        A( row, col( i - 1, j, w ) )     = - ( 1.0 - iaR ) * U_eta * Xd / ( 2 * dX * zeta0 );
        //A( row, col( i + 1, j, w ) )     =   U_eta * Xd / ( 2 * dX * zeta0 );
        //A( row, col( i - 1, j, w ) )     = - U_eta * Xd / ( 2 * dX * zeta0 );

        // - (1 - i / (alpha*Rx^(1/2)) ) * ( 1 / zeta0 ) * U_{hzeta} * w_{eta}
        A( row, col( i, j + 1, w ) )     = - ( 1.0 - iaR ) * U_hzeta * Yd / ( 2 * dY * zeta0 );
        A( row, col( i, j - 1, w ) )     =   ( 1.0 - iaR ) * U_hzeta * Yd / ( 2 * dY * zeta0 );
        //A( row, col( i, j + 1, w ) )     = - U_hzeta * Yd / ( 2 * dY * zeta0 );
        //A( row, col( i, j - 1, w ) )     =   U_hzeta * Yd / ( 2 * dY * zeta0 );

        // - (1 - i / (alpha*Rx^(1/2)) ) * ( 1 / zeta0 ) * U_{hzeta eta} * w
        A( row, col( i, j, w ) )         = - ( 1.0 - iaR ) * U_eta_hzeta / zeta0;
        //A( row, col( i, j, w ) )         = - U_eta_hzeta / zeta0;

        // - (1 - i / (alpha*Rx^(1/2)) ) * U_{eta eta} * v
        A( row, col( i, j, v ) )          = - ( 1.0 - iaR ) * U_eta_eta;
        //A( row, col( i, j, v ) )          = - U_eta_eta;

/*
        // Rich's version
        // + ( 1 / zeta0 ) * s_{hzeta eta}
        A( row, col( i + 1, j + 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j - 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i + 1, j - 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j + 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );

        // - 2 * U_{hzeta} * w_{eta} / zeta0
        A( row, col( i, j + 1, w ) )     = - 2 * U_hzeta * Yd / ( 2 * dY * zeta0 );
        A( row, col( i, j - 1, w ) )     =   2 * U_hzeta * Yd / ( 2 * dY * zeta0 );

        // - 2 * U_{eta} * v_{eta}
        A( row, col( i, j + 1, w ) )     = - 2 * U_eta * Yd / ( 2 * dY );
        A( row, col( i, j - 1, w ) )     =   2 * U_eta * Yd / ( 2 * dY );

        // - 2 * U_{eta hzeta} * w / zeta0
        A( row, col( i, j, w ) )         = - 2 * U_eta_hzeta / zeta0;

        // - 2 * U_{eta eta} * v
        A( row, col( i, j, v ) )         = - 2 * U_eta_eta;
*/
        ++row;

        ////////////////
        // s equation //
        ////////////////

        // A matrix

        // ( 1 / zeta0^2 ) * s_{hzeta hzeta}
        A( row, col( i - 1, j, s ) )      =   ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) )
                                            / ( zeta0 * zeta0 );
        A( row, col( i, j, s ) )          = - 2 * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
        A( row, col( i + 1, j, s ) )      =   ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) )
                                            / ( zeta0 * zeta0 );

        // + ( i * alpha + Rx^(-1/2) )^2 * s
        //A( row, col( i, j, s ) )         += - ALPHA * ALPHA;
        A( row, col( i, j, s ) )         += Ralpha * Ralpha;

        // + ( 1 / zeta0 ) * q_{hzeta eta}
        A( row, col( i + 1, j + 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j - 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i + 1, j - 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j + 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );

        // + (1 - i / (alpha*Rx^(1/2)) ) * ( 1 / zeta0 ) * U_{hzeta} * v_{eta}
        A( row, col( i, j + 1, v ) )     =   ( 1.0 - iaR ) * U_hzeta * Yd / ( 2 * dY * zeta0 );
        A( row, col( i, j - 1, v ) )     = - ( 1.0 - iaR ) * U_hzeta * Yd / ( 2 * dY * zeta0 );
        //A( row, col( i, j + 1, v ) )     =   U_hzeta * Yd / ( 2 * dY * zeta0 );
        //A( row, col( i, j - 1, v ) )     = - U_hzeta * Yd / ( 2 * dY * zeta0 );

        // - (1 - i / (alpha*Rx^(1/2)) ) * ( 1 / zeta0 ) * U_{eta} * v_{hzeta}
        A( row, col( i + 1, j, v ) )     = - ( 1.0 - iaR ) * U_eta * Xd / ( 2 * dX * zeta0 );
        A( row, col( i - 1, j, v ) )     =   ( 1.0 - iaR ) * U_eta * Xd / ( 2 * dX * zeta0 );
        //A( row, col( i + 1, j, v ) )     = - U_eta * Xd / ( 2 * dX * zeta0 );
        //A( row, col( i - 1, j, v ) )     =   U_eta * Xd / ( 2 * dX * zeta0 );

        // - (1 - i / (alpha*Rx^(1/2)) ) * ( 1 / zeta0 ) * U_{eta hzeta} * v
        A( row, col( i, j, v ) )         = - ( 1.0 - iaR ) * U_eta_hzeta / zeta0;
        //A( row, col( i, j, v ) )         = - U_eta_hzeta / zeta0;

        // - (1 - i / (alpha*Rx^(1/2)) ) * ( 1 / zeta0^2 ) * U_{hzeta hzeta} * w
        A( row, col( i, j, w ) )         = - ( 1.0 - iaR ) * U_hzeta_hzeta / ( zeta0 * zeta0 );
        //A( row, col( i, j, w ) )         = - U_hzeta_hzeta / ( zeta0 * zeta0 );

/*
        // Rich's version
        // + ( 1 / zeta0 ) * q_{hzeta eta}
        A( row, col( i + 1, j + 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j - 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i + 1, j - 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );
        A( row, col( i - 1, j + 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * zeta0 );

        // - 2 * U_{hzeta} * w_{hzeta} / zeta0^2
        A( row, col( i + 1, j, w ) )     = - 2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );
        A( row, col( i - 1, j, w ) )     =   2 * U_hzeta * Xd / ( 2 * dX * zeta0 * zeta0 );

        // - 2 * U_{eta} * v_{hzeta} / zeta0
        A( row, col( i + 1, j, w ) )     = - 2 * U_eta * Xd / ( 2 * dX * zeta0 );
        A( row, col( i - 1, j, w ) )     =   2 * U_eta * Xd / ( 2 * dX * zeta0 );

        // - 2 * U_{eta hzeta} * v / zeta0
        A( row, col( i, j, v ) )         = - 2 * U_eta_hzeta / zeta0;

        // - 2 * U_{hzeta hzeta} * w / zeta0^2
        A( row, col( i, j, w ) )         = - 2 * U_hzeta_hzeta / ( zeta0 * zeta0 );
*/
        ++row;

      }

      // eta = eta_inf boundary ( top boundary )
      j = M ;
      eta = ETA_NODES[ j ];
      Yd = SSI.mesh_Yd( eta );

      // v = 0
      A( row, col( i, j, v ) ) = 1;

      ++row;

      // w = 0
      A( row, col( i, j, w ) ) = 1;

      ++row;

      // s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
      A( row, col( i, j, s ) )        =   1;
      A( row, col( i, j - 1, w ) )    = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) );
      A( row, col( i, j - 2, w ) )    =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) );
      A( row, col( i, j - 3, w ) )    =   iaR * Yd * Yd / ( dY * dY );

      ++row;

      // q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
      //A( row, col( i, j, q ) )        =   1;
      //A( row, col( i, j - 1, v ) )    = - 2. * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, j, q ) )        =   1;
      A( row, col( i, j - 1, v ) )    = - 6. * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, j - 2, v ) )    =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
      A( row, col( i, j - 3, v ) )    = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

      ++row;

    } // End of for loop over interior nodes

    // hzeta = hzeta_inf boundary ( right boundary )
    for ( std::size_t j = 0; j < M + 1; ++j )
    {
      std::size_t i( N );
      double hzeta( HZETA_NODES[ i ] );
      double Xd( SSI.mesh_Xd( hzeta ) );
      double Xdd( SSI.mesh_Xdd( hzeta ) );

      // eta location
      double eta( ETA_NODES[ j ] );
      double Yd( SSI.mesh_Yd( eta ) );
      double Ydd( SSI.mesh_Ydd( eta ) );

      // w = 0
      A( row, col( i, j, w ) ) = 1;

      ++row;

      // v = 0
      A( row, col( i, j, v ) ) = 1;

      ++row;

      // q - (1/zeta0^2) * (i / (alpha*Rx^(1/2))) * v_{hzeta hzeta} = 0
      A( row, col( i, j, q ) )        =   1;
      A( row, col( i - 1, j, v ) )    = - iaR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) )
                                    / ( zeta0 * zeta0 );
      A( row, col( i - 2, j, v ) )    =   iaR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) )
                                    / ( zeta0 * zeta0 );
      A( row, col( i - 3, j, v ) )    =   iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );

      ++row;

      // s - (1/zeta0^2) * (i / (alpha*Rx^(1/2))) * w_{hzeta hzeta} = 0
      //A( row, col( i, j, s ) )        =   1;
      //A( row, col( i - 1, j, w ) )    = - 2. * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
      A( row, col( i, j, s ) )        =   1;
      A( row, col( i - 1, j, w ) )    = - 6. * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
      A( row, col( i - 2, j, w ) )    =   ( 3. / 2. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );
      A( row, col( i - 3, j, w ) )    = - ( 2. / 9. ) * iaR * Xd * Xd / ( dX * dX * zeta0 * zeta0 );

      ++row;

    } // End of loop over RHS eta nodes


    // Solve the sparse eigenvalue problem
    try
    {
      system.eigensolve();
    }
    catch ( std::runtime_error )
    {
      throw Error( "OrrSommerfeld_2D solve_evp() failed through exception being raised." );
    }

    EIGENVALUES = system.eigenvalues();
    // Convert the eigenvector matrix to a 2D mesh


    //TODO need nev * 4 for v,w,q,s eigenvectors
    if ( CALC_EIGENVECTORS )
    {
      Matrix< std::complex<double> > evecs;
      evecs = system.eigenvectors();
      std::size_t nev( system.get_nconv() ); // Number of converged eigenvalues
      TwoD_node_mesh< std::complex<double> > output( HZETA_NODES, ETA_NODES, nev * 4 );
      EIGENVECTORS = output;
      // Push the data back into the unmapped domain
      for ( std::size_t n = 0; n < nev; ++n )
      {
        for ( std::size_t i = 0; i < N_hzeta; ++i )
        {
          for ( std::size_t j = 0; j < N_eta; ++j )
          {
            EIGENVECTORS( i, j, 4 * n + v ) = evecs( n, 4 * ( i * N_eta + j ) + v );
            EIGENVECTORS( i, j, 4 * n + w ) = evecs( n, 4 * ( i * N_eta + j ) + w );
            EIGENVECTORS( i, j, 4 * n + q ) = evecs( n, 4 * ( i * N_eta + j ) + q );
            EIGENVECTORS( i, j, 4 * n + s ) = evecs( n, 4 * ( i * N_eta + j ) + s );
          }
        }
      }
    }
/*
    // Check BCs
    // hzeta = 0
    std::cout << "w(hzeta=0) = " << std::endl;
    for ( std::size_t j = 0; j < N_eta; ++j )
    {
      std::cout << EIGENVECTORS( 0, j, w ) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "v_hzeta(hzeta=0) = " << std::endl;
    double hzeta( HZETA_NODES[ 0 ] );
    double Xd( SSI.mesh_Xd( hzeta ) );
    for ( std::size_t j = 0; j < N_eta; ++j )
    {
      std::cout << (- 3. * EIGENVECTORS( 0, j, v )
                    + 4. * EIGENVECTORS( 1, j, v ) - EIGENVECTORS( 2, j, v ) )
                    * Xd / ( 2. * dX ) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "w_{hzeta hzeta}(hzeta=0) = " << std::endl;
    double Xdd( SSI.mesh_Xdd( hzeta ) );
    for ( std::size_t j = 0; j < N_eta; ++j )
    {
      std::complex<double> ans;
      ans = ( 2 * Xd * Xd / ( dX * dX ) - 3 * Xdd / ( 2 * dX ) ) * EIGENVECTORS( 0, j, w );
      ans += ( - 5 * Xd * Xd / ( dX * dX ) + 4 * Xdd / ( 2 * dX ) ) * EIGENVECTORS( 1, j, w );
      ans += ( 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) ) * EIGENVECTORS( 2, j, w );
      ans += ( - 1 * Xd * Xd / ( dX * dX ) ) * EIGENVECTORS( 3, j, w );

      std::cout << ans << std::endl;
    }

    std::cout << std::endl;


    //eta = 0
    std::cout << "v(eta=0) = " << std::endl;
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      std::cout << EIGENVECTORS( i, 0, v ) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "w(eta=0) = " << std::endl;
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      std::cout << EIGENVECTORS( i, 0, w ) << std::endl;
    }
    std::cout << std::endl;

    double eta( ETA_NODES[ 0 ] );
    double Yd( SSI.mesh_Yd( eta ) );
    std::cout << "v_eta(eta=0) = " << std::endl;
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      std::complex<double> ans;
      ans = ( -3 * Yd / ( 2 * dY ) ) * EIGENVECTORS( i, 0, v );
      ans += ( 4 * Yd / ( 2 * dY ) ) * EIGENVECTORS( i, 1, v );
      ans += ( -1 * Yd / ( 2 * dY ) ) * EIGENVECTORS( i, 2, v );

      std::cout << ans << std::endl;
    }
    std::cout << std::endl;

    //eta = eta_inf
    std::cout << "v(eta=eta_inf) = " << std::endl;
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      std::cout << EIGENVECTORS( i, M, v ) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "w(eta=eta_inf) = " << std::endl;
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      std::cout << EIGENVECTORS( i, M, w ) << std::endl;
    }
    std::cout << std::endl;

    eta = ETA_NODES[ M ];
    Yd = SSI.mesh_Yd( eta );

    std::cout << "v_eta(eta=eta_inf) = " << std::endl;
    for ( std::size_t i = 0; i < N_hzeta; ++i )
    {
      std::complex<double> ans;
      ans = ( 3 * Yd / ( 2 * dY ) ) * EIGENVECTORS( i, M, v );
      ans += ( -4 * Yd / ( 2 * dY ) ) * EIGENVECTORS( i, M - 1, v );
      ans += ( 1 * Yd / ( 2 * dY ) ) * EIGENVECTORS( i, M - 2, v );

      std::cout << ans << std::endl;
    }
    std::cout << std::endl;


    //hzeta = hzeta_inf
    hzeta = HZETA_NODES[ N ];
    Xd = SSI.mesh_Xd( hzeta );

    std::cout << "w(hzeta=hzeta_inf) = " << std::endl;
    for ( std::size_t j = 0; j < N_eta; ++j )
    {
      std::cout << EIGENVECTORS( N, j, w ) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "v(hzeta=hzeta_inf) = " << std::endl;
    for ( std::size_t j = 0; j < N_eta; ++j )
    {
      std::cout << EIGENVECTORS( N, j, v ) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "w_hzeta(hzeta=hzeta_inf) = " << std::endl;
    for ( std::size_t j = 0; j < N_eta; ++j )
    {
      std::complex<double> ans;
      ans = ( 3 * Xd / ( 2 * dX ) ) * EIGENVECTORS( N, j, w );
      ans += ( -4 * Xd / ( 2 * dX ) ) * EIGENVECTORS( N - 1, j, w );
      ans += ( 1 * Xd / ( 2 * dX ) ) * EIGENVECTORS( N - 2, j, w );

      std::cout << ans << std::endl;
    }*/

    //SlepcFinalize();
  }

  void OrrSommerfeld_2D::step_in_alpha( const double& step, const double& max )
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

  void OrrSommerfeld_2D::iterate_to_neutral( const double& epsilon )
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
