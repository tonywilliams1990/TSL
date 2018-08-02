/* VWI - Here we define the VWI class which is useful for solving
         the vortex wave interaction equations.
*/

#ifndef VWI_H
#define VWI_H

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

enum class VWI_base{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };
enum class VWI_enum{ v, w, q, s, Phi, Psi, U, Theta };

namespace TSL
{
  class VWI
  {
    protected:
      SelfSimInjection SSI;          // Self similar injection object
      double ALPHA;                  // Wavenumber
      double RX;                     // Local Reynolds number
      double SIGMA;                  // Wave amplitude
      double BETA;                   // Hartree parameter
      double ZETA0;                  // Injection width
      double K;                      // Injection parameter ( +ve = blowing )
      std::complex<double> C_GUESS;  // Current eigenvalue guess
      std::string OUTPUT_PATH;

      // Mesh
      std::size_t N;                 // Number of intervals in the zeta_hat direction
      std::size_t M;                 // Number of intervals in the eta direction
      Vector<double> ETA_NODES;
      Vector<double> HZETA_NODES;
      Vector<double> X_NODES;
      Vector<double> Y_NODES;
      Vector<double> BASE_ETA_NODES;

      // Solution
      OneD_node_mesh<double> BASE_SOLUTION;            // Base flow ODE solution
      TwoD_node_mesh< std::complex<double> > Q;        // Current guess mesh
      TwoD_node_mesh< std::complex<double> > Q_OUTPUT; // Output mesh

      std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
      {
        // Return the column number for the kth variable at node (i,j)
        return 8 * ( i * ( SSI.eta_intervals() + 1 ) + j ) + k;
      }

    public:

      /// Constructor
      VWI( SelfSimInjection& ssi, double& alpha, double& rx, double& sigma )
      {
        SSI   = ssi;
        ALPHA = alpha;
        RX    = rx;
        BETA  = SSI.hartree();
        ZETA0 = SSI.injection_width();
        K     = SSI.injection();

        N              = SSI.hzeta_intervals();
        M              = SSI.eta_intervals();
        ETA_NODES      = SSI.eta_nodes();
        HZETA_NODES    = SSI.hzeta_nodes();
        X_NODES        = SSI.x_nodes();
        Y_NODES        = SSI.y_nodes();
        BASE_ETA_NODES = SSI.base_eta_nodes();


        BASE_SOLUTION = SSI.base_flow_solution();
        TwoD_node_mesh< std::complex<double> > q( X_NODES, Y_NODES, 8 );
        TwoD_node_mesh< std::complex<double> > q_output( HZETA_NODES, ETA_NODES, 12 );
        Q        = q;
        Q_OUTPUT = q_output;

        OUTPUT_PATH = SSI.output_path();

      }

      /// Destructor
	   	~VWI()
      {
        SlepcFinalize();
      }

      /* ----- Methods ----- */

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

      /// Return a pointer to current eigenvalue guess
      std::complex<double>& c_guess()
      {
        return C_GUESS;
      }

      /// Set an initial guess
      void set_guess( const TwoD_node_mesh< std::complex<double> >& guess )
      {
        //TODO error checking for size of mesh and number of variables
        //std::cout << "X_NODES.size() = " << X_NODES.size() << std::endl;
        for ( std::size_t i = 0; i < X_NODES.size(); ++i )
        {
          for ( std::size_t j = 0; j < Y_NODES.size(); ++j )
          {
            Q( i, j, 0 ) = guess( i, j, 0 );
            Q( i, j, 1 ) = guess( i, j, 1 );
            Q( i, j, 2 ) = guess( i, j, 2 );
            Q( i, j, 3 ) = guess( i, j, 3 );
            Q( i, j, 4 ) = guess( i, j, 4 );
            Q( i, j, 5 ) = guess( i, j, 5 );
            Q( i, j, 6 ) = guess( i, j, 6 );
            Q( i, j, 7 ) = guess( i, j, 7 );
          }
        }
      }

      /// Return a pointer to the current guess mesh
      TwoD_node_mesh< std::complex<double> >& current_guess()
      {
        return Q;
      }

      /// Return the solution mesh
      TwoD_node_mesh< std::complex<double> > solution() {
        return Q_OUTPUT;
      }

      /// Solve the sparse eigenvalue problem using Newton iteration
      void solve_local();

      // Virtual function for defining the transpiration function
      virtual double Phi_w_func( const double& hzeta )
      {
        throw Error( "--- Phi_w function is not defined ---" );
      }

  }; // End of class OrrSommerfeld_2D

  void VWI::solve_local()
  {
    SlepcInitialize( NULL, NULL, (char*)0, (char*)0 );
    int v     = static_cast<int>(VWI_enum::v);
    int w     = static_cast<int>(VWI_enum::w);
    int q     = static_cast<int>(VWI_enum::q);
    int s     = static_cast<int>(VWI_enum::s);
    int Phi   = static_cast<int>(VWI_enum::Phi);
    int Psi   = static_cast<int>(VWI_enum::Psi);
    int U     = static_cast<int>(VWI_enum::U);
    int Theta = static_cast<int>(VWI_enum::Theta);

    int UB      = static_cast<int>(VWI_base::UB);
    int UBd     = static_cast<int>(VWI_base::UBd);
    int PhiB    = static_cast<int>(VWI_base::PhiB);
    int ThetaB  = static_cast<int>(VWI_base::ThetaB);
    int ThetaBd = static_cast<int>(VWI_base::ThetaBd);
    int PsiB    = static_cast<int>(VWI_base::PsiB);


    // Create the sparse matrices
    std::size_t N_hzeta( N + 1 );
    std::size_t N_eta( M + 1 );
    std::size_t size( 8 * N_eta * N_hzeta + 1 );
    Vector< std::complex<double> > B( size, 0.0 );;

    // Step sizes
    const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
    const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );

    // Normalise based on the extra condition q(0,0) = 1 or q_eta(0,0) = 1
    //std::complex<double> lambda( Q( 0, 0, q ) );
    std::complex<double> lambda;
    double eta( ETA_NODES[ 0 ] );
    double Yd( SSI.mesh_Yd( eta ) );
    lambda = ( 3 * Yd / ( 2 * dY ) ) * Q( 0, 0, q )
           - ( 4 * Yd / ( 2 * dY ) ) * Q( 0, 1, q )
           + ( 1 * Yd / ( 2 * dY ) ) * Q( 0, 2, q );

    for ( std::size_t i = 0; i < N + 1; ++i )
    {
      for ( std::size_t j = 0; j < M + 1; ++j )
      {
        Q( i, j, v )  = Q( i, j, v ) / lambda;
        Q( i, j, w )  = Q( i, j, w ) / lambda;
        Q( i, j, q )  = Q( i, j, q ) / lambda;
        Q( i, j, s )  = Q( i, j, s ) / lambda;
      }
    }

    // Iterate to a solution
    double max_residual( 0.0 );             // Maximum residual
    std::size_t iteration( 0 );             // Initialise iteration counter
    std::size_t max_iterations( 20 );       // Maximum number of iterations

    std::complex<double> iaR ( 0.0, 1.0 / ( ALPHA * sqrt( RX ) ) );

    do{
      SparseMatrix< std::complex<double> > A( size, size );
      std::cout << "  * Assembling sparse matrix problem" << std::endl;

      Timer timer;                                        // Timer
      timer.start();

      std::size_t row( 0 );                               // Initialise row counter

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
        A( row, col( i, j, w ) ) = 1.;

        B[ row ] = - Q( i, j, w );
        ++row;

        // v_hzeta = 0
        A( row, col( i, j, v ) )       = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, v ) )   =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, v ) )   = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, v ) + 4. * Q( i + 1, j, v )
                                - Q( i + 2, j, v ) ) / ( 2 * dX ) );
        ++row;

        // s = 0
        A( row, col( 0, j, s ) ) = 1.;

        B[ row ] = - Q( i, j, s );
        ++row;

        // q_hzeta = 0
        A( row, col( i, j, q ) )      = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, q ) )  =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, q ) )  = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, q ) + 4. * Q( i + 1, j, q )
                                - Q( i + 2, j, q ) ) / ( 2 * dX ) );
        ++row;

        // Phi_hzeta = 0
        A( row, col( i, j, Phi ) )      = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, Phi ) )  =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, Phi ) )  = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, Phi ) + 4. * Q( i + 1, j, Phi )
                             - Q( i + 2, j, Phi ) ) / ( 2 * dX ) );
        ++row;

        // Psi = 0
        A( row, col( i, j, Psi ) )      =   1.;

        B[ row ]                        = - Q( i, j, Psi );
        ++row;

        // U_hzeta = 0
        A( row, col( i, j, U ) )        = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, U ) )    =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, U ) )    = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, U ) + 4. * Q( i + 1, j, U )
                             - Q( i + 2, j, U ) ) / ( 2 * dX ) );
        ++row;

        // Theta = 0
        A( row, col( i, j, Theta ) )    =   1.;

        B[ row ]                        = - Q( i, j, Theta );
        ++row;

      } // End of for loop over LHS eta nodes

      // Interior points between the hzeta boundaries
      for ( std::size_t i = 1; i < N; ++i )
      {
        // hzeta location
        double hzeta( HZETA_NODES[ i ] );
        double Xd( SSI.mesh_Xd( hzeta ) );
        double Xdd( SSI.mesh_Xdd( hzeta ) );

        // Wall transpiration
        double Phi_w( Phi_w_func( hzeta ) );

        // eta = 0 boundary ( bottom boundary )
        std::size_t j( 0 );
        double eta( ETA_NODES[ j ] );
        double Yd( SSI.mesh_Yd( eta ) );
        double Ydd( SSI.mesh_Ydd( eta ) );

        // v = 0
        A( row, col( i, j, v ) ) = 1.;

        B[ row ] = - Q( i, j, v );
        ++row;

        // w = 0
        A( row, col( i, j, w ) ) = 1.;

        B[ row ] = - Q( i, j, w );
        ++row;

        //s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
        A( row, col( i, j, s ) )     =   1.;
        A( row, col( i, j + 1, w ) ) = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 2, w ) ) =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 3, w ) ) =   iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s )
                   + iaR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) ) * Q( i, j + 1, w)
                   - iaR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) ) * Q( i, j + 2, w )
                   - iaR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, w );
        ++row;

        //q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
        A( row, col( i, j, q ) )     =   1.;
        A( row, col( i, j + 1, v ) ) = - 6. * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 2, v ) ) =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 3, v ) ) = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q )
                   + iaR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j + 1, v )
                   - iaR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 2, v )
                   + iaR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, v );
        ++row;

        // Phi = Phi_w
        A( row, col( i, j, Phi ) )        =   1.;

        B[ row ]                          = - Q( i, j, Phi ) + Phi_w;
        ++row;

        // Psi = 0
        A( row, col( i, j, Psi ) )        =   1.;

        B[ row ]                          = - Q( i, j, Psi );
        ++row;

        // U = 0
        A( row, col( i, j, U ) )          =   1;

        B[ row ]                          = - Q( i, j, U );
        ++row;

        // Theta = Psi_eta  - Phi_zeta
        A( row, col( i, j, Theta ) )   =   1;
        A( row, col( i, j, Psi ) )     =   3 * Yd / ( 2 * dY );
        A( row, col( i, j + 1, Psi ) ) = - 4 * Yd / ( 2 * dY );
        A( row, col( i, j + 2, Psi ) ) =   1 * Yd / ( 2 * dY );
        A( row, col( i + 1, j, Phi ) ) =   1 * Xd / ( 2 * dX );
        A( row, col( i - 1, j, Phi ) ) = - 1 * Xd / ( 2 * dX );

        B[ row ] = - Q( i, j, Theta )
                   + Yd * (( -3. * Q( i, j, Psi ) + 4. * Q( i, j + 1, Psi )
                            - Q( i, j + 2, Psi ) ) / ( 2 * dY ))
                   - Xd * ( Q( i + 1, j, Phi ) - Q( i - 1, j, Phi ) ) / ( 2 * dX );
        ++row;

        // Main interior grid points
        for ( std::size_t j = 1; j < M; ++j )
        {
          // eta location
          double eta( ETA_NODES[ j ] );
          double Yd( SSI.mesh_Yd( eta ) );
          double Ydd( SSI.mesh_Ydd( eta ) );
          Vector<double> Base( BASE_SOLUTION.get_interpolated_vars( eta ) );

          // Laplacian coefficients for finite-differencing

          // OS
          // X(i,j-1)
          double laplace_1_os =  ( Yd*Yd / (dY*dY) - Ydd / (2.*dY) ) ;
          // X(i-1,j)
          double laplace_3_os = ( Xd*Xd / (dX*dX) - Xdd / (2.*dX) ) / ( ZETA0 * ZETA0 );
          // X(i,j)
          double laplace_4_os = - 2 *( Yd*Yd / (dY*dY) + Xd*Xd / ( ZETA0 * ZETA0 * dX * dX ) )
                             - ALPHA * ALPHA;
          // X(i+1,j)
          double laplace_5_os = ( Xdd / (2.*dX) + Xd*Xd / (dX*dX) ) / ( ZETA0 * ZETA0 );
          // X(i,j+1)
          double laplace_7_os = ( Yd*Yd / (dY*dY) + Ydd / (2.*dY) );

          // Self-sim
          // X(i,j-1)
          double laplace_1 =  ( Yd*Yd/(dY*dY) - Ydd/ (2.*dY) ) ;
          // X(i-1,j)
          double laplace_3 = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) ) / ( ZETA0 * ZETA0 );
          // X(i,j)
          double laplace_4 = -2.*( Yd*Yd / (dY*dY) + Xd*Xd/( ZETA0 * ZETA0 * dX * dX ) );
          // X(i+1,j)
          double laplace_5 = ( Xdd/(2.*dX) + Xd*Xd/(dX*dX) ) / ( ZETA0 * ZETA0 );
          // X(i,j+1)
          double laplace_7 = ( Yd*Yd/(dY*dY) + Ydd/ (2.*dY) );

          // Guessed/known components and various derivative values
          Vector< std::complex<double> > Guess( Q.get_nodes_vars( i, j ) );
          Vector< std::complex<double> > Guess_eta( ( Q.get_nodes_vars( i, j + 1 )
                                    - Q.get_nodes_vars( i, j - 1 ) ) * ( Yd /( 2 * dY )) );
          Vector< std::complex<double> > Guess_hzeta( ( Q.get_nodes_vars( i + 1, j )
                                      - Q.get_nodes_vars( i - 1, j ) )
                                      * ( Xd /( 2 * dX )) );
          Vector< std::complex<double> > Guess_laplace( Q.get_nodes_vars( i, j - 1 ) * laplace_1
                                     +  Q.get_nodes_vars( i - 1, j ) * laplace_3
                                     +  Q.get_nodes_vars( i, j ) * laplace_4
                                     +  Q.get_nodes_vars( i + 1, j ) * laplace_5
                                     +  Q.get_nodes_vars( i, j + 1 ) * laplace_7 );

          Vector< std::complex<double> > Guess_eta_eta( ( Q.get_nodes_vars( i, j + 1 ) * laplace_7
                                                        + Q.get_nodes_vars( i, j ) * ( - 2. * Yd * Yd / ( dY * dY ) )
                                                        + Q.get_nodes_vars( i, j - 1 ) * laplace_1 ) );
          Vector< std::complex<double> > Guess_hzeta_hzeta( ( Q.get_nodes_vars( i + 1, j ) * laplace_5
                                                        + Q.get_nodes_vars( i, j ) * ( - 2. * Xd * Xd / ( ZETA0 * ZETA0 * dX * dX ))
                                                        + Q.get_nodes_vars( i - 1, j ) * laplace_3 ) );

          Vector< std::complex<double> > Guess_eta_hzeta( ( Q.get_nodes_vars( i + 1, j + 1 )
                                                          + Q.get_nodes_vars( i - 1, j - 1 )
                                                          - Q.get_nodes_vars( i + 1, j - 1 )
                                                          - Q.get_nodes_vars( i - 1, j + 1 ) ) * ( Xd * Yd / ( 4 * dX * dY) ) );

          double UBdd = BETA * ( Base[ UB ] * Base[ UB ] - 1. ) - Base[ PhiB ] * Base[ UBd ];

          ///////////////////////////////
          // 2D OrrSommerfeld equation //
          ///////////////////////////////

          ////////////////
          // v equation //
          ////////////////

          // ( i / alpha * Rx^(1/2) ) * Laplacian of v
          A( row, col( i, j - 1, v ) ) = iaR * laplace_1_os;
          A( row, col( i - 1, j, v ) ) = iaR * laplace_3_os;
          A( row, col( i, j, v ) )     = iaR * laplace_4_os;
          A( row, col( i + 1, j, v ) ) = iaR * laplace_5_os;
          A( row, col( i, j + 1, v ) ) = iaR * laplace_7_os;

          // + ( UB + UG - c_g ) * v
          A( row, col( i, j, v ) )    +=   Base[ UB ] + Guess[ U ] - C_GUESS;

          // + ((1 - beta) / Rx) * v
          A( row, col( i, j, v ) )    +=   ( 1.0 - BETA ) / RX;

          // + v_g * U
          A( row, col( i, j, U ) )     =   Q( i, j, v);

          // - v_g * c
          A( row, size - 1 )           = - Q( i, j, v );

          // - q
          A( row, col( i, j, q ) )     = - 1.;

          B[ row ] =  - iaR * ( Q( i, j - 1, v ) * laplace_1_os
                             +  Q( i - 1, j, v ) * laplace_3_os
                             +  Q( i, j, v )     * laplace_4_os
                             +  Q( i + 1, j, v ) * laplace_5_os
                             +  Q( i, j + 1, v ) * laplace_7_os )
                      - ( Base[ UB ] + Guess[ U ] - C_GUESS
                             + ( ( 1.0 - BETA ) / RX ) ) * Q( i, j, v )
                      + Q( i, j, q );
          ++row;

          ////////////////
          // w equation //
          ////////////////

          // ( i / alpha * Rx^(1/2) ) * Laplacian of v
          A( row, col( i, j - 1, w ) ) = iaR * laplace_1_os;
          A( row, col( i - 1, j, w ) ) = iaR * laplace_3_os;
          A( row, col( i, j, w ) )     = iaR * laplace_4_os;
          A( row, col( i + 1, j, w ) ) = iaR * laplace_5_os;
          A( row, col( i, j + 1, w ) ) = iaR * laplace_7_os;

          // + ( UB + UG - c_g ) * w
          A( row, col( i, j, w ) )    +=  Base[ UB ] + Guess[ U ] - C_GUESS;

          // + ((1 - beta) / Rx) * w
          A( row, col( i, j, w ) )    +=  ( 1. - BETA ) / RX;

          // + w_g * U
          A( row, col( i, j, U ) )     =   Q( i, j, w);

          // - w_g * c
          A( row, size - 1 )           = - Q( i, j, w );

          // - s
          A( row, col( i, j, s ) )     = - 1.;

          B[ row ] =  - iaR * ( Q( i, j - 1, w ) * laplace_1_os
                             +  Q( i - 1, j, w ) * laplace_3_os
                             +  Q( i, j, w )     * laplace_4_os
                             +  Q( i + 1, j, w ) * laplace_5_os
                             +  Q( i, j + 1, w ) * laplace_7_os )
                      - ( Base[ UB ] + Guess[ U ] - C_GUESS
                             + ( ( 1. - BETA ) / RX ) ) * Q( i, j, w )
                      + Q( i, j, s );
          ++row;

          ////////////////
          // q equation //
          ////////////////

          // q_{eta eta}
          A( row, col( i, j - 1, q ) )      =   Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ;
          A( row, col( i, j, q ) )          = - 2 * Yd * Yd / ( dY * dY );
          A( row, col( i, j + 1, q ) )      =   Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ;

          // - alpha^2 * q
          A( row, col( i, j, q ) )         += - ALPHA * ALPHA;

          // + ( 1 / zeta0 ) * s_{hzeta eta}
          A( row, col( i + 1, j + 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, s ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, s ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // + ( 1 / zeta0 ) * (UB' + UG_{eta}) * w_{hzeta}
          A( row, col( i + 1, j, w ) )      =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, w ) )      = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );

          // + ( 1 / zeta0 ) * w_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      =   Guess_hzeta[ w ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, U ) )      = - Guess_hzeta[ w ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * U_{hzeta} * w_{eta}
          A( row, col( i, j + 1, w ) )      = - Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, w ) )      =   Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * w_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      = - Guess_eta[ w ] * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, U ) )      =   Guess_eta[ w ] * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * UG_{hzeta eta} * w
          A( row, col( i, j, w ) )          = - Guess_eta_hzeta[ U ] / ZETA0;

          // - w_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  = - Guess[ w ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, U ) )  = - Guess[ w ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, U ) )  =   Guess[ w ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, U ) )  =   Guess[ w ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // - ( UB'' + UG_{eta eta} ) * v
          A( row, col( i, j, v ) )          = - ( UBdd + Guess_eta_eta[ U ] );

          // - v_g * U_{eta eta}
          A( row, col( i, j + 1, U ) )      = - Guess[ v ] * laplace_7;
          A( row, col( i, j, U ) )          =   Guess[ v ] * 2. * Yd * Yd / ( dY * dY );
          A( row, col( i, j - 1, U ) )      = - Guess[ v ] * laplace_1;

          B[ row ] = - ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) * Q( i, j - 1, q )
                     - ( - 2 * Yd * Yd / ( dY * dY ) ) * Q( i, j, q )
                     - ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) * Q( i, j + 1, q )
                     + ALPHA * ALPHA * Q( i, j, q )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j + 1, s )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j - 1, s )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j - 1, s )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j + 1, s )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i + 1, j, w )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i - 1, j, w )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j + 1, w )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j - 1, w )
                     + ( Guess_eta_hzeta[ U ] / ZETA0 ) * Q( i, j, w )
                     + ( UBdd + Guess_eta_eta[ U ] ) * Q( i, j, v );
          ++row;

          ////////////////
          // s equation //
          ////////////////

          // ( 1 / zeta0^2 ) * s_{hzeta hzeta}
          A( row, col( i - 1, j, s ) )      =   ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) )
                                              / ( ZETA0 * ZETA0 );
          A( row, col( i, j, s ) )          = - 2 * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
          A( row, col( i + 1, j, s ) )      =   ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) )
                                              / ( ZETA0 * ZETA0 );

          // - alpha^2 * s
          A( row, col( i, j, s ) )         += - ALPHA * ALPHA;

          // + ( 1 / zeta0 ) * q_{hzeta eta}
          A( row, col( i + 1, j + 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, q ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, q ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // + ( 1 / zeta0 ) * UG_{hzeta} * v_{eta}
          A( row, col( i, j + 1, v ) )     =   Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, v ) )     = - Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );

          // + ( 1 / zeta0 ) * v_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      =   Guess_eta[ v ] * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, U ) )      = - Guess_eta[ v ] * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * (UBd + UG_{eta}) * v_{hzeta}
          A( row, col( i + 1, j, v ) )     = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, v ) )     =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * v_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      = - Guess_hzeta[ v ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, U ) )      =   Guess_hzeta[ v ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * UG_{eta hzeta} * v
          A( row, col( i, j, v ) )         = - Guess_eta_hzeta[ U ] / ZETA0;

          // - v_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  = - Guess[ v ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, U ) )  = - Guess[ v ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, U ) )  =   Guess[ v ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, U ) )  =   Guess[ v ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // - ( 1 / zeta0^2 ) * UG_{hzeta hzeta} * w
          A( row, col( i, j, w ) )         = - Guess_hzeta_hzeta[ U ] / ( ZETA0 * ZETA0 );

          // - w_g * U_{hzeta hzeta} / zeta0^2
          A( row, col( i + 1, j, U ) )      = - Guess[ w ] * laplace_5;
          A( row, col( i, j, U ) )          =   Guess[ w ] * 2. * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
          A( row, col( i - 1, j, U ) )      = - Guess[ w ] * laplace_3;

          B[ row ] = - ( ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 1, j, s )
                     - ( - 2 * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i, j, s )
                     - ( ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i + 1, j, s )
                     + ALPHA * ALPHA * Q( i, j, s )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j + 1, q )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j - 1, q )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j - 1, q )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j + 1, q )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j + 1, v )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j - 1, v )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i + 1, j, v )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i - 1, j, v )
                     + ( Guess_eta_hzeta[ U ] / ZETA0 ) * Q( i, j, v )
                     + ( Guess_hzeta_hzeta[ U ] / ( ZETA0 * ZETA0 ) ) * Q( i, j, w );
                     //TODO replace these with Guess[ ... ] instead of Q( i, j, ...)
          ++row;

          ///////////////////////////////
          // Self-sim streak equations //
          ///////////////////////////////

          //////////////////
          // Phi equation //
          //////////////////

          // Laplacian of Phi
          A( row, col( i, j - 1, Phi ) )      = laplace_1;
          A( row, col( i - 1, j, Phi ) )      = laplace_3;
          A( row, col( i, j, Phi ) )          = laplace_4;
          A( row, col( i + 1, j, Phi ) )      = laplace_5;
          A( row, col( i, j + 1, Phi ) )      = laplace_7;
          // -(2-beta)*U_eta
          A( row, col( i, j + 1, U ) )        = - ( 2. - BETA ) * Yd / ( 2 * dY );
          A( row, col( i, j - 1, U ) )        =   ( 2. - BETA ) * Yd / ( 2 * dY );
          // Theta_hzeta
          A( row, col( i + 1, j, Theta ) )    =   Xd / ( 2 * dX );
          A( row, col( i - 1, j, Theta ) )    = - Xd / ( 2 * dX );

          // Residual
          B[ row ]      = - Guess_laplace[ Phi ] + ( 2. - BETA ) * Guess_eta[ U ]
                          - Guess_hzeta[ Theta ];

          ++row;

          //////////////////
          // Psi equation //
          //////////////////

          // Laplacian of Psi
          A( row, col( i, j - 1, Psi ) )      = laplace_1;
          A( row, col( i - 1, j, Psi ) )      = laplace_3;
          A( row, col( i, j, Psi ) )          = laplace_4;
          A( row, col( i + 1, j, Psi ) )      = laplace_5;
          A( row, col( i, j + 1, Psi ) )      = laplace_7;

          // -(2-beta)*U_hzeta / (zeta0^2)
          A( row, col( i + 1, j, U ) )        = - ( 2. - BETA ) * Xd
                                                / ( 2. * dX * ZETA0 * ZETA0 );
          A( row, col( i - 1, j, U ) )        =   ( 2. - BETA ) * Xd
                                                / ( 2. * dX * ZETA0 * ZETA0 );

          // -Theta_eta
          A( row, col( i, j + 1, Theta ) )    = - Yd / ( 2 * dY );
          A( row, col( i, j - 1, Theta ) )    =   Yd / ( 2 * dY );

          // Residual
          B[ row ]      = - Guess_laplace[ Psi ] + ( 2. - BETA )
                          * ( Guess_hzeta[ U ] )
                          / ( ZETA0 * ZETA0 )
                          + Guess_eta[ Theta ];

          ++row;

          ////////////////
          // U equation //
          ////////////////

          // Laplacian of U
          A( row, col( i, j - 1, U ) )        = laplace_1;
          A( row, col( i - 1, j, U ) )        = laplace_3;
          A( row, col( i, j, U ) )            = laplace_4;
          A( row, col( i + 1, j, U ) )        = laplace_5;
          A( row, col( i, j + 1, U ) )        = laplace_7;

          // -2 * beta * ( UB + UG ) * U
          A( row, col( i, j, U ) )           += - 2.* BETA * ( Base[ UB ]
                                                + Guess[ U ] );

          // ( hzeta * PsiB + PsiG ) * U_hzeta
          A( row, col( i + 1, j, U ) )       +=   ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                                * Xd / ( 2 * dX );
          A( row, col( i - 1, j, U ) )       += - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                                * Xd / ( 2 * dX );

          // [ PhiB + PhiG ] * U_eta
          A( row, col( i, j + 1, U ) )       +=   ( Base[ PhiB ] + Guess[ Phi ] )
                                                * Yd / ( 2 * dY );
          A( row, col( i, j - 1, U ) )       += - ( Base[ PhiB ] + Guess[ Phi ] )
                                                * Yd / ( 2 * dY );

          // [ UG_hzeta ] * Psi
          A( row, col( i, j, Psi ) )          =   Guess_hzeta[ U ];

          // ( UB' + UG_eta ) * Phi
          A( row, col( i, j, Phi ) )          =   Base[ UBd ] + Guess_eta[ U ];

          // Residual
          B[ row ]        = - Guess_laplace[ U ]
                            + BETA * ( 2. * Base[ UB ] + Guess[ U ] ) * Guess[ U ]
                            - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                            * ( Guess_hzeta[ U ] )
                            - ( Base[ PhiB ] + Guess[ Phi ] ) * Guess_eta[ U ]
                            - Base[ UBd ] * Guess[ Phi ];

          ++row;

          ////////////////////
          // Theta equation //
          ////////////////////

          // Laplacian of Theta
          A( row, col( i, j - 1, Theta ) )     = laplace_1;
          A( row, col( i - 1, j, Theta ) )     = laplace_3;
          A( row, col( i, j, Theta ) )         = laplace_4;
          A( row, col( i + 1, j, Theta ) )     = laplace_5;
          A( row, col( i, j + 1, Theta ) )     = laplace_7;

          // -2 * (1-beta) * (UB+UG) * [hzeta] * U_eta
          A( row, col( i, j + 1, U ) )         = - 2. * ( 1. - BETA )
                                                   * ( Base[ UB ] + Guess[ U ] )
                                                   * ( hzeta )
                                                   * Yd / ( 2 * dY );
          A( row, col( i, j - 1, U ) )         =   2. * ( 1. - BETA )
                                                   * ( Base[ UB ] + Guess[ U ] )
                                                   * ( hzeta )
                                                   * Yd / ( 2 * dY );

          // -2 * (1-beta) * (UB' + UG) * ( hzeta ) * U
          A( row, col( i, j, U ) )             = - 2. * ( 1. - BETA )
                                                   * ( Base[ UBd ] + Guess_eta[ U ] )
                                                   * ( hzeta );

          // (2 * (1-beta) * eta * UG_hzeta / (zeta0^2)) * U
          A( row, col( i, j, U ) )            +=  2. * ( 1. - BETA )
                                                  * eta * Guess_hzeta[ U ]
                                                  / ( ZETA0 * ZETA0 );

          // 2 * (1-beta) * eta * (UB + UG) * U_hzeta / ( zeta0^2 )
          A( row, col( i + 1, j, U ) )         =  2. * ( 1. - BETA ) * eta
                                                  * ( Base[ UB ] + Guess[ U ] )
                                                  * Xd / ( 2 * dX * ZETA0 * ZETA0 );
          A( row, col( i - 1, j, U ) )         = -2. * ( 1. - BETA ) * eta
                                                  * ( Base[ UB ] + Guess[ U ] )
                                                  * Xd / ( 2 * dX * ZETA0 * ZETA0 );

          // ( PhiB + PhiG ) * Theta_eta
          A( row, col( i, j + 1, Theta ) )    +=  ( Base[ PhiB ] + Guess[ Phi ] ) * Yd
                                                  / ( 2 * dY );
          A( row, col( i, j - 1, Theta ) )    += -( Base[ PhiB ] + Guess[ Phi ] ) * Yd
                                                  / ( 2 * dY );

          // (hzeta * ThetaB' + ThetaG_eta ) * Phi
          A( row, col( i, j, Phi ) )           =   hzeta * Base[ ThetaBd ]
                                                 + Guess_eta[ Theta ];

          // (hzeta * PsiB + PsiG ) * Theta_hzeta
          A( row, col( i + 1, j, Theta ) )    +=  ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                                 * Xd / ( 2 * dX );
          A( row, col( i - 1, j, Theta ) )    += -( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                                 * Xd / ( 2 * dX );

          // [ThetaB + ThetaG_hzeta] * Psi
          A( row, col( i, j, Psi ) )           =  Base[ ThetaB ] + Guess_hzeta[ Theta ];

          // (2-beta) * ( UB + UG ) * Theta
          A( row, col( i, j, Theta ) )        +=   ( 2. - BETA ) * ( Base[ UB ]
                                                 + Guess[ U ] );

          // (2-beta) * ( hzeta * ThetaB + ThetaG ) * U
          A( row, col( i, j, U ) )            +=   ( 2. - BETA ) * ( hzeta * Base[ ThetaB ]
                                                 + Guess[ Theta ] );

          // 4 * real( v_g_zeta - w_g_eta ) * v_eta

          // 4 * real( v_g_zeta - w_g_eta ) * w_zeta

          // 4 * real( v_g_eta + w_g_zeta ) * v_zeta

          // - 4 * real( v_g_eta + w_g_zeta ) * w_eta

          // 2 * real( 2 * v_g_eta_zeta + w_g_zeta_zeta - w_g_eta_eta ) * v

          // 4 * real( v_g ) * v_eta_zeta

          // 2 * real( v_g ) * w_zeta_zeta

          // - 2 * real( v_g ) * w_eta_eta

          // 2 * real( v_g_zeta_zeta - v_g_eta_eta - 2 * w_g_zeta_eta ) * w

          // 2 * real( w_g ) * v_zeta_zeta

          // - 2 * real( w_g ) * v_eta_eta

          // - 4 * real( w_g ) * w_zeta_eta

          // Forcing (RHS)
          std::complex<double> F_2;
          F_2 =   2. * ( Guess_hzeta[ v ] - Guess_eta[ w ] )
                     * ( conj( Guess_eta[ v ] ) + conj( Guess_hzeta[ w ] ) )
                + Guess[ v ] * ( 2. * conj( Guess_eta_hzeta[ v ] ) + conj( Guess_hzeta_hzeta[ w ] ) - conj( Guess_eta_eta[ w ] ) )
                + Guess[ w ] * ( conj( Guess_hzeta_hzeta[ v ] ) - conj( Guess_eta_eta[ v ] ) - 2. * conj( Guess_eta_hzeta[ w ] ) );
          F_2 = F_2 + conj( F_2 );

          // Residual
          B[ row ]      = - Guess_laplace[ Theta ]
                          + 2. * ( 1. - BETA )
                          * ( hzeta * ( Base[ UB ] + Guess[ U ] )
                          * Guess_eta[ U ] + hzeta * Base[ UBd ] * Guess[ U ]
                          - eta * ( Base[ UB ] + Guess[ U ] )
                          * ( Guess_hzeta[ U ] )
                          / ( ZETA0 * ZETA0 ) )
                          - ( Base[ PhiB ] + Guess[ Phi ] ) * Guess_eta[ Theta ]
                          - hzeta * Base[ ThetaBd ] * Guess[ Phi ]
                          - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                          * ( Guess_hzeta[ Theta ] ) - Guess[ Psi ] * Base[ ThetaB ]
                          - ( 2. - BETA ) * ( ( Base[ UB ] + Guess[ U ] )
                          * Guess[ Theta ] + hzeta * Base[ ThetaB ] * Guess[ U ] );
                          //+ std::pow( RX, -1.0/6.0 ) * SIGMA * SIGMA * F_2;
          ++row;

        }

        // eta = eta_inf boundary ( top boundary )
        j = M ;
        eta = ETA_NODES[ j ];
        Yd = SSI.mesh_Yd( eta );

        // v = 0
        A( row, col( i, j, v ) ) = 1;

        B[ row ] = - Q( i, j, v );
        ++row;

        // w = 0
        A( row, col( i, j, w ) ) = 1;

        B[ row ] = - Q( i, j, w );
        ++row;

        // s - (i / (alpha*Rx^(1/2))) * w_{eta eta} = 0
        A( row, col( i, j, s ) )        =   1;
        A( row, col( i, j - 1, w ) )    = - iaR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j - 2, w ) )    =   iaR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j - 3, w ) )    =   iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s )
                   + iaR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) ) * Q( i, j - 1, w)
                   - iaR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) ) * Q( i, j - 2, w )
                   - iaR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, w );
        ++row;

        // q - (i / (alpha*Rx^(1/2))) * v_{eta eta} = 0
        A( row, col( i, j, q ) )        =   1;
        A( row, col( i, j - 1, v ) )    = - 6. * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j - 2, v ) )    =   ( 3. / 2. ) * iaR * Yd * Yd / ( dY * dY );
        A( row, col( i, j - 3, v ) )    = - ( 2. / 9. ) * iaR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q )
                   + iaR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j - 1, v )
                   - iaR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 2, v )
                   + iaR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, v );
        ++row;

        const double rsq =  eta * eta + ZETA0 * ZETA0 * hzeta * hzeta ;
        // Phi_eta*( eta^2 + zeta_0^2*hzeta^2) + [ 2*eta - (eta^2 + zeta_0^2*hzeta^2)/eta ]*Phi = 0
        A( row, col( i, j, Phi ) )        =   3 * Yd * rsq / ( 2 * dY );
        A( row, col( i, j - 1, Phi ) )    = - 4 * Yd * rsq / ( 2 * dY );
        A( row, col( i, j - 2, Phi ) )    =   1 * Yd * rsq / ( 2 * dY );
        A( row, col( i, j, Phi ) )       +=   2 * eta - ( rsq / eta );

        B[ row ] = - ( 3. * Q( i, j, Phi ) - 4. * Q( i, j-1, Phi )
                       + Q( i, j-2, Phi ) ) * Yd * rsq / ( 2 * dY )
                   + ( ( rsq / eta ) - 2. * eta ) * Q( i, j, Phi );
        ++row;

        // Psi_eta*( eta^2 + zeta_0^2*hzeta^2) + 2 * eta * Psi = 0
        A( row, col( i, j, Psi ) )        =   3 * Yd * rsq / (2*dY);
        A( row, col( i, j - 1, Psi ) )    = - 4 * Yd * rsq / (2*dY);
        A( row, col( i, j - 2, Psi ) )    =   1 * Yd * rsq / (2*dY);
        A( row, col( i, j, Psi ) )       +=   2 * eta;

        B[ row ] =  - (( 3. * Q( i, j, Psi ) - 4. * Q( i, j-1, Psi )
                        + Q( i, j-2, Psi ) ) * Yd * rsq / ( 2 * dY ))
                        - 2. * eta * Q( i, j, Psi );
        ++row;

        // U = 0
        A( row, col( i, j, U ) )            =   1;

        B[ row ]                            = - Q( i, j, U );
        ++row;

        // Theta = 0
        A( row, col( i, j, Theta ) )        =   1;

        B[ row ]                            = - Q( i, j, Theta );
        ++row;

      } // End of loop over interior nodes

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

        const double rsq =  eta * eta + ZETA0 * ZETA0 * hzeta * hzeta ;

        // w = 0
        A( row, col( i, j, w ) ) = 1;

        B[ row ] = - Q( i, j, w );
        ++row;

        // v = 0
        A( row, col( i, j, v ) ) = 1;

        B[ row ] = - Q( i, j, v );
        ++row;

        // q - (1/zeta0^2) * (i / (alpha*Rx^(1/2))) * v_{hzeta hzeta} = 0
        A( row, col( i, j, q ) )        =   1;
        A( row, col( i - 1, j, v ) )    = - iaR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) )
                                      / ( ZETA0 * ZETA0 );
        A( row, col( i - 2, j, v ) )    =   iaR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) )
                                      / ( ZETA0 * ZETA0 );
        A( row, col( i - 3, j, v ) )    =   iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );

        B[ row ] = - Q( i, j, q )
                   + ( iaR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 1, j, v )
                   - ( iaR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 2, j, v )
                   - ( iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 3, j, v );
        ++row;

        // s - (1/zeta0^2) * (i / (alpha*Rx^(1/2))) * w_{hzeta hzeta} = 0
        A( row, col( i, j, s ) )        =   1;
        A( row, col( i - 1, j, w ) )    = - 6. * iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
        A( row, col( i - 2, j, w ) )    =   ( 3. / 2. ) * iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
        A( row, col( i - 3, j, w ) )    = - ( 2. / 9. ) * iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );

        B[ row ] = - Q( i, j, s )
                   + ( 6. * iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 1, j, w )
                   - ( ( 3. / 2. ) * iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 2, j, w )
                   + ( ( 2. / 9. ) * iaR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 3, j, w );
        ++row;

        // (eta^2 + zeta_0^2 * hzeta^2) * Phi_hzeta + 2 * zeta_0^2 * hzeta * Phi = 0
        A( row, col( i, j, Phi ) )          =   3 * Xd * rsq / ( 2 * dX );
        A( row, col( i - 1, j, Phi ) )      = - 4 * Xd * rsq / ( 2 * dX );
        A( row, col( i - 2, j, Phi ) )      =   1 * Xd * rsq / ( 2 * dX );
        A( row, col( i, j, Phi ) )         +=   2 * ZETA0 * ZETA0 * hzeta;

        B[ row ]        = - rsq * ( 3. * Q( i, j, Phi) - 4. * Q( i - 1, j, Phi)
                          + Q( i - 2, j, Phi) ) * Xd / ( 2 * dX )
                          - 2. * ZETA0 * ZETA0 * hzeta * Q( i, j, Phi );
        ++row;

        // (eta^2 + zeta_0^2 * hzeta^2)*Psi_hzeta + (2*zeta_0^2*hzeta-(eta^2 + zeta_0^2*hzeta^2)/hzeta)*Psi = 0
        A( row, col( i, j, Psi ) )          =   3 * Xd * rsq / ( 2 * dX );
        A( row, col( i - 1, j, Psi ) )      = - 4 * Xd * rsq / ( 2 * dX );
        A( row, col( i - 2, j, Psi ) )      =   1 * Xd * rsq / ( 2 * dX );
        A( row, col( i, j, Psi ) )         +=   2 * ZETA0 * ZETA0 * hzeta - (rsq / hzeta);


        B[ row ]  = - (rsq * ( 3. * Q( i, j, Psi ) - 4. * Q( i - 1, j, Psi )
                      + Q( i - 2, j, Psi) ) * Xd / ( 2 * dX ))
                      - 2. * ZETA0 * ZETA0 * hzeta  * Q( i, j, Psi)
                      + (rsq / hzeta)  * Q( i, j, Psi);
        ++row;

        // hzeta * U_hzeta + 2 * U = 0
        A( row, col( i, j, U ) )            =   hzeta * 3. * Xd / ( 2 * dX ) + 2.;
        A( row, col( i - 1, j, U ) )        = - hzeta * 4. * Xd / ( 2 * dX );
        A( row, col( i - 2, j, U ) )        =   hzeta * 1. * Xd / ( 2 * dX );

        B[ row  ] = - ( hzeta * Xd * ( 3. * Q( i, j, U ) - 4. * Q( i - 1, j, U )
                      + Q( i - 2, j, U) ) / ( 2 * dX ) ) - 2. * Q( i, j, U );
        ++row;

        // hzeta * Theta_hzeta + Theta = 0
        A( row, col( i, j, Theta )  )       =   hzeta * 3. * Xd / ( 2 * dX ) + 1.;
        A( row, col( i - 1, j, Theta ) )    = - hzeta * 4. * Xd / ( 2 * dX );
        A( row, col( i - 2, j, Theta ) )    =   hzeta * 1. * Xd / ( 2 * dX );

        B[ row ] = - ( hzeta * Xd * ( 3. * Q( i, j, Theta ) - 4. * Q( i - 1, j, Theta )
                          + Q( i - 2, j, Theta ) ) / ( 2 * dX ) ) - Q( i, j, Theta ) ;
        ++row;

      } // End of loop over RHS eta nodes

      // q_eta(0,0) = 1 (extra condition)
      double eta( ETA_NODES[ 0 ] );
      double Yd( SSI.mesh_Yd( eta ) );
      double Ydd( SSI.mesh_Ydd( eta ) );
      A( row, col( 0, 0, q ) ) = - 3 * Yd / ( 2 * dY );
      A( row, col( 0, 1, q ) ) =   4 * Yd / ( 2 * dY );
      A( row, col( 0, 2, q ) ) = - 1 * Yd / ( 2 * dY );
      B[ row ] = 1.0 + ( 3 * Yd / ( 2 * dY ) ) * Q( 0, 0, q )
                     - ( 4 * Yd / ( 2 * dY ) ) * Q( 0, 1, q )
                     + ( 1 * Yd / ( 2 * dY ) ) * Q( 0, 2, q );
      ++row;

      max_residual = B.norm_inf();
      std::cout << "***                                              Maximum residual = "
                << max_residual << std::endl;

      // Solve the sparse system
      Vector< std::complex<double> > x;
      x = A.solve( B );
      B = x;

      timer.print();
      timer.stop();

      // Update the known values using the correction which we just found
      for ( std::size_t i = 0; i < N + 1; ++i )
      {
        for ( std::size_t j = 0; j < M + 1; ++j )
        {
          Q( i, j, v )     += B[ col( i, j, v ) ];
          Q( i, j, w )     += B[ col( i, j, w ) ];
          Q( i, j, q )     += B[ col( i, j, q ) ];
          Q( i, j, s )     += B[ col( i, j, s ) ];
          Q( i, j, Phi )   += B[ col( i, j, Phi ) ];
          Q( i, j, Psi )   += B[ col( i, j, Psi ) ];
          Q( i, j, U )     += B[ col( i, j, U ) ];
          Q( i, j, Theta ) += B[ col( i, j, Theta ) ];
        }
      }

      C_GUESS += B[ size - 1 ];

      std::cout << "***    Iteration = " << iteration
                << "    Maximum correction = " << B.norm_inf() << std::endl;

      ++iteration;
    }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) );

    if ( iteration >= max_iterations )
    {
      std::string problem;
      problem = "STOPPED AFTER TOO MANY ITERATIONS\n";
      throw Error( problem );
    }

    // Push the data back into the unmapped domain
    for ( std::size_t i = 0; i < N + 1; ++i )
    {
      double hzeta=HZETA_NODES[i];
      for ( std::size_t j = 0; j < M + 1; ++j )
      {
        double eta=ETA_NODES[j];
        // First 4 values are the wave
        Q_OUTPUT( i, j, 0 ) = Q( i, j, v);
        Q_OUTPUT( i, j, 1 ) = Q( i, j, w);
        Q_OUTPUT( i, j, 2 ) = Q( i, j, q);
        Q_OUTPUT( i, j, 3 ) = Q( i, j, s);
        // next 4 values output are the streak without the underlying base flow
        Q_OUTPUT( i, j, 4 ) = Q( i, j, Phi);
        Q_OUTPUT( i, j, 5 ) = Q( i, j, Psi);
        Q_OUTPUT( i, j, 6 ) = Q( i, j, U);
        Q_OUTPUT( i, j, 7 ) = Q( i, j, Theta);
        // final 4 values are the "full" solution, but still with the zeta0 scaling
        Q_OUTPUT( i, j, 8 ) =   Q( i, j, Phi)
                            + BASE_SOLUTION.get_interpolated_vars( eta )[PhiB];
        Q_OUTPUT( i, j, 9 ) = Q( i, j, Psi)
                            + hzeta * BASE_SOLUTION.get_interpolated_vars( eta )[PsiB];
        Q_OUTPUT( i, j, 10 ) =   Q( i, j, U)
                            + BASE_SOLUTION.get_interpolated_vars( eta )[UB];
        Q_OUTPUT( i, j, 11 ) = Q( i, j, Theta)
                            + hzeta * BASE_SOLUTION.get_interpolated_vars( eta )[ThetaB];
      }
    }

  } // End of solve_local method

} // End of namespace TSL

#endif
