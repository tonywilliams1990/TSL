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
#include <fstream>

#include "Vector.h"
#include "SparseMatrix.h"
#include "SelfSimInjection.h"
#include "SparseEigenSystem.h"
#include "TwoD_node_mesh.h"
#include "Error.h"

#include "Eigenvalue"

enum class VWI_base{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };
enum class VWI_enum{ v_r, v_i, w_r, w_i, q_r, q_i, s_r, s_i, Phi, Psi, U, Theta };

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
      bool SPEED_UP;                 // Reuse the factorised matrix for speed gains
      std::string OUTPUT_PATH;       // Output path string
      bool SOLVED;                   // True if equations have been solved

      // Mesh
      double HZETA_RIGHT;         // Size of the domain in the zeta_hat direction
      double ETA_TOP;             // Size of the domain in the eta direction
      std::size_t N;                 // Number of intervals in the zeta_hat direction
      std::size_t M;                 // Number of intervals in the eta direction
      Vector<double> ETA_NODES;
      Vector<double> HZETA_NODES;
      Vector<double> X_NODES;
      Vector<double> Y_NODES;
      Vector<double> BASE_ETA_NODES;

      // Solution
      OneD_node_mesh<double> BASE_SOLUTION;            // Base flow ODE solution
      TwoD_node_mesh<double> Q;        // Current guess mesh
      TwoD_node_mesh<double> Q_OUTPUT; // Output mesh

      std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
      {
        // Return the column number for the kth variable at node (i,j)
        return 12 * ( i * ( SSI.eta_intervals() + 1 ) + j ) + k;
      }

    public:

      /// Constructor
      VWI( SelfSimInjection& ssi, double& alpha, double& rx, double& sigma )
      {
        SSI       = ssi;
        ALPHA     = alpha;
        RX        = rx;
        SIGMA     = sigma;
        BETA      = SSI.hartree();
        ZETA0     = SSI.injection_width();
        K         = SSI.injection();
        SPEED_UP  = false;

        HZETA_RIGHT    = SSI.hzeta_right();
        ETA_TOP        = SSI.eta_top();
        N              = SSI.hzeta_intervals();
        M              = SSI.eta_intervals();
        ETA_NODES      = SSI.eta_nodes();
        HZETA_NODES    = SSI.hzeta_nodes();
        X_NODES        = SSI.x_nodes();
        Y_NODES        = SSI.y_nodes();
        BASE_ETA_NODES = SSI.base_eta_nodes();


        BASE_SOLUTION = SSI.base_flow_solution();
        TwoD_node_mesh<double> q( X_NODES, Y_NODES, 12 );
        TwoD_node_mesh<double> q_output( HZETA_NODES, ETA_NODES, 16 );
        Q        = q;
        Q_OUTPUT = q_output;

        OUTPUT_PATH = SSI.output_path();

      }

      /// Destructor
	   	~VWI()
      {
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

      /// Return a pointer to the wave amplitude
      double& Sigma()
      {
        return SIGMA;
      }

      /// Return a pointer to the injection velocity
      double& injection()
      {
        return K;
      }

      /// Return a pointer to current eigenvalue guess
      std::complex<double>& c_guess()
      {
        return C_GUESS;
      }

      /// Set an initial guess
      void set_guess( const TwoD_node_mesh<double>& guess )
      {
        //TODO error checking for size of mesh and number of variables
        //std::cout << "X_NODES.size() = " << X_NODES.size() << std::endl;
        for ( std::size_t i = 0; i < X_NODES.size(); ++i )
        {
          for ( std::size_t j = 0; j < Y_NODES.size(); ++j )
          {
            for ( std::size_t n = 0; n < 12; ++n )
            {
              Q( i, j, n )  = guess( i, j, n );
            }
          }
        }

        // Normalise based on the extra condition q_eta(0,0) = 1e-5
        int q_r   = static_cast<int>(VWI_enum::q_r);
        int q_i   = static_cast<int>(VWI_enum::q_i);
        const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
        double lambda_r, lambda_i;
        double eta( ETA_NODES[ 0 ] );
        double Yd( SSI.mesh_Yd( eta ) );
        lambda_r = ( 3 * Yd / ( 2 * dY ) ) * Q( 0, 0, q_r )
                 - ( 4 * Yd / ( 2 * dY ) ) * Q( 0, 1, q_r )
                 + ( 1 * Yd / ( 2 * dY ) ) * Q( 0, 2, q_r );
        lambda_i = ( 3 * Yd / ( 2 * dY ) ) * Q( 0, 0, q_i )
                 - ( 4 * Yd / ( 2 * dY ) ) * Q( 0, 1, q_i )
                 + ( 1 * Yd / ( 2 * dY ) ) * Q( 0, 2, q_i );

        lambda_r /= 1e-5;
        double lambda( sqrt( lambda_r * lambda_r + lambda_i * lambda_i ) );
        //lambda /= 1e-5;

        for ( std::size_t i = 0; i < N + 1; ++i )
        {
          for ( std::size_t j = 0; j < M + 1; ++j )
          {
            for ( std::size_t n = 0; n < 8; ++n )
            {
              Q( i, j, n ) = Q( i, j, n ) / lambda;
            }
          }
        }

        // Normalise based on the extra condition int v = 1 along centreline
        /*int v_r   = static_cast<int>(VWI_enum::v_r);
        int v_i   = static_cast<int>(VWI_enum::v_i);
        double lambda_r, lambda_i;
        const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
        double eta, Yd;

        // Sum
        for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
        {
          //dx = ( NODES[ node + 1 ] - NODES[ node ] );
          eta = ETA_NODES[ j ];
          Yd = SSI.mesh_Yd( eta );
          lambda_r += 0.5 * dY * ( Q( 0, j, v_r ) + Q( 0, j + 1, v_r ) ) / Yd;
          lambda_i += 0.5 * dY * ( Q( 0, j, v_i ) + Q( 0, j + 1, v_i ) ) / Yd;
        }

        double lambda( sqrt( lambda_r * lambda_r + lambda_i * lambda_i ) );

        for ( std::size_t i = 0; i < N + 1; ++i )
        {
          for ( std::size_t j = 0; j < M + 1; ++j )
          {
            for ( std::size_t n = 0; n < 8; ++n )
            {
              Q( i, j, n )  = Q( i, j, n ) / lambda;
            }
          }
        }*/

      }

      /// Return a pointer to the current guess mesh
      TwoD_node_mesh<double>& current_guess()
      {
        return Q;
      }

      /// Return the solution mesh
      TwoD_node_mesh<double> solution() {
        return Q_OUTPUT;
      }

      /// Reset the solution mesh from an external mesh
      void set_solution( TwoD_node_mesh<double> sol ){
        if ( SOLVED ){
          Q_OUTPUT = sol;
          for ( std::size_t i = 0; i < X_NODES.size(); ++i )
          {
            for ( std::size_t j = 0; j < Y_NODES.size(); ++j )
            {
              for ( std::size_t n = 0; n < 12; ++n )
              {
                Q( i, j, n )  = Q_OUTPUT( i, j, n );
              }
            }
          }
        }
        else { throw Error( "set_solution() error equations have not been solved." ); }
      }

      /// Set solved bool
      void set_solved( bool solved ){ SOLVED = solved; }

      /// Solve the sparse eigenvalue problem using Newton iteration
      void solve_local();

      /// Virtual function for defining the transpiration function
      virtual double Phi_w_func( const double& hzeta )
      {
        throw Error( "--- Phi_w function is not defined ---" );
      }

      /// Use a factorised matrix to solve the equations more quickly (but less accurately)
      bool& speed_up()
      {
        return SPEED_UP;
      }

      /// Set the output path
      void set_output_path()
      {
        std::ostringstream ss;
        ss << "./DATA/VWI/";
        OUTPUT_PATH = ss.str();
      }

      /// Make the output directory
      void make_output_directory()
      {
        int status = mkdir( OUTPUT_PATH.c_str(), S_IRWXU );
        if ( status == 0 ) {
        std::cout << "  * Output directory " + OUTPUT_PATH +
                " has been made successfully." << std::endl;
        }
      }

      /// Make eigenvalues directory
      void make_eigenvalues_directory()
      {
        OUTPUT_PATH += "eigenvalues/";
        int status = mkdir( OUTPUT_PATH.c_str(), S_IRWXU );
        if ( status == 0 ) {
        std::cout << "  * Eigenvalue directory " + OUTPUT_PATH +
                " has been made successfully." << std::endl;
        }
        set_output_path();
      }

      /// Output the solution mesh
      void output(){
        if ( SOLVED ) {
          // Convert param to string
          std::stringstream ss;
          ss << "K_" << K << "_R_" << RX << "_Sigma_" << SIGMA << "_" << N + 1
             << "x" << M + 1 << "_" << HZETA_RIGHT << "_" << ETA_TOP << ".dat";
          std::string str = ss.str();
          Q_OUTPUT.dump_gnu( OUTPUT_PATH + str );
        }
        else { throw Error( "output() error equations have not been solved." ); }
      }

      /// Output the eigenvalue
      void output_eigenvalue(){
        if ( SOLVED ) {
          // Convert param to string
          std::stringstream ss;
          ss << "eval_K_" << K << "_R_" << RX << "_Sigma_" << SIGMA << "_" << N + 1
             << "x" << M + 1 << "_" << HZETA_RIGHT << "_" << ETA_TOP << ".dat";
          std::string str = ss.str();
          Vector<double> c( 2, 0.0 );
          c[0] = real( C_GUESS );
          c[1] = imag( C_GUESS );
          c.output( OUTPUT_PATH + "eigenvalues/" + str, 12 );
        }
        else { throw Error( "output_eigenvalue() error equations have not been solved." ); }
      }

      void solve_check_exists()
      {
        TwoD_node_mesh<double> sol( HZETA_NODES, ETA_NODES, 16 );
        std::stringstream ss;
        ss << "K_" << K << "_R_" << RX << "_Sigma_" << SIGMA << "_" << N + 1
           << "x" << M + 1 << "_" << HZETA_RIGHT << "_" << ETA_TOP << ".dat";
        std::string str = ss.str();
        // Don't bother solving it all again if the solution file already exists
        bool exists, eval_exists;
        exists = Utility::file_exists( OUTPUT_PATH + str );
        eval_exists = Utility::file_exists( OUTPUT_PATH + "eigenvalues/eval_" + str );

        try
        {
          if ( !exists || !eval_exists ){
            solve_local();
            output();
            output_eigenvalue();
           }
          if ( exists ){
            std::cout << "--- Reading solution from file" << std::endl;
            sol.read( OUTPUT_PATH + str );
            std::cout << "--- Finished reading" << std::endl;
            set_solved( true );
            set_solution( sol );
          }
          if ( eval_exists ){
            std::ifstream infile(OUTPUT_PATH + "eigenvalues/eval_" + str);
            double a;
            int i( 0 );
            while (infile >> a )
            {
              //std::cout << "a = " << a << std::endl;
              if ( i== 0 ){ C_GUESS.real( a ); }
              if ( i== 1 ){ C_GUESS.imag( a ); }
              ++i;
            }
          }
        }
        catch ( std::runtime_error )
        {
          std::cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED (SSI) \033[0m\n";
          assert( false );
        }

      }

      void refine_mesh( std::size_t& n, std::size_t& m, std::size_t& mb  )
      {
        TwoD_node_mesh<double> solution( Q_OUTPUT );
        SSI.hzeta_intervals() = n;
        SSI.eta_intervals() = m;
        SSI.base_intervals() = mb;
        SSI.mesh_setup();
        SSI.solve_base_flow();
        BASE_SOLUTION  = SSI.base_flow_solution();
        N              = SSI.hzeta_intervals();
        M              = SSI.eta_intervals();
        ETA_NODES      = SSI.eta_nodes();
        HZETA_NODES    = SSI.hzeta_nodes();
        X_NODES        = SSI.x_nodes();
        Y_NODES        = SSI.y_nodes();
        BASE_ETA_NODES = SSI.base_eta_nodes();

        Q_OUTPUT.remesh1( HZETA_NODES, ETA_NODES );
        Q.remesh1( X_NODES, Y_NODES );
        set_solution( Q_OUTPUT );
        solve_local();
      }

  }; // End of class OrrSommerfeld_2D

  void VWI::solve_local()
  {
    //SlepcInitialize( NULL, NULL, (char*)0, (char*)0 );
    int v_r   = static_cast<int>(VWI_enum::v_r);
    int w_r   = static_cast<int>(VWI_enum::w_r);
    int q_r   = static_cast<int>(VWI_enum::q_r);
    int s_r   = static_cast<int>(VWI_enum::s_r);
    int v_i   = static_cast<int>(VWI_enum::v_i);
    int w_i   = static_cast<int>(VWI_enum::w_i);
    int q_i   = static_cast<int>(VWI_enum::q_i);
    int s_i   = static_cast<int>(VWI_enum::s_i);
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
    std::size_t size( 12 * N_eta * N_hzeta + 2 );
    Vector<double> B( size, 0.0 );;

    // Step sizes
    const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
    const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );

    // Iterate to a solution
    double max_residual( 0.0 );             // Maximum residual
    std::size_t iteration( 0 );             // Initialise iteration counter
    std::size_t max_iterations( 30 );       // Maximum number of iterations

    double aR ( 1.0 / ( ALPHA * sqrt( RX ) ) );
    double Rsigma2( std::pow( RX, -1.0/6.0 ) * SIGMA * SIGMA );

    if( SPEED_UP ){ max_iterations = 500; }
    // Eigen objects (only used if SPEED_UP=true)
    Eigen::SparseMatrix<double, Eigen::ColMajor, long long> A_Eigen( size, size );
    Eigen::SparseLU< Eigen::SparseMatrix<double, Eigen::ColMajor, long long> > solver;

    do{
      SparseMatrix<double> A( size, size );
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

        // w_r = 0
        A( row, col( i, j, w_r ) ) = 1.;

        B[ row ] = - Q( i, j, w_r );
        ++row;

        // w_i = 0
        A( row, col( i, j, w_i ) ) = 1.;

        B[ row ] = - Q( i, j, w_i );
        ++row;

        // v_r_hzeta = 0
        A( row, col( i, j, v_r ) )       = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, v_r ) )   =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, v_r ) )   = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, v_r ) + 4. * Q( i + 1, j, v_r )
                                - Q( i + 2, j, v_r ) ) / ( 2 * dX ) );
        ++row;

        // v_i_hzeta = 0
        A( row, col( i, j, v_i ) )       = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, v_i ) )   =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, v_i ) )   = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, v_i ) + 4. * Q( i + 1, j, v_i )
                                - Q( i + 2, j, v_i ) ) / ( 2 * dX ) );
        ++row;

        // s_r = 0
        A( row, col( 0, j, s_r ) ) = 1.;

        B[ row ] = - Q( i, j, s_r );
        ++row;

        // s_i = 0
        A( row, col( 0, j, s_i ) ) = 1.;

        B[ row ] = - Q( i, j, s_i );
        ++row;

        // q_r_hzeta = 0
        A( row, col( i, j, q_r ) )      = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, q_r ) )  =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, q_r ) )  = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, q_r ) + 4. * Q( i + 1, j, q_r )
                                - Q( i + 2, j, q_r ) ) / ( 2 * dX ) );
        ++row;

        // q_i_hzeta = 0
        A( row, col( i, j, q_i ) )      = - 3 * Xd / ( 2 * dX );
        A( row, col( i + 1, j, q_i ) )  =   4 * Xd / ( 2 * dX );
        A( row, col( i + 2, j, q_i ) )  = - 1 * Xd / ( 2 * dX );

        B[ row ] = -( Xd * ( - 3. * Q( i, j, q_i ) + 4. * Q( i + 1, j, q_i )
                                - Q( i + 2, j, q_i ) ) / ( 2 * dX ) );
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

        // v_r = 0
        A( row, col( i, j, v_r ) ) = 1.;

        B[ row ] = - Q( i, j, v_r );
        ++row;

        // v_i = 0
        A( row, col( i, j, v_i ) ) = 1.;

        B[ row ] = - Q( i, j, v_i );
        ++row;

        // w_r = 0
        A( row, col( i, j, w_r ) ) = 1.;

        B[ row ] = - Q( i, j, w_r );
        ++row;

        // w_i = 0
        A( row, col( i, j, w_i ) ) = 1.;

        B[ row ] = - Q( i, j, w_i );
        ++row;

        //s_r + (1 / (alpha*Rx^(1/2))) * w_i_{eta eta} = 0
        A( row, col( i, j, s_r ) )     =   1.;
        A( row, col( i, j + 1, w_i ) ) =   aR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 2, w_i ) ) = - aR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 3, w_i ) ) = - aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s_r )
                   - aR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) ) * Q( i, j + 1, w_i )
                   + aR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) ) * Q( i, j + 2, w_i )
                   + aR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, w_i );
        ++row;

        //s_i - (1 / (alpha*Rx^(1/2))) * w_r_{eta eta} = 0
        A( row, col( i, j, s_i ) )     =   1.;
        A( row, col( i, j + 1, w_r ) ) = - aR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 2, w_r ) ) =   aR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j + 3, w_r ) ) =   aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s_i )
                   + aR * ( - 5 * Yd * Yd / ( dY * dY ) + 4 * Ydd / ( 2 * dY ) ) * Q( i, j + 1, w_r )
                   - aR * ( - 4 * Yd * Yd / ( dY * dY ) + 1 * Ydd / ( 2 * dY ) ) * Q( i, j + 2, w_r )
                   - aR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, w_r );
        ++row;

        //q_r + (1 / (alpha*Rx^(1/2))) * v_i_{eta eta} = 0
        A( row, col( i, j, q_r ) )     =   1.;
        A( row, col( i, j + 1, v_i ) ) =   6. * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 2, v_i ) ) = - ( 3. / 2. ) * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 3, v_i ) ) =   ( 2. / 9. ) * aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q_r )
                   - aR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j + 1, v_i )
                   + aR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 2, v_i )
                   - aR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, v_i );
        ++row;

        //q_i - (1 / (alpha*Rx^(1/2))) * v_r_{eta eta} = 0
        A( row, col( i, j, q_i ) )     =   1.;
        A( row, col( i, j + 1, v_r ) ) = - 6. * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 2, v_r ) ) =   ( 3. / 2. ) * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j + 3, v_r ) ) = - ( 2. / 9. ) * aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q_i )
                   + aR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j + 1, v_r )
                   - aR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 2, v_r )
                   + aR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j + 3, v_r );
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
          Vector<double> Guess( Q.get_nodes_vars( i, j ) );
          Vector<double> Guess_eta( ( Q.get_nodes_vars( i, j + 1 )
                                    - Q.get_nodes_vars( i, j - 1 ) ) * ( Yd /( 2 * dY )) );
          Vector<double> Guess_hzeta( ( Q.get_nodes_vars( i + 1, j )
                                      - Q.get_nodes_vars( i - 1, j ) )
                                      * ( Xd /( 2 * dX )) );
          Vector<double> Guess_laplace( Q.get_nodes_vars( i, j - 1 ) * laplace_1
                                     +  Q.get_nodes_vars( i - 1, j ) * laplace_3
                                     +  Q.get_nodes_vars( i, j ) * laplace_4
                                     +  Q.get_nodes_vars( i + 1, j ) * laplace_5
                                     +  Q.get_nodes_vars( i, j + 1 ) * laplace_7 );

          Vector<double> Guess_eta_eta( ( Q.get_nodes_vars( i, j + 1 ) * laplace_7
                                        + Q.get_nodes_vars( i, j ) * ( - 2. * Yd * Yd / ( dY * dY ) )
                                        + Q.get_nodes_vars( i, j - 1 ) * laplace_1 ) );
          Vector<double> Guess_hzeta_hzeta( ( Q.get_nodes_vars( i + 1, j ) * laplace_5
                                        + Q.get_nodes_vars( i, j ) * ( - 2. * Xd * Xd / ( ZETA0 * ZETA0 * dX * dX ))
                                        + Q.get_nodes_vars( i - 1, j ) * laplace_3 ) );

          Vector<double> Guess_eta_hzeta( ( Q.get_nodes_vars( i + 1, j + 1 )
                                          + Q.get_nodes_vars( i - 1, j - 1 )
                                          - Q.get_nodes_vars( i + 1, j - 1 )
                                          - Q.get_nodes_vars( i - 1, j + 1 ) ) * ( Xd * Yd / ( 4 * dX * dY) ) );

          double UBdd = BETA * ( Base[ UB ] * Base[ UB ] - 1. ) - Base[ PhiB ] * Base[ UBd ];

          ///////////////////////////////
          // 2D OrrSommerfeld equation //
          ///////////////////////////////

          ///////////////////////
          // v equation (real) //
          ///////////////////////

          // - ( 1 / alpha * Rx^(1/2) ) * Laplacian of v_i
          A( row, col( i, j - 1, v_i ) ) = - aR * laplace_1_os;
          A( row, col( i - 1, j, v_i ) ) = - aR * laplace_3_os;
          A( row, col( i, j, v_i ) )     = - aR * laplace_4_os;
          A( row, col( i + 1, j, v_i ) ) = - aR * laplace_5_os;
          A( row, col( i, j + 1, v_i ) ) = - aR * laplace_7_os;

          // + ( UB + UG - c_r_g ) * v_r
          A( row, col( i, j, v_r ) )    +=   Base[ UB ] + Guess[ U ] - C_GUESS.real();

          // + ((1 - beta) / Rx) * v_r
          A( row, col( i, j, v_r ) )    +=   ( 1.0 - BETA ) / RX;

          // + v_r_g * U
          A( row, col( i, j, U ) )     =   Q( i, j, v_r );

          // - v_r_g * c_r
          A( row, size - 2 )           = - Q( i, j, v_r );

          // + c_i_g * v_i
          A( row, col( i, j, v_i ) )  +=   C_GUESS.imag();

          // + v_i_g * c_i
          A( row, size - 1 )           =   Q( i, j, v_i );

          // - q_r
          A( row, col( i, j, q_r ) )     = - 1.;

          B[ row ] =   aR * ( Q( i, j - 1, v_i ) * laplace_1_os
                           +  Q( i - 1, j, v_i ) * laplace_3_os
                           +  Q( i, j, v_i )     * laplace_4_os
                           +  Q( i + 1, j, v_i ) * laplace_5_os
                           +  Q( i, j + 1, v_i ) * laplace_7_os )
                     - ( Base[ UB ] + Guess[ U ] - C_GUESS.real()
                             + ( ( 1.0 - BETA ) / RX ) ) * Q( i, j, v_r )
                     - C_GUESS.imag() * Q( i, j, v_i )
                     + Q( i, j, q_r );
          ++row;

          ////////////////////////////
          // v equation (imaginary) //
          ////////////////////////////

          // ( 1 / alpha * Rx^(1/2) ) * Laplacian of v_r
          A( row, col( i, j - 1, v_r ) ) = aR * laplace_1_os;
          A( row, col( i - 1, j, v_r ) ) = aR * laplace_3_os;
          A( row, col( i, j, v_r ) )     = aR * laplace_4_os;
          A( row, col( i + 1, j, v_r ) ) = aR * laplace_5_os;
          A( row, col( i, j + 1, v_r ) ) = aR * laplace_7_os;

          // + ( UB + UG - c_r_g ) * v_i
          A( row, col( i, j, v_i ) ) +=   Base[ UB ] + Guess[ U ] - C_GUESS.real();

          // + ((1 - beta) / Rx) * v_i
          A( row, col( i, j, v_i ) ) +=   ( 1.0 - BETA ) / RX;

          // + v_i_g * U
          A( row, col( i, j, U ) )    =   Q( i, j, v_i );

          // - v_i_g * c_r
          A( row, size - 2 )          = - Q( i, j, v_i );

          // - c_i_g * v_r
          A( row, col( i, j, v_r ) ) += - C_GUESS.imag();

          // - v_r_g * c_i
          A( row, size - 1 )          = - Q( i, j, v_r );

          // - q_i
          A( row, col( i, j, q_i ) )  = - 1.;

          B[ row ] =  - aR * ( Q( i, j - 1, v_r ) * laplace_1_os
                            +  Q( i - 1, j, v_r ) * laplace_3_os
                            +  Q( i, j, v_r )     * laplace_4_os
                            +  Q( i + 1, j, v_r ) * laplace_5_os
                            +  Q( i, j + 1, v_r ) * laplace_7_os )
                      - ( Base[ UB ] + Guess[ U ] - C_GUESS.real()
                             + ( ( 1.0 - BETA ) / RX ) ) * Q( i, j, v_i )
                      + C_GUESS.imag() * Q( i, j, v_r )
                      + Q( i, j, q_i );
          ++row;

          ///////////////////////
          // w equation (real) //
          ///////////////////////

          // - ( 1 / alpha * Rx^(1/2) ) * Laplacian of w_i
          A( row, col( i, j - 1, w_i ) ) = - aR * laplace_1_os;
          A( row, col( i - 1, j, w_i ) ) = - aR * laplace_3_os;
          A( row, col( i, j, w_i ) )     = - aR * laplace_4_os;
          A( row, col( i + 1, j, w_i ) ) = - aR * laplace_5_os;
          A( row, col( i, j + 1, w_i ) ) = - aR * laplace_7_os;

          // + ( UB + UG - c_r_g ) * w_r
          A( row, col( i, j, w_r ) )    +=   Base[ UB ] + Guess[ U ] - C_GUESS.real();

          // + ((1 - beta) / Rx) * w_r
          A( row, col( i, j, w_r ) )    +=   ( 1.0 - BETA ) / RX;

          // + w_r_g * U
          A( row, col( i, j, U ) )     =   Q( i, j, w_r );

          // - w_r_g * c_r
          A( row, size - 2 )           = - Q( i, j, w_r );

          // + c_i_g * w_i
          A( row, col( i, j, w_i ) )  +=   C_GUESS.imag();

          // + w_i_g * c_i
          A( row, size - 1 )           =   Q( i, j, w_i );

          // - s_r
          A( row, col( i, j, s_r ) )     = - 1.;

          B[ row ] =   aR * ( Q( i, j - 1, w_i ) * laplace_1_os
                           +  Q( i - 1, j, w_i ) * laplace_3_os
                           +  Q( i, j, w_i )     * laplace_4_os
                           +  Q( i + 1, j, w_i ) * laplace_5_os
                           +  Q( i, j + 1, w_i ) * laplace_7_os )
                     - ( Base[ UB ] + Guess[ U ] - C_GUESS.real()
                             + ( ( 1.0 - BETA ) / RX ) ) * Q( i, j, w_r )
                     - C_GUESS.imag() * Q( i, j, w_i )
                     + Q( i, j, s_r );
          ++row;

          ////////////////////////////
          // w equation (imaginary) //
          ////////////////////////////

          // ( 1 / alpha * Rx^(1/2) ) * Laplacian of v_r
          A( row, col( i, j - 1, w_r ) ) = aR * laplace_1_os;
          A( row, col( i - 1, j, w_r ) ) = aR * laplace_3_os;
          A( row, col( i, j, w_r ) )     = aR * laplace_4_os;
          A( row, col( i + 1, j, w_r ) ) = aR * laplace_5_os;
          A( row, col( i, j + 1, w_r ) ) = aR * laplace_7_os;

          // + ( UB + UG - c_r_g ) * w_i
          A( row, col( i, j, w_i ) ) +=   Base[ UB ] + Guess[ U ] - C_GUESS.real();

          // + ((1 - beta) / Rx) * w_i
          A( row, col( i, j, w_i ) ) +=   ( 1.0 - BETA ) / RX;

          // + w_i_g * U
          A( row, col( i, j, U ) )    =   Q( i, j, w_i );

          // - w_i_g * c_r
          A( row, size - 2 )          = - Q( i, j, w_i );

          // - c_i_g * w_r
          A( row, col( i, j, w_r ) ) += - C_GUESS.imag();

          // - w_r_g * c_i
          A( row, size - 1 )          = - Q( i, j, w_r );

          // - s_i
          A( row, col( i, j, s_i ) )  = - 1.;

          B[ row ] =  - aR * ( Q( i, j - 1, w_r ) * laplace_1_os
                            +  Q( i - 1, j, w_r ) * laplace_3_os
                            +  Q( i, j, w_r )     * laplace_4_os
                            +  Q( i + 1, j, w_r ) * laplace_5_os
                            +  Q( i, j + 1, w_r ) * laplace_7_os )
                      - ( Base[ UB ] + Guess[ U ] - C_GUESS.real()
                             + ( ( 1.0 - BETA ) / RX ) ) * Q( i, j, w_i )
                      + C_GUESS.imag() * Q( i, j, w_r )
                      + Q( i, j, s_i );
          ++row;

          ///////////////////////
          // q equation (real) //
          ///////////////////////

          // q_r_{eta eta}
          A( row, col( i, j - 1, q_r ) ) =   Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ;
          A( row, col( i, j, q_r ) )     = - 2 * Yd * Yd / ( dY * dY );
          A( row, col( i, j + 1, q_r ) ) =   Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ;

          // - alpha^2 * q_r
          A( row, col( i, j, q_r ) )    += - ALPHA * ALPHA;

          // + ( 1 / zeta0 ) * s_r_{hzeta eta}
          A( row, col( i + 1, j + 1, s_r ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, s_r ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, s_r ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, s_r ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // + ( 1 / zeta0 ) * (UB' + UG_{eta}) * w_r_{hzeta}
          A( row, col( i + 1, j, w_r ) )      =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, w_r ) )      = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );

          // + ( 1 / zeta0 ) * w_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      =   Guess_hzeta[ w_r ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, U ) )      = - Guess_hzeta[ w_r ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * U_{hzeta} * w_r_{eta}
          A( row, col( i, j + 1, w_r ) )      = - Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, w_r ) )      =   Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * w_r_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      = - Guess_eta[ w_r ] * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, U ) )      =   Guess_eta[ w_r ] * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * UG_{hzeta eta} * w_r
          A( row, col( i, j, w_r ) )          = - Guess_eta_hzeta[ U ] / ZETA0;

          // - w_r_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  = - Guess[ w_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, U ) )  = - Guess[ w_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, U ) )  =   Guess[ w_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, U ) )  =   Guess[ w_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // - ( UB'' + UG_{eta eta} ) * v_r
          A( row, col( i, j, v_r ) )          = - ( UBdd + Guess_eta_eta[ U ] );

          // - v_r_g * U_{eta eta}
          A( row, col( i, j + 1, U ) )      = - Guess[ v_r ] * laplace_7;
          A( row, col( i, j, U ) )          =   Guess[ v_r ] * 2. * Yd * Yd / ( dY * dY );
          A( row, col( i, j - 1, U ) )      = - Guess[ v_r ] * laplace_1;

          B[ row ] = - ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) * Q( i, j - 1, q_r )
                     - ( - 2 * Yd * Yd / ( dY * dY ) ) * Q( i, j, q_r )
                     - ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) * Q( i, j + 1, q_r )
                     + ALPHA * ALPHA * Q( i, j, q_r )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j + 1, s_r )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j - 1, s_r )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j - 1, s_r )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j + 1, s_r )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i + 1, j, w_r )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i - 1, j, w_r )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j + 1, w_r )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j - 1, w_r )
                     + ( Guess_eta_hzeta[ U ] / ZETA0 ) * Q( i, j, w_r )
                     + ( UBdd + Guess_eta_eta[ U ] ) * Q( i, j, v_r );
          ++row;

          ////////////////////////////
          // q equation (imaginary) //
          ////////////////////////////

          // q_i_{eta eta}
          A( row, col( i, j - 1, q_i ) )    =   Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ;
          A( row, col( i, j, q_i ) )        = - 2 * Yd * Yd / ( dY * dY );
          A( row, col( i, j + 1, q_i ) )    =   Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ;

          // - alpha^2 * q_i
          A( row, col( i, j, q_i ) )         += - ALPHA * ALPHA;

          // + ( 1 / zeta0 ) * s_i_{hzeta eta}
          A( row, col( i + 1, j + 1, s_i ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, s_i ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, s_i ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, s_i ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // + ( 1 / zeta0 ) * (UB' + UG_{eta}) * w_i_{hzeta}
          A( row, col( i + 1, j, w_i ) )      =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, w_i ) )      = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );

          // + ( 1 / zeta0 ) * w_i_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      =   Guess_hzeta[ w_i ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, U ) )      = - Guess_hzeta[ w_i ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * U_{hzeta} * w_i_{eta}
          A( row, col( i, j + 1, w_i ) )      = - Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, w_i ) )      =   Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * w_i_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      = - Guess_eta[ w_i ] * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, U ) )      =   Guess_eta[ w_i ] * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * UG_{hzeta eta} * w_i
          A( row, col( i, j, w_i ) )          = - Guess_eta_hzeta[ U ] / ZETA0;

          // - w_i_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  = - Guess[ w_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, U ) )  = - Guess[ w_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, U ) )  =   Guess[ w_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, U ) )  =   Guess[ w_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // - ( UB'' + UG_{eta eta} ) * v_i
          A( row, col( i, j, v_i ) )          = - ( UBdd + Guess_eta_eta[ U ] );

          // - v_i_g * U_{eta eta}
          A( row, col( i, j + 1, U ) )      = - Guess[ v_i ] * laplace_7;
          A( row, col( i, j, U ) )          =   Guess[ v_i ] * 2. * Yd * Yd / ( dY * dY );
          A( row, col( i, j - 1, U ) )      = - Guess[ v_i ] * laplace_1;

          B[ row ] = - ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) * Q( i, j - 1, q_i )
                     - ( - 2 * Yd * Yd / ( dY * dY ) ) * Q( i, j, q_i )
                     - ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) * Q( i, j + 1, q_i )
                     + ALPHA * ALPHA * Q( i, j, q_i )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j + 1, s_i )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j - 1, s_i )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j - 1, s_i )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j + 1, s_i )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i + 1, j, w_i )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i - 1, j, w_i )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j + 1, w_i )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j - 1, w_i )
                     + ( Guess_eta_hzeta[ U ] / ZETA0 ) * Q( i, j, w_i )
                     + ( UBdd + Guess_eta_eta[ U ] ) * Q( i, j, v_i );
          ++row;

          ///////////////////////
          // s equation (real) //
          ///////////////////////

          // ( 1 / zeta0^2 ) * s_r_{hzeta hzeta}
          A( row, col( i - 1, j, s_r ) )  =   ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) )
                                              / ( ZETA0 * ZETA0 );
          A( row, col( i, j, s_r ) )      = - 2 * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
          A( row, col( i + 1, j, s_r ) )  =   ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) )
                                              / ( ZETA0 * ZETA0 );

          // - alpha^2 * s_r
          A( row, col( i, j, s_r ) )       += - ALPHA * ALPHA;

          // + ( 1 / zeta0 ) * q_r_{hzeta eta}
          A( row, col( i + 1, j + 1, q_r ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, q_r ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, q_r ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, q_r ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // + ( 1 / zeta0 ) * UG_{hzeta} * v_r_{eta}
          A( row, col( i, j + 1, v_r ) )     =   Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, v_r ) )     = - Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );

          // + ( 1 / zeta0 ) * v_r_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      =   Guess_eta[ v_r ] * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, U ) )      = - Guess_eta[ v_r ] * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * (UBd + UG_{eta}) * v_r_{hzeta}
          A( row, col( i + 1, j, v_r ) )     = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, v_r ) )     =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * v_r_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      = - Guess_hzeta[ v_r ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, U ) )      =   Guess_hzeta[ v_r ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * UG_{eta hzeta} * v_r
          A( row, col( i, j, v_r ) )         = - Guess_eta_hzeta[ U ] / ZETA0;

          // - v_r_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  = - Guess[ v_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, U ) )  = - Guess[ v_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, U ) )  =   Guess[ v_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, U ) )  =   Guess[ v_r ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // - ( 1 / zeta0^2 ) * UG_{hzeta hzeta} * w_r
          A( row, col( i, j, w_r ) )         = - Guess_hzeta_hzeta[ U ] / ( ZETA0 * ZETA0 );

          // - w_r_g * U_{hzeta hzeta} / zeta0^2
          A( row, col( i + 1, j, U ) )      = - Guess[ w_r ] * laplace_5;
          A( row, col( i, j, U ) )          =   Guess[ w_r ] * 2. * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
          A( row, col( i - 1, j, U ) )      = - Guess[ w_r ] * laplace_3;

          B[ row ] = - ( ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 1, j, s_r )
                     - ( - 2 * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i, j, s_r )
                     - ( ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i + 1, j, s_r )
                     + ALPHA * ALPHA * Q( i, j, s_r )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j + 1, q_r )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j - 1, q_r )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j - 1, q_r )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j + 1, q_r )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j + 1, v_r )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j - 1, v_r )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i + 1, j, v_r )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i - 1, j, v_r )
                     + ( Guess_eta_hzeta[ U ] / ZETA0 ) * Q( i, j, v_r )
                     + ( Guess_hzeta_hzeta[ U ] / ( ZETA0 * ZETA0 ) ) * Q( i, j, w_r );
                     //TODO replace these with Guess[ ... ] instead of Q( i, j, ...)
          ++row;

          ////////////////////////////
          // s equation (imaginary) //
          ////////////////////////////

          // ( 1 / zeta0^2 ) * s_i_{hzeta hzeta}
          A( row, col( i - 1, j, s_i ) )      =   ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) )
                                              / ( ZETA0 * ZETA0 );
          A( row, col( i, j, s_i ) )          = - 2 * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
          A( row, col( i + 1, j, s_i ) )      =   ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) )
                                              / ( ZETA0 * ZETA0 );

          // - alpha^2 * s_i
          A( row, col( i, j, s_i ) )         += - ALPHA * ALPHA;

          // + ( 1 / zeta0 ) * q_i_{hzeta eta}
          A( row, col( i + 1, j + 1, q_i ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, q_i ) )  =   Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, q_i ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, q_i ) )  = - Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // + ( 1 / zeta0 ) * UG_{hzeta} * v_i_{eta}
          A( row, col( i, j + 1, v_i ) )     =   Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, v_i ) )     = - Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 );

          // + ( 1 / zeta0 ) * v_i_g_{eta} * U_{hzeta}
          A( row, col( i + 1, j, U ) )      =   Guess_eta[ v_i ] * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, U ) )      = - Guess_eta[ v_i ] * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * (UBd + UG_{eta}) * v_i_{hzeta}
          A( row, col( i + 1, j, v_i ) )     = - ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );
          A( row, col( i - 1, j, v_i ) )     =   ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 );

          // - ( 1 / zeta0 ) * v_i_g_hzeta * U_eta
          A( row, col( i, j + 1, U ) )      = - Guess_hzeta[ v_i ] * Yd / ( 2 * dY * ZETA0 );
          A( row, col( i, j - 1, U ) )      =   Guess_hzeta[ v_i ] * Yd / ( 2 * dY * ZETA0 );

          // - ( 1 / zeta0 ) * UG_{eta hzeta} * v_i
          A( row, col( i, j, v_i ) )         = - Guess_eta_hzeta[ U ] / ZETA0;

          // - v_i_g * ( 1 / zeta0 ) * U_{hzeta eta}
          A( row, col( i + 1, j + 1, U ) )  = - Guess[ v_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j - 1, U ) )  = - Guess[ v_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i + 1, j - 1, U ) )  =   Guess[ v_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );
          A( row, col( i - 1, j + 1, U ) )  =   Guess[ v_i ] * Xd * Yd / ( 4 * dX * dY * ZETA0 );

          // - ( 1 / zeta0^2 ) * UG_{hzeta hzeta} * w_i
          A( row, col( i, j, w_i ) )         = - Guess_hzeta_hzeta[ U ] / ( ZETA0 * ZETA0 );

          // - w_i_g * U_{hzeta hzeta} / zeta0^2
          A( row, col( i + 1, j, U ) )      = - Guess[ w_i ] * laplace_5;
          A( row, col( i, j, U ) )          =   Guess[ w_i ] * 2. * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
          A( row, col( i - 1, j, U ) )      = - Guess[ w_i ] * laplace_3;

          B[ row ] = - ( ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 1, j, s_i )
                     - ( - 2 * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i, j, s_i )
                     - ( ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i + 1, j, s_i )
                     + ALPHA * ALPHA * Q( i, j, s_i )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j + 1, q_i )
                     - ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j - 1, q_i )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i + 1, j - 1, q_i )
                     + ( Xd * Yd / ( 4 * dX * dY * ZETA0 ) ) * Q( i - 1, j + 1, q_i )
                     - ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j + 1, v_i )
                     + ( Guess_hzeta[ U ] * Yd / ( 2 * dY * ZETA0 ) ) * Q( i, j - 1, v_i )
                     + ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i + 1, j, v_i )
                     - ( ( Base[ UBd ] + Guess_eta[ U ] ) * Xd / ( 2 * dX * ZETA0 ) ) * Q( i - 1, j, v_i )
                     + ( Guess_eta_hzeta[ U ] / ZETA0 ) * Q( i, j, v_i )
                     + ( Guess_hzeta_hzeta[ U ] / ( ZETA0 * ZETA0 ) ) * Q( i, j, w_i );
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

//Term 1
          // - Rx^-1/6 * Sigma^2 * 4 * ( v_r_g_eta + w_r_g_zeta ) * v_r_zeta
          A( row, col( i + 1, j, v_r ) ) += - Rsigma2 * 4 * ( Guess_eta[ v_r ]
                                            + Guess_hzeta[ w_r ] ) * Xd / ( 2 * dX );
          A( row, col( i - 1, j, v_r ) ) +=   Rsigma2 * 4 * ( Guess_eta[ v_r ]
                                            + Guess_hzeta[ w_r ] ) * Xd / ( 2 * dX );

          // - Rx^-1/6 * Sigma^2 * 4 * v_r_g_zeta * v_r_eta
          A( row, col( i, j + 1, v_r ) ) += - Rsigma2 * 4 * Guess_hzeta[ v_r ]
                                            * Yd / ( 2 * dY );
          A( row, col( i, j - 1, v_r ) ) +=   Rsigma2 * 4 * Guess_hzeta[ v_r ]
                                            * Yd / ( 2 * dY );

          // - Rx^-1/6 * Sigma^2 * 4 * v_r_g_zeta * w_r_zeta
          A( row, col( i + 1, j, w_r ) ) += - Rsigma2 * 4 * Guess_hzeta[ v_r ]
                                            * Xd / ( 2 * dX );
          A( row, col( i - 1, j, w_r ) ) +=   Rsigma2 * 4 * Guess_hzeta[ v_r ]
                                            * Xd / ( 2 * dX );
// Term 2
          // - Rx^-1/6 * Sigma^2 * 4 * ( v_i_g_eta + w_i_g_zeta ) * v_i_zeta
          A( row, col( i + 1, j, v_i ) ) += - Rsigma2 * 4 * ( Guess_eta[ v_i ]
                                            + Guess_hzeta[ w_i ] ) * Xd / ( 2 * dX );
          A( row, col( i - 1, j, v_i ) ) +=   Rsigma2 * 4 * ( Guess_eta[ v_i ]
                                            + Guess_hzeta[ w_i ] ) * Xd / ( 2 * dX );

          // - Rx^-1/6 * Sigma^2 * 4 * v_i_g_zeta * v_i_eta
          A( row, col( i, j + 1, v_i ) ) += - Rsigma2 * 4 * Guess_hzeta[ v_i ]
                                            * Yd / ( 2 * dY );
          A( row, col( i, j - 1, v_i ) ) +=   Rsigma2 * 4 * Guess_hzeta[ v_i ]
                                            * Yd / ( 2 * dY );

          // - Rx^-1/6 * Sigma^2 * 4 * v_i_g_zeta * w_i_zeta
          A( row, col( i + 1, j, w_i ) ) += - Rsigma2 * 4 * Guess_hzeta[ v_i ]
                                            * Xd / ( 2 * dX );
          A( row, col( i - 1, j, w_i ) ) +=   Rsigma2 * 4 * Guess_hzeta[ v_i ]
                                            * Xd / ( 2 * dX );
// Term 3
          // - Rx^-1/6 * Sigma^2 * 2 * ( 2 * v_r_g_eta_zeta + w_r_g_zeta_zeta - w_r_g_eta_eta ) * v_r
          A( row, col( i, j, v_r ) ) += - Rsigma2 * 2 * ( 2 * Guess_eta_hzeta[ v_r ]
                                      + Guess_hzeta_hzeta[ w_r ] - Guess_eta_eta[ w_r ] );

          // - Rx^-1/6 * Sigma^2 * 4 * v_r_g * v_r_eta_zeta
          A( row, col( i + 1, j + 1, v_r ) ) += - Rsigma2 * 4 * ( Guess[ v_r ] )
                                              * Xd * Yd / ( 4 * dX * dY );
          A( row, col( i - 1, j - 1, v_r ) ) += - Rsigma2 * 4 * ( Guess[ v_r ] )
                                              * Xd * Yd / ( 4 * dX * dY );
          A( row, col( i + 1, j - 1, v_r ) ) +=   Rsigma2 * 4 * ( Guess[ v_r ] )
                                              * Xd * Yd / ( 4 * dX * dY );
          A( row, col( i - 1, j + 1, v_r ) ) +=   Rsigma2 * 4 * ( Guess[ v_r ] )
                                              * Xd * Yd / ( 4 * dX * dY );

          // - Rx^-1/6 * Sigma^2 * 2 * v_r_g * w_r_zeta_zeta
          A( row, col( i - 1, j, w_r ) )  += - Rsigma2 * 2 * Guess[ v_r ]
                                           * ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) );
          A( row, col( i, j, w_r ) )      +=   Rsigma2 * 2 * Guess[ v_r ]
                                           * 2 * Xd * Xd / ( dX * dX );
          A( row, col( i + 1, j, w_r ) )  += - Rsigma2 * 2 * Guess[ v_r ]
                                           *  ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) );

         //   Rx^-1/6 * Sigma^2 * 2 * v_r_g * w_r_eta_eta
         A( row, col( i, j - 1, w_r ) )  +=   Rsigma2 * 2 * Guess[ v_r ]
                                          * ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) ;
         A( row, col( i, j, w_r ) )      += - Rsigma2 * 2 * Guess[ v_r ]
                                          * 2 * Yd * Yd / ( dY * dY );
         A( row, col( i, j + 1, w_r ) )  +=   Rsigma2 * 2 * Guess[ v_r ]
                                          * ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) ;

// Term 4
        // - Rx^-1/6 * Sigma^2 * 2 * ( 2 * v_i_g_eta_zeta + w_i_g_zeta_zeta - w_i_g_eta_eta ) * v_i
        A( row, col( i, j, v_i ) ) += - Rsigma2 * 2 * ( 2 * Guess_eta_hzeta[ v_i ]
                                    + Guess_hzeta_hzeta[ w_i ] - Guess_eta_eta[ w_i ] );

        // - Rx^-1/6 * Sigma^2 * 4 * v_i_g * v_i_eta_zeta
        A( row, col( i + 1, j + 1, v_i ) ) += - Rsigma2 * 4 * ( Guess[ v_i ] )
                                            * Xd * Yd / ( 4 * dX * dY );
        A( row, col( i - 1, j - 1, v_i ) ) += - Rsigma2 * 4 * ( Guess[ v_i ] )
                                            * Xd * Yd / ( 4 * dX * dY );
        A( row, col( i + 1, j - 1, v_i ) ) +=   Rsigma2 * 4 * ( Guess[ v_i ] )
                                            * Xd * Yd / ( 4 * dX * dY );
        A( row, col( i - 1, j + 1, v_i ) ) +=   Rsigma2 * 4 * ( Guess[ v_i ] )
                                            * Xd * Yd / ( 4 * dX * dY );

        // - Rx^-1/6 * Sigma^2 * 2 * v_i_g * w_i_zeta_zeta
        A( row, col( i - 1, j, w_i ) )  += - Rsigma2 * 2 * Guess[ v_i ]
                                         * ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) );
        A( row, col( i, j, w_i ) )      +=   Rsigma2 * 2 * Guess[ v_i ]
                                         * 2 * Xd * Xd / ( dX * dX );
        A( row, col( i + 1, j, w_i ) )  += - Rsigma2 * 2 * Guess[ v_i ]
                                         *  ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) );

       //   Rx^-1/6 * Sigma^2 * 2 * v_i_g * w_i_eta_eta
       A( row, col( i, j - 1, w_i ) )  +=   Rsigma2 * 2 * Guess[ v_i ]
                                        * ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) ;
       A( row, col( i, j, w_i ) )      += - Rsigma2 * 2 * Guess[ v_i ]
                                        * 2 * Yd * Yd / ( dY * dY );
       A( row, col( i, j + 1, w_i ) )  +=   Rsigma2 * 2 * Guess[ v_i ]
                                        * ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) ;
// Term 5
       // + Rx^-1/6 * Sigma^2 * 4 * ( w_r_g_zeta + v_r_g_eta ) * w_r_eta
       A( row, col( i, j + 1, w_r ) ) +=   Rsigma2 * 4 * ( Guess_hzeta[ w_r ]
                                         + Guess_eta[ v_r ] ) * Yd / ( 2 * dY );
       A( row, col( i, j - 1, w_r ) ) += - Rsigma2 * 4 * ( Guess_hzeta[ w_r ]
                                         + Guess_eta[ v_r ] ) * Yd / ( 2 * dY );
       // + Rx^-1/6 * Sigma^2 * 4 * w_r_g_eta * w_r_zeta
       A( row, col( i + 1, j, w_r ) ) +=   Rsigma2 * 4 * Guess_eta[ w_r ] * Xd
                                         / ( 2 * dX );
       A( row, col( i - 1, j, w_r ) ) += - Rsigma2 * 4 * Guess_eta[ w_r ] * Xd
                                         / ( 2 * dX );

       // + Rx^-1/6 * Sigma^2 * 4 * w_r_g_eta * v_r_eta
       A( row, col( i, j + 1, v_r ) ) +=   Rsigma2 * 4 * Guess_eta[ w_r ] * Yd
                                         / ( 2 * dY );
       A( row, col( i, j - 1, v_r ) ) += - Rsigma2 * 4 * Guess_eta[ w_r ] * Yd
                                         / ( 2 * dY );
// Term 6
       // + Rx^-1/6 * Sigma^2 * 4 * ( w_i_g_zeta + v_i_g_eta ) * w_i_eta
       A( row, col( i, j + 1, w_i ) ) +=   Rsigma2 * 4 * ( Guess_hzeta[ w_i ]
                                         + Guess_eta[ v_i ] ) * Yd / ( 2 * dY );
       A( row, col( i, j - 1, w_i ) ) += - Rsigma2 * 4 * ( Guess_hzeta[ w_i ]
                                         + Guess_eta[ v_i ] ) * Yd / ( 2 * dY );
       // + Rx^-1/6 * Sigma^2 * 4 * w_i_g_eta * w_i_zeta
       A( row, col( i + 1, j, w_i ) ) +=   Rsigma2 * 4 * Guess_eta[ w_i ] * Xd
                                         / ( 2 * dX );
       A( row, col( i - 1, j, w_i ) ) += - Rsigma2 * 4 * Guess_eta[ w_i ] * Xd
                                         / ( 2 * dX );

       // + Rx^-1/6 * Sigma^2 * 4 * w_r_g_eta * v_r_eta
       A( row, col( i, j + 1, v_i ) ) +=   Rsigma2 * 4 * Guess_eta[ w_i ] * Yd
                                         / ( 2 * dY );
       A( row, col( i, j - 1, v_i ) ) += - Rsigma2 * 4 * Guess_eta[ w_i ] * Yd
                                         / ( 2 * dY );
// Term 7
      // + Rx^-1/6 * Sigma^2 * 2 * ( 2 * w_r_g_eta_zeta + v_r_g_eta_eta - v_r_g_zeta_zeta ) * w_r
      A( row, col( i, j, w_r ) ) +=   Rsigma2 * 2 * ( 2 * Guess_eta_hzeta[ w_r ]
                                  + Guess_eta_eta[ v_r ] - Guess_hzeta_hzeta[ v_r ] );

      // + Rx^-1/6 * Sigma^2 * 4 * w_r_g * w_r_eta_zeta
      A( row, col( i + 1, j + 1, w_r ) ) +=   Rsigma2 * 4 * ( Guess[ w_r ] )
                                          * Xd * Yd / ( 4 * dX * dY );
      A( row, col( i - 1, j - 1, w_r ) ) +=   Rsigma2 * 4 * ( Guess[ w_r ] )
                                          * Xd * Yd / ( 4 * dX * dY );
      A( row, col( i + 1, j - 1, w_r ) ) += - Rsigma2 * 4 * ( Guess[ w_r ] )
                                          * Xd * Yd / ( 4 * dX * dY );
      A( row, col( i - 1, j + 1, w_r ) ) += - Rsigma2 * 4 * ( Guess[ w_r ] )
                                          * Xd * Yd / ( 4 * dX * dY );

      // + Rx^-1/6 * Sigma^2 * 2 * w_r_g * v_r_eta_eta
      A( row, col( i, j - 1, v_r ) )  +=   Rsigma2 * 2 * Guess[ w_r ]
                                       * ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) ;
      A( row, col( i, j, v_r ) )      += - Rsigma2 * 2 * Guess[ w_r ]
                                       * 2 * Yd * Yd / ( dY * dY );
      A( row, col( i, j + 1, v_r ) )  +=   Rsigma2 * 2 * Guess[ w_r ]
                                       * ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) ;

      // - Rx^-1/6 * Sigma^2 * 2 * w_r_g * v_r_zeta_zeta
      A( row, col( i - 1, j, v_r ) )  += - Rsigma2 * 2 * Guess[ w_r ]
                                       * ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) );
      A( row, col( i, j, v_r ) )      +=   Rsigma2 * 2 * Guess[ w_r ]
                                       * 2 * Xd * Xd / ( dX * dX );
      A( row, col( i + 1, j, v_r ) )  += - Rsigma2 * 2 * Guess[ w_r ]
                                       *  ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) );
// Term 8
     // + Rx^-1/6 * Sigma^2 * 2 * ( 2 * w_i_g_eta_zeta + v_i_g_eta_eta - v_i_g_zeta_zeta ) * w_i
     A( row, col( i, j, w_i ) ) +=   Rsigma2 * 2 * ( 2 * Guess_eta_hzeta[ w_i ]
                                 + Guess_eta_eta[ v_i ] - Guess_hzeta_hzeta[ v_i ] );

     // + Rx^-1/6 * Sigma^2 * 4 * w_i_g * w_i_eta_zeta
     A( row, col( i + 1, j + 1, w_i ) ) +=   Rsigma2 * 4 * ( Guess[ w_i ] )
                                         * Xd * Yd / ( 4 * dX * dY );
     A( row, col( i - 1, j - 1, w_i ) ) +=   Rsigma2 * 4 * ( Guess[ w_i ] )
                                         * Xd * Yd / ( 4 * dX * dY );
     A( row, col( i + 1, j - 1, w_i ) ) += - Rsigma2 * 4 * ( Guess[ w_i ] )
                                         * Xd * Yd / ( 4 * dX * dY );
     A( row, col( i - 1, j + 1, w_i ) ) += - Rsigma2 * 4 * ( Guess[ w_i ] )
                                         * Xd * Yd / ( 4 * dX * dY );

     // + Rx^-1/6 * Sigma^2 * 2 * w_i_g * v_i_eta_eta
     A( row, col( i, j - 1, v_i ) )  +=   Rsigma2 * 2 * Guess[ w_i ]
                                      * ( Yd * Yd / ( dY * dY ) - Ydd / ( 2 * dY ) ) ;
     A( row, col( i, j, v_i ) )      += - Rsigma2 * 2 * Guess[ w_i ]
                                      * 2 * Yd * Yd / ( dY * dY );
     A( row, col( i, j + 1, v_i ) )  +=   Rsigma2 * 2 * Guess[ w_i ]
                                      * ( Yd * Yd / ( dY * dY ) + Ydd / ( 2 * dY ) ) ;

     // - Rx^-1/6 * Sigma^2 * 2 * w_i_g * v_i_zeta_zeta
     A( row, col( i - 1, j, v_i ) )  += - Rsigma2 * 2 * Guess[ w_i ]
                                      * ( Xd * Xd / ( dX * dX ) - Xdd / ( 2 * dX ) );
     A( row, col( i, j, v_i ) )      +=   Rsigma2 * 2 * Guess[ w_i ]
                                      * 2 * Xd * Xd / ( dX * dX );
     A( row, col( i + 1, j, v_i ) )  += - Rsigma2 * 2 * Guess[ w_i ]
                                      *  ( Xd * Xd / ( dX * dX ) + Xdd / ( 2 * dX ) );

        double Forcing;
        Forcing  = 4 * Guess_hzeta[ v_r ] * ( Guess_eta[ v_r ] + Guess_hzeta[ w_r ] );
        Forcing += 4 * Guess_hzeta[ v_i ] * ( Guess_eta[ v_i ] + Guess_hzeta[ w_i ] );
        Forcing += 2 * Guess[ v_r ] * ( 2 * Guess_eta_hzeta[ v_r ]
                     + Guess_hzeta_hzeta[ w_r ] - Guess_eta_eta[ w_r ] );
        Forcing += 2 * Guess[ v_i ] * ( 2 * Guess_eta_hzeta[ v_i ]
                     + Guess_hzeta_hzeta[ w_i ] - Guess_eta_eta[ w_i ] );

        Forcing -= 4 * Guess_eta[ w_r ] * ( Guess_hzeta[ w_r ] + Guess_eta[ v_r ] );
        Forcing -= 4 * Guess_eta[ w_i ] * ( Guess_hzeta[ w_i ] + Guess_eta[ v_i ] );
        Forcing -= 2 * Guess[ w_r ] * ( 2 * Guess_eta_hzeta[ w_r ]
                     + Guess_eta_eta[ v_r ] - Guess_hzeta_hzeta[ v_r ] );
        Forcing -= 2 * Guess[ w_i ] * ( 2 * Guess_eta_hzeta[ w_i ]
                     + Guess_eta_eta[ v_i ] - Guess_hzeta_hzeta[ v_i ] );

        //std::cout << " Forcing = " << Forcing << std::endl;


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
                        * Guess[ Theta ] + hzeta * Base[ ThetaB ] * Guess[ U ] )
                        + Rsigma2 * Forcing;
        ++row;

        }

        // eta = eta_inf boundary ( top boundary )
        j = M ;
        eta = ETA_NODES[ j ];
        Yd = SSI.mesh_Yd( eta );

        // v_r = 0
        A( row, col( i, j, v_r ) ) = 1;

        B[ row ] = - Q( i, j, v_r );
        ++row;

        // v_i = 0
        A( row, col( i, j, v_i ) ) = 1;

        B[ row ] = - Q( i, j, v_i );
        ++row;

        // w_r = 0
        A( row, col( i, j, w_r ) ) = 1;

        B[ row ] = - Q( i, j, w_r );
        ++row;

        // w_i = 0
        A( row, col( i, j, w_i ) ) = 1;

        B[ row ] = - Q( i, j, w_i );
        ++row;

        // s_r + (1 / (alpha*Rx^(1/2))) * w_i_{eta eta} = 0
        A( row, col( i, j, s_r ) )     =   1;
        A( row, col( i, j - 1, w_i ) ) =   aR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j - 2, w_i ) ) = - aR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j - 3, w_i ) ) = - aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s_r )
                   - aR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) ) * Q( i, j - 1, w_i )
                   + aR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) ) * Q( i, j - 2, w_i )
                   + aR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, w_i );
        ++row;

        // s_i - (1 / (alpha*Rx^(1/2))) * w_r_{eta eta} = 0
        A( row, col( i, j, s_i ) )     =   1;
        A( row, col( i, j - 1, w_r ) ) = - aR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) );
        A( row, col( i, j - 2, w_r ) ) =   aR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) );
        A( row, col( i, j - 3, w_r ) ) =   aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, s_i )
                   + aR * ( - 5 * Yd * Yd / ( dY * dY ) - 4 * Ydd / ( 2 * dY ) ) * Q( i, j - 1, w_r )
                   - aR * ( - 4 * Yd * Yd / ( dY * dY ) - 1 * Ydd / ( 2 * dY ) ) * Q( i, j - 2, w_r )
                   - aR * ( Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, w_r );
        ++row;

        // q_r + (1 / (alpha*Rx^(1/2))) * v_i_{eta eta} = 0
        A( row, col( i, j, q_r ) )     =   1;
        A( row, col( i, j - 1, v_i ) ) =   6. * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j - 2, v_i ) ) = - ( 3. / 2. ) * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j - 3, v_i ) ) =   ( 2. / 9. ) * aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q_r )
                   - aR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j - 1, v_i )
                   + aR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 2, v_i )
                   - aR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, v_i );
        ++row;

        // q_i - (1 / (alpha*Rx^(1/2))) * v_r_{eta eta} = 0
        A( row, col( i, j, q_i ) )     =   1;
        A( row, col( i, j - 1, v_r ) ) = - 6. * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j - 2, v_r ) ) =   ( 3. / 2. ) * aR * Yd * Yd / ( dY * dY );
        A( row, col( i, j - 3, v_r ) ) = - ( 2. / 9. ) * aR * Yd * Yd / ( dY * dY );

        B[ row ] = - Q( i, j, q_i )
                   + aR * ( 6. * Yd * Yd / ( dY * dY ) ) * Q( i, j - 1, v_r )
                   - aR * ( ( 3. / 2. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 2, v_r )
                   + aR * ( ( 2. / 9. ) * Yd * Yd / ( dY * dY ) ) * Q( i, j - 3, v_r );
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

        // w_r = 0
        A( row, col( i, j, w_r ) ) = 1;

        B[ row ] = - Q( i, j, w_r );
        ++row;

        // w_i = 0
        A( row, col( i, j, w_i ) ) = 1;

        B[ row ] = - Q( i, j, w_i );
        ++row;

        // v_r = 0
        A( row, col( i, j, v_r ) ) = 1;

        B[ row ] = - Q( i, j, v_r );
        ++row;

        // v_i = 0
        A( row, col( i, j, v_i ) ) = 1;

        B[ row ] = - Q( i, j, v_i );
        ++row;

        // q_r + (1/zeta0^2) * (1 / (alpha*Rx^(1/2))) * v_i_{hzeta hzeta} = 0
        A( row, col( i, j, q_r ) )        =   1;
        A( row, col( i - 1, j, v_i ) )    =   aR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) )
                                      / ( ZETA0 * ZETA0 );
        A( row, col( i - 2, j, v_i ) )    = - aR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) )
                                      / ( ZETA0 * ZETA0 );
        A( row, col( i - 3, j, v_i ) )    = - aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );

        B[ row ] = - Q( i, j, q_r )
                   - ( aR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 1, j, v_i )
                   + ( aR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 2, j, v_i )
                   + ( aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 3, j, v_i );
        ++row;

        // q_i - (1/zeta0^2) * (1 / (alpha*Rx^(1/2))) * v_r_{hzeta hzeta} = 0
        A( row, col( i, j, q_i ) )     =   1;
        A( row, col( i - 1, j, v_r ) ) = - aR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) )
                                          / ( ZETA0 * ZETA0 );
        A( row, col( i - 2, j, v_r ) ) =   aR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) )
                                          / ( ZETA0 * ZETA0 );
        A( row, col( i - 3, j, v_r ) ) =   aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );

        B[ row ] = - Q( i, j, q_i )
                   + ( aR * ( - 5 * Xd * Xd / ( dX * dX ) - 4 * Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 1, j, v_r )
                   - ( aR * ( - 4 * Xd * Xd / ( dX * dX ) - 1 * Xdd / ( 2 * dX ) ) / ( ZETA0 * ZETA0 ) ) * Q( i - 2, j, v_r )
                   - ( aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 3, j, v_r );
        ++row;

        // s_r + (1/zeta0^2) * (1 / (alpha*Rx^(1/2))) * w_i_{hzeta hzeta} = 0
        A( row, col( i, j, s_r ) )     =   1;
        A( row, col( i - 1, j, w_i ) ) =   6. * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
        A( row, col( i - 2, j, w_i ) ) = - ( 3. / 2. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
        A( row, col( i - 3, j, w_i ) ) =   ( 2. / 9. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );

        B[ row ] = - Q( i, j, s_r )
                   - ( 6. * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 1, j, w_i )
                   + ( ( 3. / 2. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 2, j, w_i )
                   - ( ( 2. / 9. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 3, j, w_i );
        ++row;

        // s_i - (1/zeta0^2) * (1 / (alpha*Rx^(1/2))) * w_r_{hzeta hzeta} = 0
        A( row, col( i, j, s_i ) )     =   1;
        A( row, col( i - 1, j, w_r ) ) = - 6. * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
        A( row, col( i - 2, j, w_r ) ) =   ( 3. / 2. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );
        A( row, col( i - 3, j, w_r ) ) = - ( 2. / 9. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 );

        B[ row ] = - Q( i, j, s_i )
                   + ( 6. * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 1, j, w_r )
                   - ( ( 3. / 2. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 2, j, w_r )
                   + ( ( 2. / 9. ) * aR * Xd * Xd / ( dX * dX * ZETA0 * ZETA0 ) ) * Q( i - 3, j, w_r );
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

      // q_r_eta(0,0) = 1e-5 (extra condition)
      double eta( ETA_NODES[ 0 ] );
      double Yd( SSI.mesh_Yd( eta ) );
      double Ydd( SSI.mesh_Ydd( eta ) );
      A( row, col( 0, 0, q_r ) ) = - 3 * Yd / ( 2 * dY );
      A( row, col( 0, 1, q_r ) ) =   4 * Yd / ( 2 * dY );
      A( row, col( 0, 2, q_r ) ) = - 1 * Yd / ( 2 * dY );
      B[ row ] = 1e-5 + ( 3 * Yd / ( 2 * dY ) ) * Q( 0, 0, q_r )
                      - ( 4 * Yd / ( 2 * dY ) ) * Q( 0, 1, q_r )
                      + ( 1 * Yd / ( 2 * dY ) ) * Q( 0, 2, q_r );
      ++row;

      // q_i_eta(0,0) = 0 (extra condition)
      A( row, col( 0, 0, q_i ) ) = - 3 * Yd / ( 2 * dY );
      A( row, col( 0, 1, q_i ) ) =   4 * Yd / ( 2 * dY );
      A( row, col( 0, 2, q_i ) ) = - 1 * Yd / ( 2 * dY );
      B[ row ] =       ( 3 * Yd / ( 2 * dY ) ) * Q( 0, 0, q_i )
                     - ( 4 * Yd / ( 2 * dY ) ) * Q( 0, 1, q_i )
                     + ( 1 * Yd / ( 2 * dY ) ) * Q( 0, 2, q_i );
      ++row;

      // Extra condition (integral on centreline) ???
      /*double int_r, int_i;
      double eta, Yd;

      for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
      {
        eta = ETA_NODES[ j ];
        Yd = SSI.mesh_Yd( eta );
        int_r += 0.5 * dY * ( Q( 0, j, v_r ) + Q( 0, j + 1, v_r ) ) / Yd;
        int_i += 0.5 * dY * ( Q( 0, j, v_i ) + Q( 0, j + 1, v_i ) ) / Yd;
      }

      // int v_r = 1 on centreline
      for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
      {
        eta = ETA_NODES[ j ];
        Yd = SSI.mesh_Yd( eta );
        A( row, col( 0, j, v_r ) )      += 0.5 * dY / Yd;
        A( row, col( 0, j + 1, v_r ) )  += 0.5 * dY / Yd;
      }

      B[ row ] = 1.0 - int_r;
      ++row;

      // int v_i = 0 on centreline
      for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
      {
        eta = ETA_NODES[ j ];
        Yd = SSI.mesh_Yd( eta );
        A( row, col( 0, j, v_i ) )      += 0.5 * dY / Yd;
        A( row, col( 0, j + 1, v_i ) )  += 0.5 * dY / Yd;
      }

      B[ row ] = - int_i;
      ++row;*/

      max_residual = B.norm_inf();
      std::cout << "***                                              Maximum residual = "
                << max_residual << std::endl;

      // Solve the sparse system
      Vector<double> x;
      if( SPEED_UP )
      {
        // Convert things to Eigen to see if we can get some speed benefits
        // Only decompose the matrix on the first iteration as it can be reused after
        // This means more iterations but they should be quicker
        if ( iteration == 0 )
        {
          A_Eigen = A.convert_to_Eigen();
          solver.compute( A_Eigen );
        }

        Eigen::Matrix<double, -1, 1> B_Eigen( 4 * ( M + 1 ) * ( N + 1 ) );
        B_Eigen = B.convert_to_Eigen_Matrix();

        Eigen::Matrix<double, -1, 1> x_Eigen( 4 * ( M + 1 ) * ( N + 1 ) );
        x_Eigen = solver.solve( B_Eigen );
        x = x.convert_to_Vector( x_Eigen );
      } else {
        x = A.solve( B );
      }
      B = x;

      timer.print();
      timer.stop();

      // Update the known values using the correction which we just found
      for ( std::size_t i = 0; i < N + 1; ++i )
      {
        for ( std::size_t j = 0; j < M + 1; ++j )
        {
          Q( i, j, v_r )   += B[ col( i, j, v_r ) ];
          Q( i, j, v_i )   += B[ col( i, j, v_i ) ];
          Q( i, j, w_r )   += B[ col( i, j, w_r ) ];
          Q( i, j, w_i )   += B[ col( i, j, w_i ) ];
          Q( i, j, q_r )   += B[ col( i, j, q_r ) ];
          Q( i, j, q_i )   += B[ col( i, j, q_i ) ];
          Q( i, j, s_r )   += B[ col( i, j, s_r ) ];
          Q( i, j, s_i )   += B[ col( i, j, s_i ) ];
          Q( i, j, Phi )   += B[ col( i, j, Phi ) ];
          Q( i, j, Psi )   += B[ col( i, j, Psi ) ];
          Q( i, j, U )     += B[ col( i, j, U ) ];
          Q( i, j, Theta ) += B[ col( i, j, Theta ) ];
        }
      }

      double c_g_r_temp( C_GUESS.real() );
      double c_g_i_temp( C_GUESS.imag() );
      C_GUESS.real( c_g_r_temp + B[ size - 2 ] );
      C_GUESS.imag( c_g_i_temp + B[ size - 1 ] );

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
        // First 8 values are the wave
        Q_OUTPUT( i, j, 0 )  = Q( i, j, v_r );
        Q_OUTPUT( i, j, 1 )  = Q( i, j, v_i );
        Q_OUTPUT( i, j, 2 )  = Q( i, j, w_r );
        Q_OUTPUT( i, j, 3 )  = Q( i, j, w_i );
        Q_OUTPUT( i, j, 4 )  = Q( i, j, q_r );
        Q_OUTPUT( i, j, 5 )  = Q( i, j, q_i );
        Q_OUTPUT( i, j, 6 )  = Q( i, j, s_r );
        Q_OUTPUT( i, j, 7 )  = Q( i, j, s_i );
        // next 4 values output are the streak without the underlying base flow
        Q_OUTPUT( i, j, 8 )  = Q( i, j, Phi );
        Q_OUTPUT( i, j, 9 )  = Q( i, j, Psi );
        Q_OUTPUT( i, j, 10 ) = Q( i, j, U );
        Q_OUTPUT( i, j, 11 ) = Q( i, j, Theta );
        // final 4 values are the "full" solution, but still with the zeta0 scaling
        Q_OUTPUT( i, j, 12 ) =   Q( i, j, Phi )
                            + BASE_SOLUTION.get_interpolated_vars( eta )[PhiB];
        Q_OUTPUT( i, j, 13 ) = Q( i, j, Psi )
                            + hzeta * BASE_SOLUTION.get_interpolated_vars( eta )[PsiB];
        Q_OUTPUT( i, j, 14 ) =   Q( i, j, U )
                            + BASE_SOLUTION.get_interpolated_vars( eta )[UB];
        Q_OUTPUT( i, j, 15 ) = Q( i, j, Theta )
                            + hzeta * BASE_SOLUTION.get_interpolated_vars( eta )[ThetaB];
      }
    }

    SOLVED = true;

  } // End of solve_local method

} // End of namespace TSL

#endif
