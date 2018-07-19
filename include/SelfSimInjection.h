/* SelfSimInjection - Here we define the SelfSimInjection class used for solving
      the self similar short scale injection problem.
*/

#ifndef SELFSIMINJECTION_H
#define SELFSIMINJECTION_H

#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <sstream>
#include <string>

#include "Vector.h"
#include "Residual.h"
#include "Newton.h"
#include "Equation.h"
#include "Arclength.h"
#include "Error.h"

// Enumerations
enum{ f, fd, fdd, g, gd, gdd };                               // Base ODE
enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };                 // Base ODE
enum class SelfSim{ Phi, Psi, U, Theta };                     // PDE

namespace TSL
{
  class SelfSimInjection
  {
  protected:
      double HZETA_RIGHT;         // Size of the domain in the zeta_hat direction
      double ETA_TOP;             // Size of the domain in the eta direction
      std::size_t N;              // Number of intervals in the zeta_hat direction
      std::size_t M;              // Number of intervals in the eta direction
      std::size_t MB;             // Number of eta intervals in the base flow ODE
      double BETA;                // Hartree parameter
      double KB;                  // Base flow transpiration ( +ve = blowing )
      double ZETA0;               // Injection width
      double K;                   // Injection parameter ( +ve = blowing )
      std::string MESH;           // Mesh definition (UNIFORM or NONUNIFORM)
      std::string BASE_FLOW;      // 2D Falkner-Skan or 3D alternative base flow
      bool SPEED_UP;              // Reuse the factorised matrix for speed gains
      std::string OUTPUT_PATH;    // Output path string
      bool SOLVED;                // True if equations have been solved
      bool OUTPUT;                // Output solution if True

      // Wave forcing terms
      bool FORCING;                                       // True if forcing on
      double SIGMA;                                       // Wave amplitude
      double RX;                                          // Local Reynolds number
      double ALPHA;                                       // Wavenumber
      TwoD_node_mesh< std::complex<double> > v_WAVE;      // Transverse wave velocity
      TwoD_node_mesh< std::complex<double> > w_WAVE;      // Spanwise wave velocity
      //TwoD_node_mesh< std::complex<double> > v_WAVE_CONJ; // Transverse wave velocity conjugate
      //TwoD_node_mesh< std::complex<double> > w_WAVE_CONJ; // Spanwise wave velocity conjugate

      // Non-uniform mesh parameters
      double a1, a2, b1, b2;

      // Transpiration function parameters
      double GAMMA;         // Steepness of the injection
      std::size_t N_TRANSP;    // Number of blowing and sucking regions (1 = standard blowing)

      // Nodal positions in the remapped domain (spanned by X,Y)
      Vector<double> X_NODES, Y_NODES;
      // Nodal positions in the original zeta-eta domain
      Vector<double> ETA_NODES, HZETA_NODES;
      Vector<double> BASE_ETA_NODES;
      OneD_node_mesh<double> BASE_SOLUTION;     // Base flow ODE solution
      TwoD_node_mesh<double> Q;                 // Current guess mesh
      TwoD_node_mesh<double> Q_output;          // Output mesh
      OneD_node_mesh<double> WALL_SHEAR;        // Wall shear
      double A_PARAM;                           // Mass flux parameter
      double U_ETA;                             // Shear stress at origin
      double ETA_HALF;                          // Value of eta on hzeta=0 at which U=1/2
      double INT_U_SQUARED;                     // Integral of U^2 over cross section

      /* ----- Methods (protected)----- */
      void setup(){
        if ( OUTPUT )
        {
          make_output_directory();
        }
        mesh_setup();
        solve_base_flow();
      }

      double calculate_forcing_1( std::size_t i, std::size_t j )
      {
        const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
        const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );

        // hzeta location
        double hzeta( HZETA_NODES[ i ] );
        double Xd( mesh_Xd( hzeta ) );
        double Xdd( mesh_Xdd( hzeta ) );

        // eta location
        double eta( ETA_NODES[ j ] );
        double Yd( mesh_Yd( eta ) );
        double Ydd( mesh_Ydd( eta ) );

        double F_1;

        if ( FORCING ){
          double v_real, v_imag, w_real, w_imag;
          double u_real_eta, v_imag_eta_eta, w_imag_zeta_eta;
          double u_imag_eta, v_real_eta_eta, w_real_zeta_eta;
          double u_real_zeta, v_imag_eta_zeta, w_imag_zeta_zeta;
          double u_imag_zeta, v_real_eta_zeta, w_real_zeta_zeta;

          v_real = v_WAVE( i, j, 0 ).real();
          v_imag = v_WAVE( i, j, 0 ).imag();
          w_real = w_WAVE( i, j, 0 ).real();
          w_imag = w_WAVE( i, j, 0 ).imag();

          v_imag_eta_eta = Ydd * ( v_WAVE( i, j + 1, 0 ).imag() - v_WAVE( i, j - 1, 0 ).imag() ) / ( 2 * dY );
          v_imag_eta_eta += Yd * Yd * ( v_WAVE( i, j + 1, 0 ).imag() - 2.0 * v_WAVE( i, j, 0 ).imag() + v_WAVE( i, j - 1, 0 ).imag() ) / ( dY * dY );
          w_imag_zeta_eta = Xd * Yd * ( w_WAVE( i + 1, j + 1, 0 ).imag() + w_WAVE( i - 1, j - 1, 0 ).imag()
                               - w_WAVE( i + 1, j - 1, 0 ).imag() - w_WAVE( i - 1, j + 1, 0 ).imag()  ) / ( 4. * dX * dY );
          u_real_eta = - ( v_imag_eta_eta + w_imag_zeta_eta ) / ALPHA;

          v_real_eta_eta = Ydd * ( v_WAVE( i, j + 1, 0 ).real() - v_WAVE( i, j - 1, 0 ).real() ) / ( 2 * dY );
          v_real_eta_eta += Yd * Yd * ( v_WAVE( i, j + 1, 0 ).real() - 2.0 * v_WAVE( i, j, 0 ).real() + v_WAVE( i, j - 1, 0 ).real() ) / ( dY * dY );
          w_real_zeta_eta = Xd * Yd * ( w_WAVE( i + 1, j + 1, 0 ).real() + w_WAVE( i - 1, j - 1, 0 ).real()
                               - w_WAVE( i + 1, j - 1, 0 ).real() - w_WAVE( i - 1, j + 1, 0 ).real()  ) / ( 4. * dX * dY );
          u_imag_eta = ( v_real_eta_eta + w_real_zeta_eta ) / ALPHA;

          v_imag_eta_zeta = Xd * Yd * ( v_WAVE( i + 1, j + 1, 0 ).imag() + v_WAVE( i - 1, j - 1, 0 ).imag()
                               - v_WAVE( i + 1, j - 1, 0 ).imag() - v_WAVE( i - 1, j + 1, 0 ).imag()  ) / ( 4. * dX * dY );
          w_imag_zeta_zeta = Xdd * ( w_WAVE( i + 1, j, 0 ).imag() - w_WAVE( i - 1, j, 0 ).imag() ) / ( 2 * dX );
          w_imag_zeta_zeta += Xd * Xd * ( w_WAVE( i + 1, j, 0 ).imag() - 2.0 * w_WAVE( i, j, 0 ).imag() + w_WAVE( i - 1, j, 0 ).imag() ) / ( dX * dX );
          u_real_zeta = - ( v_imag_eta_zeta + w_imag_zeta_zeta ) / ALPHA;

          v_real_eta_zeta = Xd * Yd * ( v_WAVE( i + 1, j + 1, 0 ).real() + v_WAVE( i - 1, j - 1, 0 ).real()
                               - v_WAVE( i + 1, j - 1, 0 ).real() - v_WAVE( i - 1, j + 1, 0 ).real()  ) / ( 4. * dX * dY );
          w_real_zeta_zeta = Xdd * ( w_WAVE( i + 1, j, 0 ).real() - w_WAVE( i - 1, j, 0 ).real() ) / ( 2 * dX );
          w_real_zeta_zeta += Xd * Xd * ( w_WAVE( i + 1, j, 0 ).real() - 2.0 * w_WAVE( i, j, 0 ).real() + w_WAVE( i - 1, j, 0 ).real() ) / ( dX * dX );
          u_imag_zeta = ( v_real_eta_zeta + w_real_zeta_zeta ) / ALPHA;

          F_1 = 2 * ( v_real * u_real_eta + v_imag * u_imag_eta + w_real * u_real_zeta + w_imag * u_imag_zeta  );
        } else {
          F_1 = 0;
        }
        return F_1;
      }

      double calculate_forcing_2( std::size_t i, std::size_t j ){
        const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
        const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );

        // hzeta location
        double hzeta( HZETA_NODES[ i ] );
        double Xd( mesh_Xd( hzeta ) );
        double Xdd( mesh_Xdd( hzeta ) );

        // eta location
        double eta( ETA_NODES[ j ] );
        double Yd( mesh_Yd( eta ) );
        double Ydd( mesh_Ydd( eta ) );

        double F_2;

        if( FORCING ){
          double v_real, v_imag, w_real, w_imag;
          double v_real_eta, v_imag_eta, w_real_eta, w_imag_eta;
          double v_real_zeta, v_imag_zeta, w_real_zeta, w_imag_zeta;
          double v_real_eta_eta, v_imag_eta_eta, w_real_eta_eta, w_imag_eta_eta;
          double v_real_zeta_zeta, v_imag_zeta_zeta, w_real_zeta_zeta, w_imag_zeta_zeta;
          double v_real_eta_zeta, v_imag_eta_zeta, w_real_eta_zeta, w_imag_eta_zeta;
          double u_real_eta, u_imag_eta, u_real_zeta, u_imag_zeta;

          v_real = v_WAVE( i, j, 0 ).real();
          v_imag = v_WAVE( i, j, 0 ).imag();
          w_real = w_WAVE( i, j, 0 ).real();
          w_imag = w_WAVE( i, j, 0 ).imag();

          v_real_eta = Yd * ( v_WAVE( i, j + 1, 0 ).real() - v_WAVE( i, j - 1, 0 ).real() ) / dY;
          v_imag_eta = Yd * ( v_WAVE( i, j + 1, 0 ).imag() - v_WAVE( i, j - 1, 0 ).imag() ) / dY;
          w_real_eta = Yd * ( w_WAVE( i, j + 1, 0 ).real() - w_WAVE( i, j - 1, 0 ).real() ) / dY;
          w_imag_eta = Yd * ( w_WAVE( i, j + 1, 0 ).imag() - w_WAVE( i, j - 1, 0 ).imag() ) / dY;

          v_real_zeta = Xd * ( v_WAVE( i + 1, j, 0 ).real() - v_WAVE( i - 1, j, 0 ).real() ) / dX;
          v_imag_zeta = Xd * ( v_WAVE( i + 1, j, 0 ).imag() - v_WAVE( i - 1, j, 0 ).imag() ) / dX;
          w_real_zeta = Xd * ( w_WAVE( i + 1, j, 0 ).real() - w_WAVE( i - 1, j, 0 ).real() ) / dX;
          w_imag_zeta = Xd * ( w_WAVE( i + 1, j, 0 ).imag() - w_WAVE( i - 1, j, 0 ).imag() ) / dX;

          v_real_eta_eta  = Ydd * ( v_WAVE( i, j + 1, 0 ).real() - v_WAVE( i, j - 1, 0 ).real() ) / ( 2 * dY );
          v_real_eta_eta += Yd * Yd * ( v_WAVE( i, j + 1, 0 ).real() - 2.0 * v_WAVE( i, j, 0 ).real() + v_WAVE( i, j - 1, 0 ).real() ) / ( dY * dY );
          v_imag_eta_eta  = Ydd * ( v_WAVE( i, j + 1, 0 ).imag() - v_WAVE( i, j - 1, 0 ).imag() ) / ( 2 * dY );
          v_imag_eta_eta += Yd * Yd * ( v_WAVE( i, j + 1, 0 ).imag() - 2.0 * v_WAVE( i, j, 0 ).imag() + v_WAVE( i, j - 1, 0 ).imag() ) / ( dY * dY );
          w_real_eta_eta  = Ydd * ( w_WAVE( i, j + 1, 0 ).real() - w_WAVE( i, j - 1, 0 ).real() ) / ( 2 * dY );
          w_real_eta_eta += Yd * Yd * ( w_WAVE( i, j + 1, 0 ).real() - 2.0 * w_WAVE( i, j, 0 ).real() + w_WAVE( i, j - 1, 0 ).real() ) / ( dY * dY );
          w_imag_eta_eta  = Ydd * ( w_WAVE( i, j + 1, 0 ).imag() - w_WAVE( i, j - 1, 0 ).imag() ) / ( 2 * dY );
          w_imag_eta_eta += Yd * Yd * ( w_WAVE( i, j + 1, 0 ).imag() - 2.0 * w_WAVE( i, j, 0 ).imag() + w_WAVE( i, j - 1, 0 ).imag() ) / ( dY * dY );

          v_real_zeta_zeta  = Xdd * ( v_WAVE( i + 1, j, 0 ).real() - v_WAVE( i - 1, j, 0 ).real() ) / ( 2 * dX );
          v_real_zeta_zeta += Xd * Xd * ( v_WAVE( i + 1, j, 0 ).real() - 2.0 * v_WAVE( i, j, 0 ).real() + v_WAVE( i - 1, j, 0 ).real() ) / ( dX * dX );
          v_imag_zeta_zeta  = Xdd * ( v_WAVE( i + 1, j, 0 ).imag() - v_WAVE( i - 1, j, 0 ).imag() ) / ( 2 * dX );
          v_imag_zeta_zeta += Xd * Xd * ( v_WAVE( i + 1, j, 0 ).imag() - 2.0 * v_WAVE( i, j, 0 ).imag() + v_WAVE( i - 1, j, 0 ).imag() ) / ( dX * dX );
          w_real_zeta_zeta  = Xdd * ( w_WAVE( i + 1, j, 0 ).real() - w_WAVE( i - 1, j, 0 ).real() ) / ( 2 * dX );
          w_real_zeta_zeta += Xd * Xd * ( w_WAVE( i + 1, j, 0 ).real() - 2.0 * w_WAVE( i, j, 0 ).real() + w_WAVE( i - 1, j, 0 ).real() ) / ( dX * dX );
          w_imag_zeta_zeta  = Xdd * ( w_WAVE( i + 1, j, 0 ).imag() - w_WAVE( i - 1, j, 0 ).imag() ) / ( 2 * dX );
          w_imag_zeta_zeta += Xd * Xd * ( w_WAVE( i + 1, j, 0 ).imag() - 2.0 * w_WAVE( i, j, 0 ).imag() + w_WAVE( i - 1, j, 0 ).imag() ) / ( dX * dX );

          v_real_eta_zeta = Xd * Yd * ( v_WAVE( i + 1, j + 1, 0 ).real() + v_WAVE( i - 1, j - 1, 0 ).real()
                               - v_WAVE( i + 1, j - 1, 0 ).real() - v_WAVE( i - 1, j + 1, 0 ).real()  ) / ( 4. * dX * dY );
          v_imag_eta_zeta = Xd * Yd * ( v_WAVE( i + 1, j + 1, 0 ).imag() + v_WAVE( i - 1, j - 1, 0 ).imag()
                               - v_WAVE( i + 1, j - 1, 0 ).imag() - v_WAVE( i - 1, j + 1, 0 ).imag()  ) / ( 4. * dX * dY );
          w_real_eta_zeta = Xd * Yd * ( w_WAVE( i + 1, j + 1, 0 ).real() + w_WAVE( i - 1, j - 1, 0 ).real()
                               - w_WAVE( i + 1, j - 1, 0 ).real() - w_WAVE( i - 1, j + 1, 0 ).real()  ) / ( 4. * dX * dY );
          w_imag_eta_zeta = Xd * Yd * ( w_WAVE( i + 1, j + 1, 0 ).imag() + w_WAVE( i - 1, j - 1, 0 ).imag()
                               - w_WAVE( i + 1, j - 1, 0 ).imag() - w_WAVE( i - 1, j + 1, 0 ).imag()  ) / ( 4. * dX * dY );

          u_real_eta = - ( v_imag_eta_eta + w_imag_eta_zeta ) / ALPHA;
          u_imag_eta = ( v_real_eta_eta + w_real_eta_zeta ) / ALPHA;
          u_real_zeta = - ( v_imag_eta_zeta + w_imag_zeta_zeta ) / ALPHA;
          u_imag_zeta = ( v_real_eta_zeta + w_real_zeta_zeta ) / ALPHA;

          F_2  = u_imag_zeta * v_real - u_real_zeta * v_imag;
          F_2 += u_real_eta * w_imag - u_imag_eta * w_real;
          F_2 *= ALPHA;
          F_2 += 2 * ( v_real_eta + w_real_zeta ) * ( v_real_zeta - w_real_eta );
          F_2 -= 2 * ( v_imag_eta + w_imag_zeta ) * ( w_imag_eta - v_imag_zeta );
          F_2 += v_real * ( v_real_eta_zeta - w_real_eta_eta );
          F_2 += v_imag * ( v_imag_eta_zeta - w_imag_eta_eta );
          F_2 += w_real * ( v_real_zeta_zeta - w_real_eta_zeta );
          F_2 += w_imag * ( v_imag_zeta_zeta - w_imag_eta_zeta );
          F_2 *= 2;

        } else {
          F_2 = 0;
        }
        return F_2;
      }

      /* ----- Make the output directory ----- */

      void make_output_directory()
      {
        int status = mkdir( OUTPUT_PATH.c_str(), S_IRWXU );
        if ( status == 0 ) {
        std::cout << "  * Output directory " + OUTPUT_PATH +
                " has been made successfully." << std::endl;
        }
      }

      /* ----- Calculate the base flow solution ----- */

      class equation_2D : public Equation<double>
      {
        public:
          double beta;                            // Hartree parameter
          // Falkner-Skan equation is 3rd order
          equation_2D() : Equation<double> ( 3 ) {}
          // Define the equation
          void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
          {
            F[ f ]   = u[ fd ];
            F[ fd ]  = u[ fdd ];
            F[ fdd ] = - u[ f ] * u[ fdd ] - beta * ( 1.0 - u[ fd ] * u[ fd ] );
          }
      }; // End of equation_2D class

      class plate_BC_2D : public Residual<double>
      {
        public:
          double KB;                              // Transpiration parameter

          plate_BC_2D() : Residual<double> ( 2, 3 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + KB;
            B[ 1 ] = z[ fd ];
          }
      }; // End Falkner-Skan plate_BC_2D class

      class far_BC_2D : public Residual<double>
      {
        public:
          far_BC_2D() : Residual<double> ( 1, 3 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ fd ] - 1.0;
          }
      }; // End Falkner-Skan far_BC_2D class

      class equation_3D : public Equation<double>
      {
        public:
          double beta;                     // Hartree parameter
          // The 3D alternative equation is 6th order
          equation_3D() : Equation<double> ( 6 ) {}
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
      }; // End 3D alternative equation_3D class

      class plate_BC_3D : public Residual<double>
      {
        public:
          double KB;                        // Transpiration parameter

          plate_BC_3D() : Residual<double> ( 4, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ f ] + KB;
            B[ 1 ] = z[ fd ];
            B[ 2 ] = z[ g ];
            B[ 3 ] = z[ gd ];
          }
      }; // End 3D alternative plate_BC_3D class

      class far_BC_3D : public Residual<double>
      {
        public:
          far_BC_3D() : Residual<double> ( 2, 6 ) {}

          void residual_fn( const Vector<double> &z, Vector<double> &B ) const
          {
            B[ 0 ] = z[ fd ] - 1.0;
            B[ 1 ] = z[ gd ];
          }
      }; // End 3D alternative far_BC_3D class

      void solve_base_flow()
      {
        std::cout << "*** Solving the base flow ODE ***" << std::endl;
        Vector<double> eta_nodes( BASE_ETA_NODES );
        //eta_nodes.linspace( 0.0, ETA_TOP, MB + 1 );

        if ( BASE_FLOW=="2D" )
        {
          equation_2D equation;
          plate_BC_2D plate_BC;
          far_BC_2D far_BC;
          equation.beta = 0.1;
          plate_BC.KB = 0.0;
          //ODE_BVP<double> base( &equation, ETA_NODES, &plate_BC, &far_BC );
          ODE_BVP<double> base( &equation, eta_nodes, &plate_BC, &far_BC );
          for (std::size_t j=0; j < MB + 1; ++j )
        	{
            double eta = eta_nodes[ j ];
        		base.solution()( j, f )  	= eta + exp( -eta );
            base.solution()( j, fd ) 	= 1.0 - exp( -eta );
        		base.solution()( j, fdd ) = exp( -eta );
        	}

          // Solve the system with KB = 0 then arc-length continue until KB = Param::KB
          double arc_step( 0.01 );
          double max_arc_step( 0.1 );
          base.init_arc( &plate_BC.KB, arc_step, max_arc_step );
          do
          {
            arc_step = base.arclength_solve( arc_step );
          }while( plate_BC.KB < KB );
          plate_BC.KB = KB;
          base.solve_bvp();                               // Solve once more with KB = KB

          // Solve the system with beta = 0.1 then arc-length continue until beta = BETA
          arc_step = -0.01;
          if ( BETA >= 0.1 ) { arc_step = 0.01; }

          base.init_arc( &equation.beta, arc_step, max_arc_step );
          do
          {
            arc_step = base.arclength_solve( arc_step );
          }while( equation.beta < BETA );
          equation.beta = BETA;
          base.solve_bvp();

          OneD_node_mesh<double> Base_soln( eta_nodes, 6 );

          for (std::size_t j=0; j < MB + 1; ++j )
        	{
        		Base_soln( j, UB )      =   base.solution()( j, fd );
            Base_soln( j, UBd )     =   base.solution()( j, fdd );
            Base_soln( j, PhiB )    =   base.solution()( j, f );
            Base_soln( j, ThetaB )  =   ( 1.0 - BETA ) * base.solution()( j, fdd );
            Base_soln( j, ThetaBd ) =   ( 1.0 - BETA ) * ( - base.solution()( j, f ) *
                                        base.solution()( j, fdd ) - BETA * ( 1.0 -
                                        base.solution()( j, fd ) * base.solution()( j, fd ) ) );
            Base_soln( j, PsiB )    =   ( 1.0 - BETA ) * base.solution()( j, fd );
        	}

          BASE_SOLUTION = Base_soln;
          std::cout << "  * Base flow: 2D Falkner-Skan with transpiration" << std::endl;

        }

        if ( BASE_FLOW=="3D" )
        {
          equation_3D equation;
          plate_BC_3D plate_BC;
          far_BC_3D far_BC;
          equation.beta = 0.1;
          plate_BC.KB = 0.0;
          ODE_BVP<double> base( &equation, eta_nodes, &plate_BC, &far_BC );

          for (std::size_t j=0; j < MB + 1; ++j )
        	{
        		double eta = eta_nodes[ j ];
        		base.solution()( j, f )  	= eta + exp( -eta );
            base.solution()( j, fd ) 	= 1.0 - exp( -eta );
        		base.solution()( j, fdd )  = exp( -eta );
            base.solution()( j, g )  	= 0.35 * (1.0 - exp( -eta ));
            base.solution()( j, gd ) 	= 1 - exp( -eta ) - exp( -1 / (eta * eta) );
        		base.solution()( j, gdd )  = exp( -eta ) - 0.5 * tanh( eta ) + 0.5 * tanh( eta - 2.0 );
        	}

          // Solve the system with KB = 0 then arc-length continue until KB = Param::KB
          double arc_step( 0.01 );
          double max_arc_step( 0.1 );
          base.init_arc( &plate_BC.KB, arc_step, max_arc_step );
          do
          {
            arc_step = base.arclength_solve( arc_step );
          }while( plate_BC.KB < KB );
          plate_BC.KB = KB;
          base.solve_bvp();                               // Solve once more with KB = KB

          // Solve the system with beta = 0.1 then arc-length continue until beta = BETA
          arc_step = -0.01;
          if ( BETA >= 0.1 ) { arc_step = 0.01; }

          base.init_arc( &equation.beta, arc_step, max_arc_step );
          do
          {
            arc_step = base.arclength_solve( arc_step );
          }while( equation.beta < BETA );
          equation.beta = BETA;
          base.solve_bvp();

          OneD_node_mesh<double> Base_soln( eta_nodes, 6 );

          for (std::size_t j=0; j < MB + 1; ++j )
        	{
        		Base_soln( j, UB )      =   base.solution()( j, fd );
            Base_soln( j, UBd )     =   base.solution()( j, fdd );
            Base_soln( j, PhiB )    =   base.solution()( j, f )
                                      + ( 2.0 - BETA ) * base.solution()( j, g );
            Base_soln( j, ThetaB )  =   ( 1.0 - BETA ) * base.solution()( j, fdd )
                                      - ( 2.0 - BETA ) * base.solution()( j, gdd );
            Base_soln( j, ThetaBd ) =   ( 1.0 - BETA ) * ( -(base.solution()( j, f ) +
                                        (2.0 - BETA ) * base.solution()( j, g )) *
                                        base.solution()( j, fdd ) - BETA * ( 1.0 -
                                        base.solution()( j, fd ) * base.solution()( j, fd ) ) )
                                      - ( 2.0 - BETA ) * ( -(base.solution()( j, f ) +
                                        ( 2.0 - BETA ) * base.solution()( j, g )) *
                                        base.solution()( j, gdd ) - BETA * ( 1.0 -
                                        base.solution()( j, gd ) * base.solution()( j, gd ) ) -
                                        2.0 * (1.0 - BETA ) * (base.solution()( j, fd ) -
                                        base.solution()( j, gd )) * base.solution()( j, gd ) );
            Base_soln( j, PsiB )    =   ( 1.0 - BETA ) * base.solution()( j, fd )
                                      - ( 2.0 - BETA ) * base.solution()( j, gd );
        	}

          BASE_SOLUTION = Base_soln;
          std::cout << "  * Base flow: 3D alternative with transpiration" << std::endl;
        }

        std::cout << "  * This number should be zero for the 2D ODE solution and"
                  << " non-zero for the 3D solution: "
                  << ( 1. - BETA ) * BASE_SOLUTION.integral2(UB)
                                   - BASE_SOLUTION.integral2(PsiB) << std::endl;
        std::cout << "  * UB'(eta=0) =" << BASE_SOLUTION( 0, UBd ) << std::endl;
        std::cout << "  * We have solved the ODE problem." << std::endl;
      }

    public:

      /* ----- Constructors and destructors ----- */

      /// Constructor
      SelfSimInjection( ) : HZETA_RIGHT( 16.0 ), ETA_TOP( 128.0 ), N( 200 ),
                            M( 200 ), MB( 20000 ), BETA( 0.0 ), KB( 0.0 ),
                            ZETA0( 1.0 ), K( 0.0 ), MESH( "UNIFORM" ),
                            BASE_FLOW( "2D" ), SPEED_UP( false ), SOLVED( false ),
                            SIGMA( 0.0 ), RX( 1.0 ), FORCING( false ), OUTPUT( true )
      {
        a1 = 60.0;//10.0;//0.1;
        a2 = 4.0;//0.5;
        b1 = 60.0;//10.0;//0.3;
        b2 = 4.0;//0.3;
        GAMMA = 20.0;
        N_TRANSP = 1;
      }

      /// Destructor
	   	~SelfSimInjection(){ }

      /* ----- Methods for accessing parameters ----- */

      /// Return a handle to the size of the domain in the zeta_hat direction
      double& hzeta_right() { return HZETA_RIGHT; }

      /// Return a handle to the size of the domain in the eta direction
      double& eta_top() { return ETA_TOP; }

      /// Return a handle to the number of intervals in the zeta_hat direction
      std::size_t& hzeta_intervals() { return N; }

      /// Return a handle to the number of intervals in the eta direction
      std::size_t& eta_intervals() { return M; }

      /// Return a handle to the number of eta intervals in the base flow ODE
      std::size_t& base_intervals() { return MB; }

      /// Return a handle to the Hartree parameter
      double& hartree() { return BETA; }

      /// Return a handle to the base flow injection
      double& base_injection() { return KB; }

      /// Return a handle to the injection width
      double& injection_width() { return ZETA0; }

      /// Return a handle to the injection parameter
      double& injection() { return K; }

      /// Return a handle to the wave amplitude
      double& wave_amplitude() { return SIGMA; }

      /// Return a handle to the local Reynold's Number
      double& local_Reynolds() { return RX; }

      /// Return a handle to the wavenumber
      double& wavenumber() { return ALPHA; }

      /// Set whether to output the solution
      void set_output( bool output ){ OUTPUT = output; }

      /// Set the transverse wave velocity
      void set_v_wave( TwoD_node_mesh< std::complex<double> > v_wave )
      {
        v_WAVE = v_wave;
        //v_WAVE_CONJ = v_wave.conjugate();
      }

      /// Set the spanwise wave velocity
      void set_w_wave( TwoD_node_mesh< std::complex<double> > w_wave )
      {
        w_WAVE = w_wave;
        //w_WAVE_CONJ = w_wave.conjugate();
      }

      /// Turn forcing on and off
      void forcing( bool f ){ FORCING = f; }

      /// Return the mesh definition string
      std::string mesh() { return MESH; }

      /// Method to set the SOLVED variable (useful if reading solution from file)
      void set_solved( bool solved ){ SOLVED = solved; }

      /// Return whether the equations have been solved or not
      bool solved() { return SOLVED; }

      /// Method to set the mesh definition string
      void set_mesh( std::string mesh )
      {
        if ( mesh == "UNIFORM" || mesh == "NONUNIFORM" )
        {
          MESH = mesh;
        }
        else
        {
          throw Error( "Definition is not recognised in set_mesh method." );
        }
      }

      /// Return the base flow definition string
      std::string base_flow() { return BASE_FLOW; }

      /// Return the output path string
      std::string output_path() { return OUTPUT_PATH; }

      /// Method to set the base flow definition string
      void set_base_flow( std::string base_flow )
      {
        if ( base_flow == "2D" || base_flow == "3D" )
        {
          BASE_FLOW = base_flow;
        }
        else
        {
          throw Error( "Definition is not recognised in set_base_flow method." );
        }
      }

      /// Method for setting whether to reuse the factorised Matrix
      void speed_up( bool speed_up ){ SPEED_UP = speed_up; }

      /// Return the base flow solution
      OneD_node_mesh<double> base_flow_solution() { return BASE_SOLUTION; }

      /// Return BASE_ETA_NODES
      Vector<double> base_eta_nodes(){ return BASE_ETA_NODES; }

      /// Return the ETA_NODES
      Vector<double> eta_nodes(){ return ETA_NODES; }

      /// Return the HZETA_NODES
      Vector<double> hzeta_nodes(){ return HZETA_NODES; }

      /// Return the X_NODES
      Vector<double> x_nodes(){ return X_NODES; }

      /// Return the Y_NODES
      Vector<double> y_nodes(){ return Y_NODES; }

      /* ----- Mesh functions ----- */

      double mesh_X( const double& zeta )
      {
        if ( MESH=="UNIFORM" ) { return zeta; }
        //if ( MESH=="NONUNIFORM" ) { return std::pow(zeta + a1, a2); }
        if ( MESH=="NONUNIFORM" ) { return a1 + zeta - a1 * std::exp( - zeta / a2 ); }
        else { throw Error( "Mesh string is not defined." ); }
      }

      double mesh_Xd( const double& zeta )
      {
        if ( MESH=="UNIFORM" ) { return 1.0; }
        //if ( MESH=="NONUNIFORM" ) { return a2 * std::pow(zeta + a1, a2 - 1); }
        if ( MESH=="NONUNIFORM" ) { return 1 + ( a1 / a2 ) * std::exp( - zeta / a2 ); }
        else { throw Error( "Mesh string is not defined." ); }
      }

      double mesh_Xdd( const double& zeta )
      {
        if ( MESH=="UNIFORM" ) { return 0.0; }
        //if ( MESH=="NONUNIFORM" ) { return a2 * (a2 - 1) * std::pow(zeta + a1, a2 - 2); }
        if ( MESH=="NONUNIFORM" ) { return - ( a1 / ( a2 * a2 ) ) * std::exp( - zeta / a2 ); }
        else { throw Error( "Mesh string is not defined." ); }
      }

      double mesh_Y( const double& eta )
      {
        if ( MESH=="UNIFORM" ) { return eta; }
        //if ( MESH=="NONUNIFORM" ) { return std::pow(eta + b1, b2); }
        if ( MESH=="NONUNIFORM" ) { return b1 + eta - b1 * std::exp( - eta / b2 ); }
        else { throw Error( "Mesh string is not defined." ); }
      }

      double mesh_Yd( const double& eta )
      {
        if ( MESH=="UNIFORM" ) { return 1.0; }
        //if ( MESH=="NONUNIFORM" ) { return b2 * std::pow(eta + b1, b2 - 1); }
        if ( MESH=="NONUNIFORM" ) { return 1 + ( b1 / b2 ) * std::exp( - eta / b2 ); }
        else { throw Error( "Mesh string is not defined." ); }
      }

      double mesh_Ydd( const double& eta )
      {
        if ( MESH=="UNIFORM" ) { return 0.0; }
        //if ( MESH=="NONUNIFORM" ) { return b2 * (b2 - 1) * std::pow(eta + b1, b2 - 2); }
        if ( MESH=="NONUNIFORM" ) { return - ( b1 / ( b2 * b2 ) ) * std::exp( - eta / b2 ); }
        else { throw Error( "Mesh string is not defined." ); }
      }

      class invert_eta : public Residual<double>
      {
        // USED TO INVERT THE NON-UNIFORM MESH MAPPING
        public:
        double Y0;
        std::string mesh;

        invert_eta() : Residual<double>( 1 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &f ) const
        {
          SelfSimInjection SSI;
          if ( mesh=="UNIFORM" ){ SSI.set_mesh( "UNIFORM" ); }
          if ( mesh=="NONUNIFORM" ){ SSI.set_mesh( "NONUNIFORM" ); }
          f[ 0 ] = SSI.mesh_Y( z[0] ) - Y0;
        }
      };

      class invert_zeta : public Residual<double>
      {
        // USED TO INVERT THE NON-UNIFORM MESH MAPPING
        public:
        double X0;
        std::string mesh;

        invert_zeta() : Residual<double>( 1 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &f ) const
        {
          SelfSimInjection SSI;
          if ( mesh=="UNIFORM" ){ SSI.set_mesh( "UNIFORM" ); }
          if ( mesh=="NONUNIFORM" ){ SSI.set_mesh( "NONUNIFORM" ); }
          f[ 0 ] = SSI.mesh_X( z[0] ) - X0;
        }
      };

      /* ----- Setup the mesh ----- */

      void mesh_setup()
      {
        // Define the remapped (non-uniform mesh) domain
        double left   = mesh_X( 0.0 );
        double right  = mesh_X( HZETA_RIGHT );
        double bottom = mesh_Y( 0.0 );
        double top    = mesh_Y( ETA_TOP );

        // Base flow ODE eta nodes
        BASE_ETA_NODES.linspace( 0.0, ETA_TOP, MB + 1 );

        // Nodal positions in the remapped domain (spanned by X,Y)
        X_NODES.linspace( left, right, N + 1 );
        Y_NODES.linspace( bottom, top, M + 1 );

        // Vectors for writing data on the original zeta-eta domain
        HZETA_NODES.linspace( 0.0, HZETA_RIGHT, N + 1 );
        ETA_NODES.linspace( 0.0, ETA_TOP, M + 1 );

        // To find eta=eta(Y) and zeta=zeta(X) we will use Newton iteration
        invert_eta find_eta;
        find_eta.mesh = MESH;
        Newton<double> newton_eta( &find_eta );
        for ( unsigned j = 0; j < M + 1; ++j )
        {
          unsigned kmin(0); double min(99);
          for ( unsigned k = 0; k < M + 1; ++k )
          {
            if ( std::abs( mesh_Y( ETA_NODES[k] ) - Y_NODES[j] ) < min )
            {
              min = std::abs( mesh_Y( ETA_NODES[k] ) - Y_NODES[j] );
              kmin = k;
            }
          }
          find_eta.Y0 = Y_NODES[ j ];
          Vector<double> guess( 1, 1.0 );
          guess[ 0 ] = ETA_NODES[ kmin ];
          newton_eta.iterate( guess );
          ETA_NODES[j] = guess[ 0 ];
        }

        invert_zeta find_zeta;
        find_zeta.mesh = MESH;
        Newton<double> newton_zeta( &find_zeta );
        //
        for ( unsigned i = 0; i < N + 1; ++i )
        {
          unsigned kmin(0); double min(99);
          for ( unsigned k = 0; k < N + 1; ++k )
          {
            if ( std::abs( mesh_X( HZETA_NODES[k] ) - X_NODES[i] ) < min )
            {
              min = std::abs( mesh_X( HZETA_NODES[k] ) - X_NODES[i] );
              kmin = k;
            }
          }
          find_zeta.X0 = X_NODES[ i ];
          Vector<double> guess( 1, 1.0 );
          guess[ 0 ] = HZETA_NODES[ kmin ];
          newton_zeta.iterate( guess );
          HZETA_NODES[ i ] = guess[ 0 ];
        }

      }

      /* ----- Methods for solving the perturbation problem ----- */

      std::size_t col( const std::size_t& i, const std::size_t& j, const std::size_t& k )
      {
        // Return the column number for the kth variable at node (i,j)
        return 4 * ( i * ( M + 1 ) + j ) + k;
      }

      // Virtual function for defining the transpiration function
      virtual double Phi_w_func( const double& hzeta )
      {
        throw Error( "Phi_w function is not defined" );
      }

      double Phi_w_function( const double& hzeta )
      {
        // Return the transpiration function
        /*if ( N_TRANSP < 1 )
        {
          return 0.0;
        }
        else
        {
          double sum( 0.0 );
          int sign;
          for (std::size_t i=1; i<N_TRANSP; ++i)
          {
            sign = i % 2 ? -1 : 1; // equivalent to (-1)^i since i % 2 = 0 = false and i % 2 = 1 = true
            sum += sign * tanh( GAMMA * ( hzeta - ((1.*i)/N_TRANSP) ) );
          }
          sign = N_TRANSP % 2 ? -1 : 1; // (-1)^N
          return - K * 0.5 *( 1 + 2 * sum + sign * tanh( GAMMA * ( hzeta - 1. ) ) );
        }*/

        /*return - Param::K * 0.5 * ( tanh( Param::gamma * ( hzeta - 1. ) )
               - tanh( Param::gamma * ( hzeta - 2. ) ) );
        */

        return - K * exp( - hzeta * hzeta / 4 ) / sqrt( M_PI );
      }

      double Phi_w_hzeta_function( const double& hzeta )
      {
        // Return derivative of transpiration wrt hzeta
        /*if ( N_TRANSP < 1 )
        {
          return 0.0;
        }
        else
        {
          double sum( 0.0 );
          double sech_squared;
          int sign;
          for (std::size_t i=1; i<N_TRANSP; ++i)
          {
            sign = i % 2 ? -1 : 1; // equivalent to (-1)^i
            sech_squared = pow( cosh( GAMMA * ( hzeta - ((1.*i)/N_TRANSP) ) ) , -2. );
            sum += sign * sech_squared;
          }
          sign = N_TRANSP % 2 ? -1 : 1; // (-1)^N
          sech_squared = pow( cosh( GAMMA * ( hzeta - 1. ) ) , -2. );
          return - K * 0.5 * GAMMA * ( 2 * sum + sign * sech_squared );
        }*/

        /*double sech_squared = pow( cosh( Param::gamma * ( hzeta - 1. ) ) , -2. );
        double sech_squared_2 = pow( cosh( Param::gamma * ( hzeta - 2. ) ) , -2. );
        return - Param::K * 0.5 * Param::gamma * ( sech_squared - sech_squared_2  );*/

        return 0.5 * hzeta * K * exp( - hzeta * hzeta / 4 ) / sqrt( M_PI );
      }

      void solve_perturbation_eqns()
      {
        std::cout << "*** Solving the perturbation equations ***" << std::endl;
        int Phi = static_cast<int>(SelfSim::Phi);
        int Psi = static_cast<int>(SelfSim::Psi);
        int U = static_cast<int>(SelfSim::U);
        int Theta = static_cast<int>(SelfSim::Theta);

        TwoD_node_mesh<double> q( X_NODES, Y_NODES, 4 );
        TwoD_node_mesh<double> q_output( HZETA_NODES, ETA_NODES, 8 );
        Q = q;
        Q_output = q_output;

        // step sizes in the remapped domain : these should be constants
        const double dY( Y_NODES[ 1 ] - Y_NODES[ 0 ] );
        const double dX( X_NODES[ 1 ] - X_NODES[ 0 ] );

        // Vector for the RHS of the matrix problem
        Vector<double> B( 4 * ( M + 1 ) * ( N + 1 ), 0.0 );
        //TODO use Q to initialise B ? (only if forcing?)

        /* Iterate to a solution */
        double max_residual( 0.0 );             // Maximum residual
        std::size_t iteration( 0 );             // Initialise iteration counter
        std::size_t max_iterations( 20 );       // Maximum number of iterations
        if( SPEED_UP ){ max_iterations = 200; }
        // Eigen objects (only used if SPEED_UP=true)
        Eigen::SparseMatrix<double, Eigen::ColMajor, long long> A_Eigen( 4 * ( M + 1 ) * ( N + 1 ), 4 * ( M + 1 ) * ( N + 1 ) );
        Eigen::SparseLU< Eigen::SparseMatrix<double, Eigen::ColMajor, long long> > solver;

        do {

          SparseMatrix<double> A( 4 * ( M + 1 ) * ( N + 1 ), 4 * ( M + 1 ) * ( N + 1 ) );
          std::cout << "  * Assembling sparse matrix problem" << std::endl;

          Timer timer;                                        // Timer
          timer.start();

          std::size_t row( 0 );                               // Initialise row counter

          /* hzeta = 0 boundary ( left boundary ) */
          std::size_t i( 0 );

          for ( std::size_t j = 0; j < M + 1 ; ++j )
          {
            double hzeta( HZETA_NODES[ 0 ] );
            double Xd( mesh_Xd( hzeta ) );
            double eta( ETA_NODES[ j ] );
            Vector<double> Base( BASE_SOLUTION.get_interpolated_vars( eta ) );

            // Phi_hzeta = 0
            A( row, col( i, j, Phi ) )      = -3.*Xd/(2*dX);
            A( row, col( i + 1, j, Phi ) )  =  4.*Xd/(2*dX);
            A( row, col( i + 2, j, Phi ) )  = -1.*Xd/(2*dX);

            B[ row ]                        = -( Xd*( -3*Q(i,j,Phi) + 4*Q(i+1,j,Phi)
                                              -Q(i+2,j,Phi) )/(2*dX) );
            ++row;

            // Psi = 0
            A( row, col( i, j, Psi ) )      =   1;

            B[ row ]                        = - Q( i, j, Psi );
            ++row;

            // U_hzeta = 0
            A( row, col( i, j, U ) )        = -3.*Xd/(2*dX);
            A( row, col( i + 1, j, U ) )    =  4.*Xd/(2*dX);
            A( row, col( i + 2, j, U ) )    = -1.*Xd/(2*dX);

            B[ row ]                        = -( Xd*( -3*Q(i,j,U) + 4*Q(i+1,j,U)
                                              -Q(i+2,j,U) )/(2*dX) );
            ++row;

            // Theta = 0
            A( row, col( i, j, Theta ) )    =   1;
            
            B[ row ]                        = - Q( i, j, Theta );
            ++row;

          } // end for loop over LHS eta nodes

          /* Interior points between the hzeta boundaries */
          for ( std::size_t i = 1; i < N; ++i )
          {
            // hzeta location
            double hzeta( HZETA_NODES[ i ] );
            double Xd( mesh_Xd( hzeta ) );
            double Xdd( mesh_Xdd( hzeta ) );
            // Wall transpiration
            //double Phi_w( Phi_w_function( hzeta ) );
            double Phi_w( Phi_w_func( hzeta ) );
            //double Phi_w_hzeta( Phi_w_hzeta_function( hzeta ) );

            /* eta = 0 boundary ( bottom boundary ) */
            std::size_t j( 0 );
            double eta( ETA_NODES[ j ] );
            double Yd( mesh_Yd( eta ) );

            // Phi = Phi_w
            A( row, col( i, j, Phi ) )        =  1.;
            B[ row ]                          = -Q( i, j, Phi ) + Phi_w;
            ++row;
            // Psi = 0
            A( row, col( i, j, Psi ) )        =  1.;
            B[ row ]                          = -Q( i, j, Psi );
            ++row;
            // U = 0
            A( row, col( i, j, U ) )          =  1.;
            B[ row ]                          = -Q( i, j, U );
            ++row;
            // Theta - Psi_eta = -( 1 / ( zeta0^2 ) ) * Phi_w_hzeta
            /*A( row, col( i, j, Theta ) )      =  1.;
            A( row, col( i, j, Psi ) )        =  3.*Yd / ( 2 * dY );
            A( row, col( i, j + 1, Psi ) )    = -4.*Yd / ( 2 * dY );
            A( row, col( i, j + 2, Psi ) )    =  1.*Yd / ( 2 * dY );
            B[ row ]                          = -Q(i,j,Theta) + Yd*( -3*Q(i,j,Psi)
                                                + 4*Q(i,j+1,Psi) - Q(i,j+2,Psi) ) / (2*dY)
                                                - ( 1. / ( ZETA0 * ZETA0 ) )
                                                * Phi_w_hzeta;*/
            // theta = psi_eta  - phi_zeta on eta=0
            A( row, col( i, j, Theta ) ) = 1;
            A( row, col( i, j, Psi ) ) = Yd*3./(2*dY);
            A( row, col( i, j + 1, Psi ) ) = -Yd*4./(2*dY);
            A( row, col( i, j + 2, Psi ) ) = Yd/(2*dY);
            A( row, col( i + 1, 0, Phi ) ) = Xd/(2*dX);
            A( row, col( i - 1, 0, Phi ) ) = -Xd/(2*dX);
            B[ row ] = -Q( i, 0, Theta ) + Yd*( -3*Q(i,0,Psi) + 4*Q(i,1,Psi) - Q(i,2,Psi) ) / (2*dY)
              - ( Q(i+1,0,Phi) -  Q(i-1,0,Phi) )*Xd/(2*dX);

            ++row;

            /* Main interior grid points */
            for ( std::size_t j = 1; j < M; ++j )
            {
              // eta location
              double eta( ETA_NODES[ j ] );
              double Yd( mesh_Yd( eta ) );
              double Ydd( mesh_Ydd( eta ) );
              // Base solution
              Vector<double> Base( BASE_SOLUTION.get_interpolated_vars( eta ) );

              // Laplacian coefficients for finite-differencing

              // X(i,j-1)
              double laplace_1 =  ( Yd*Yd/(dY*dY) - Ydd/ (2.*dY) ) ;
              // X(i-1,j)
              double laplace_3 = ( Xd*Xd/(dX*dX) - Xdd/(2.*dX) ) / ( ZETA0 * ZETA0 );
              // X(i,j)
              double laplace_4 = -2.*( Yd*Yd / (dY*dY)
                                 + Xd*Xd/( ZETA0 * ZETA0 * dX * dX ) );
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
              A( row, col( i, j + 1, U ) )        = -( 2. - BETA )*Yd/( 2 * dY );
              A( row, col( i, j - 1, U ) )        =  ( 2. - BETA )*Yd/( 2 * dY );
              // Theta_hzeta
              A( row, col( i + 1, j, Theta ) )    =  Xd / ( 2 * dX );
              A( row, col( i - 1, j, Theta ) )    = -Xd / ( 2 * dX );

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
              A( row, col( i, j, Phi ) )          =    Base[ UBd ] + Guess_eta[ U ];

              // Forcing 1
              double F_1;
              F_1 = calculate_forcing_1( i, j );

              // Residual
              B[ row ]        = - Guess_laplace[ U ]
                                + BETA * ( 2. * Base[ UB ] + Guess[ U ] ) * Guess[ U ]
                                - ( hzeta * Base[ PsiB ] + Guess[ Psi ] )
                                * ( Guess_hzeta[ U ] )
                                - ( Base[ PhiB ] + Guess[ Phi ] ) * Guess_eta[ U ]
                                - Base[UBd] * Guess[Phi]
                                + std::pow( RX, -2.0/3.0 ) * SIGMA * SIGMA * F_1;
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
              A( row, col( i, j, U ) )            +=   ( 2. - BETA )
                                                     * ( hzeta * Base[ Theta ]
                                                     + Guess[ Theta ] );

              // Forcing 2
              double F_2;
              F_2 = calculate_forcing_2( i, j );

              // Residual
              B[ row ]      = - Guess_laplace[ Theta ]
                              + 2.*( 1. - BETA )
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
                              + std::pow( RX, -1.0/6.0 ) * SIGMA * SIGMA * F_2;
              ++row;
            }

            /* eta = eta_inf boundary ( top boundary ) */
            j = M ;
            eta = ETA_NODES[ j ];
            Yd = mesh_Yd( eta );

            const double rsq =  eta * eta + ZETA0 * ZETA0 * hzeta * hzeta ;
            // Phi_eta*( eta^2 + zeta_0^2*hzeta^2) + [ 2*eta - (eta^2 + zeta_0^2*hzeta^2)/eta ]*Phi = 0
            A( row, col( i, j, Phi ) )        =   3.0 * Yd * rsq / (2*dY);
            A( row, col( i, j - 1, Phi ) )    = - 4.0 * Yd * rsq / (2*dY);
            A( row, col( i, j - 2, Phi ) )    =   1.0 * Yd * rsq / (2*dY);
            A( row, col( i, j, Phi ) )       +=   2.0 * eta - ( rsq / eta );

            B[ row ]                          = - ( 3 * Q( i, j, Phi ) - 4 * Q( i, j-1, Phi )
                                                + Q( i, j-2, Phi ) ) * Yd * rsq / ( 2 * dY )
                                                + ( ( rsq / eta ) - 2.0 * eta ) * Q( i, j, Phi );
            ++row;

            // Psi_eta*( eta^2 + zeta_0^2*hzeta^2) + 2 * eta * Psi = 0
            A( row, col( i, j, Psi ) )        =   3.0 * Yd * rsq / (2*dY);
            A( row, col( i, j - 1, Psi ) )    = - 4.0 * Yd * rsq / (2*dY);
            A( row, col( i, j - 2, Psi ) )    =   1.0 * Yd * rsq / (2*dY);
            A( row, col( i, j, Psi ) )       +=   2 * eta;

            B[ row ]                          =  - (( 3 * Q( i, j, Psi ) - 4 * Q( i, j-1, Psi )
                                                 + Q( i, j-2, Psi ) ) * Yd * rsq / ( 2 * dY ))
                                                 - 2 * eta * Q( i, j, Psi );
            ++row;

            // U = 0
            A( row, col( i, j, U ) )            =   1;
            B[ row ]                            = - ( Q( i, j, U ) );
            ++row;

            // Theta = 0
            A( row, col( i, j, Theta ) )        =   1;
            B[ row ]                            = - ( Q( i, j, Theta ) );
            ++row;

          } // End of for loop over interior nodes

          /* hzeta = hzeta_inf boundary ( right boundary ) */

          for ( std::size_t j = 0; j < M + 1; ++j )
          {
            //offset for global problem
            std::size_t i( N );
            double hzeta( HZETA_NODES[ i ] );
            double Xd( mesh_Xd( hzeta ) );
            double eta( ETA_NODES[ j ] );


            Vector<double> Base( BASE_SOLUTION.get_interpolated_vars( eta ) );

            const double rsq =  eta*eta+hzeta*hzeta;

            // (eta^2 + zeta_0^2 * hzeta^2) * Phi_hzeta + 2 * zeta_0^2 * hzeta * Phi = 0
            A( row, col( i, j, Phi ) )          =   3.0 * Xd * rsq / ( 2 * dX );
            A( row, col( i - 1, j, Phi ) )      = - 4.0 * Xd * rsq / ( 2 * dX );
            A( row, col( i - 2, j, Phi ) )      =   1.0 * Xd * rsq / ( 2 * dX );
            A( row, col( i, j, Phi ) )         +=   2 * ZETA0 * ZETA0 * hzeta;

            B[ row ]        = - rsq * ( 3 * Q( i, j, Phi) - 4 * Q( i - 1, j, Phi)
                              + Q( i - 2, j, Phi) ) * Xd / ( 2 * dX )
                              - 2 * ZETA0 * ZETA0 * hzeta * Q( i, j, Phi );
            ++row;

            // (eta^2 + zeta_0^2 * hzeta^2)*Psi_hzeta + (2*zeta_0^2*hzeta-(eta^2 + zeta_0^2*hzeta^2)/hzeta)*Psi = 0
            A( row, col( i, j, Psi ) )          =   3.0 * Xd * rsq / ( 2 * dX );
            A( row, col( i - 1, j, Psi ) )      = - 4.0 * Xd * rsq / ( 2 * dX );
            A( row, col( i - 2, j, Psi ) )      =   1.0 * Xd * rsq / ( 2 * dX );
            A( row, col( i, j, Psi ) )         +=   2 * ZETA0 * ZETA0 * hzeta - (rsq / hzeta);


            B[ row ]        = - (rsq * ( 3 * Q( i, j, Psi ) - 4 * Q( i - 1, j, Psi )
                              + Q( i - 2, j, Psi) ) * Xd / ( 2 * dX ))
                              - 2 * ZETA0 * ZETA0 * hzeta  * Q( i, j, Psi)
                              + (rsq / hzeta)  * Q( i, j, Psi);
            ++row;

            // hzeta * U_hzeta + 2 * U = 0
            A( row, col( i, j, U ) )            =   hzeta * 3. * Xd / ( 2 * dX ) + 2.;
            A( row, col( i - 1, j, U ) )        = - hzeta * 4. * Xd / ( 2 * dX );
            A( row, col( i - 2, j, U ) )        =   hzeta * 1. * Xd / ( 2 * dX );

            B[ row  ]       = - hzeta * ( 3 * Q( i, j, U ) - 4 * Q( i - 1, j, U )
                              + Q( i - 2, j, U) ) * Xd / ( 2 * dX ) - 2 * Q( i, j, U );
            ++row;

            // hzeta * Theta_hzeta + Theta = 0
            A( row, col( i, j, Theta )  )       =   hzeta * 3. * Xd / ( 2 * dX ) + 1.;
            A( row, col( i - 1, j, Theta ) )    = - hzeta * 4. * Xd / ( 2 * dX );
            A( row, col( i - 2, j, Theta ) )    =   hzeta * 1. * Xd / ( 2 * dX );

            B[ row ]        = - hzeta * ( 3 * Q( i, j, Theta ) - 4 * Q( i - 1, j, Theta )
                              + Q( i - 2, j, Theta ) ) * Xd / ( 2 * dX ) - Q( i, j, Theta ) ;
            ++row;

          } // End of loop over nodes

          max_residual = B.norm_inf();
          std::cout << "***                                              Maximum residual = "
                    << B.norm_inf() << std::endl;

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
              Q( i, j, Phi )    += B[ col( i, j, Phi ) ];
              Q( i, j, Psi )    += B[ col( i, j, Psi ) ];
              Q( i, j, U )      += B[ col( i, j, U ) ];
              Q( i, j, Theta )  += B[ col( i, j, Theta ) ];
            }
          }

          std::cout << "***    Iteration = " << iteration
               << "    Maximum correction = " << B.norm_inf() << std::endl;

          ++iteration;
        }while( ( max_residual > 1.e-8 ) && ( iteration < max_iterations ) ); // End iteration

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
            // first 4 values output are the without the underlying base flow
            Q_output( i, j, 0 ) = Q( i, j, Phi);
            Q_output( i, j, 1 ) = Q( i, j, Psi);
            Q_output( i, j, 2 ) = Q( i, j, U);
            Q_output( i, j, 3 ) = Q( i, j, Theta);
            // second 4 values are the "full" solution, but still with the zeta0 scaling
            Q_output( i, j, 4 ) =   Q( i, j, Phi)
                                + BASE_SOLUTION.get_interpolated_vars( eta )[PhiB];
            Q_output( i, j, 5 ) = Q( i, j, Psi)
                                + hzeta * BASE_SOLUTION.get_interpolated_vars( eta )[PsiB];
            Q_output( i, j, 6 ) =   Q( i, j, U)
                                + BASE_SOLUTION.get_interpolated_vars( eta )[UB];
            Q_output( i, j, 7 ) = Q( i, j, Theta)
                                + hzeta * BASE_SOLUTION.get_interpolated_vars( eta )[ThetaB];
          }
        }

        // Get the wall shear values as a function of hzeta
        OneD_node_mesh<double> wall_shear( HZETA_NODES, 1 );
        for ( std::size_t i=0; i < N + 1; ++i )
        {
          wall_shear( i, 0 ) = -( 3 * Q_output(i,0,U+4) - 4 * Q_output(i,1,U+4)
                + Q_output(i,2,U+4) ) * mesh_Yd(0.0)/(2*dY);
        }
        WALL_SHEAR = wall_shear;

        A_PARAM = Q_output( 0, M, Phi ) * ETA_TOP;
        U_ETA = -( 3 * Q_output(0,0,U+4) - 4 * Q_output(0,1,U+4)
                + Q_output(0,2,U+4) ) * mesh_Yd(0.0)/(2*dY);

        // Find value of eta on zeta=0 at which U=1/2
        std::size_t lower = 0;
        std::size_t upper = 1;
        for (std::size_t j=0; j < M; ++j)
        {
          if ( Q_output(0,j,U+4) < 0.5 && Q_output(0,j+1,U+4) > 0.5 ) { lower = j; upper=j+1; }
        }
        // linearly interpolate
        ETA_HALF =  ( 0.5 - Q_output(0,lower,U+4) ) * ( ETA_NODES[upper] - ETA_NODES[lower] )
                  / ( Q_output(0,upper,U+4) - Q_output(0,lower,U+4)  ) + ETA_NODES[lower];

        // Integral of U^2 over the cross-section
        INT_U_SQUARED = Q_output.square_integral2D( U );

        SOLVED = true;

      } // End of solve_perturbation_eqns function

      /* ----- Output solution data ----- */

      /// Return the solution mesh
      TwoD_node_mesh<double> solution() {
        if ( SOLVED ){ return Q_output; }
        else { throw Error( "solution() error equations have not been solved." ); }
      }

      /// Return the wall shear solution mesh
      OneD_node_mesh<double> solution_wall_shear() {
        if ( SOLVED ){ return WALL_SHEAR; }
        else { throw Error( "solution_wall_shear() error equations have not been solved." ); }
      }

      /// Reset the solution mesh from an external mesh
      void set_solution( TwoD_node_mesh<double> sol ){
        if ( SOLVED ){ Q_output = sol; }
        else { throw Error( "set_solution() error equations have not been solved." ); }
      }

      /// Reset the base flow solution mesh from an external mesh
      void set_base_solution( OneD_node_mesh<double> base_sol ){
        if ( SOLVED ){ BASE_SOLUTION = base_sol; }
        else { throw Error( "set_base_solution() error equations have not been solved." ); }
      }

      /// Return the mass flux parameter
      double mass_flux() { return A_PARAM; }

      /// Return the shear stress at the origin
      double shear_at_origin() { return U_ETA; }

      /// Return the value of eta on zeta=0 at which U=1/2
      double eta_half() { return ETA_HALF; }

      /// Output the base flow solution to a file
      void output_base_solution() {
        BASE_SOLUTION.output( OUTPUT_PATH + "Base_soln.dat" );
      }

      void output(){
        if ( SOLVED ) {
          // Convert param to string
          std::stringstream ss;
          ss << ZETA0;
          std::string str = ss.str();
          Q_output.dump_gnu( OUTPUT_PATH + "Qout_" + str + ".dat" );
          WALL_SHEAR.output( OUTPUT_PATH + "Wall_shear_" + str + ".dat" );
        }
        else { throw Error( "output() error equations have not been solved." ); }
      }

      void output( std::string param ){
        if ( SOLVED ) {
          // Convert param to string
          std::stringstream ss;
          if ( param == "K" )
          {
            ss << K;
          }
          if ( param == "ZETA0" )
          {
            ss << ZETA0;
          }
          std::string str = ss.str();

          Q_output.dump_gnu( OUTPUT_PATH + "Qout_" + str + ".dat" );
          WALL_SHEAR.output( OUTPUT_PATH + "Wall_shear_" + str + ".dat" );
        }
        else { throw Error( "output() error equations have not been solved." ); }
      }

      void iterate_on_zeta0( const double& step, const double& max )
      {
        std::ostringstream ss;
        ss << "./DATA/K_" << K << "_beta_" << BETA << "_" << N + 1
           << "x" << M + 1 << "_" << HZETA_RIGHT << "_" << ETA_TOP << "/";
        OUTPUT_PATH = ss.str();
        setup();
        output_base_solution();
        TrackerFile metric( OUTPUT_PATH + "A_file.dat" );
        metric.push_ptr( &ZETA0, "zeta0" );
        metric.push_ptr( &A_PARAM, "A" );
        metric.push_ptr( &U_ETA, "U_eta(0,0)");
        metric.push_ptr( &ETA_HALF, "eta at which U=1/2 on zeta=0" );
        metric.push_ptr( &INT_U_SQUARED, "integral U^2 over the cross-section");
        metric.header();

        do {
          solve_perturbation_eqns();
          metric.update();
          output( "ZETA0" );
          std::cout << "  * zeta0 = " << ZETA0 << ", A = " << A_PARAM << std::endl;
          ZETA0 += step;
        }while ( ZETA0 <= max );
      }

      void iterate_on_K( const double& step, const double& max )
      {
        std::ostringstream ss;
        ss << "./DATA/zeta0_" << ZETA0 << "_beta_" << BETA << "_" << N + 1
           << "x" << M + 1 << "_" << HZETA_RIGHT << "_" << ETA_TOP << "/";
        OUTPUT_PATH = ss.str();
        setup();
        output_base_solution();
        TrackerFile metric( OUTPUT_PATH + "A_file.dat" );
        metric.push_ptr( &ZETA0, "zeta0" );
        metric.push_ptr( &K, "K" );
        metric.push_ptr( &A_PARAM, "A" );
        metric.push_ptr( &U_ETA, "U_eta(0,0)");
        metric.push_ptr( &ETA_HALF, "eta at which U=1/2 on zeta=0" );
        metric.push_ptr( &INT_U_SQUARED, "integral U^2 over the cross-section");
        metric.header();

        do {
          solve_perturbation_eqns();
          metric.update();
          output( "K" );
          std::cout << "  * K = " << K << ", A = " << A_PARAM << std::endl;
          K += step;
        }while ( K <= max );
      }

      void set_output_path()
      {
        std::ostringstream ss;
        ss << "./DATA/K_" << K << "_zeta0_" << ZETA0 << "_beta_" << BETA << "_" << N + 1
           << "x" << M + 1 << "_" << HZETA_RIGHT << "_" << ETA_TOP << "/";
        OUTPUT_PATH = ss.str();
      }

      void solve()
      {
        if ( OUTPUT )
        {
          set_output_path();
        }
        setup();
        solve_perturbation_eqns();
      }

  }; // End of class SelfSimInjection
} // End of namespace TSL


#endif
