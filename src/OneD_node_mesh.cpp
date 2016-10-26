/* OneD_node_mesh - A class specifying a one dimensional mesh object used
                    for storing and manipulating data.
*/

#include "OneD_node_mesh.h"

namespace TSL
{

	template <>
  Vector<double> OneD_node_mesh<double, double>::get_interpolated_vars( const double& x_pos ) const
  {
    for ( unsigned node = 0; node < NODES.size() - 1; ++node )
    {
      // find bracketing nodes - incl shameless hack for evaluations at the boundary
      if ( ( ( NODES[ node ] < x_pos ) && ( NODES[ node + 1 ] > x_pos ) )   
              ||  ( std::abs( NODES[ node ] - x_pos ) < 1.e-7  ) || ( std::abs( NODES[ node + 1 ] - x_pos ) < 1.e-7 ) )
      {
        // distance from left node
        double delta_x( x_pos - NODES[ node ] );
        // empty data to return
        Vector<double> left;
        Vector<double> right;
        Vector<double> deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        deriv = ( right - left ) / ( NODES[ node + 1 ] - NODES[ node ] );
        // overwrite right
        right = left + deriv * delta_x;
        return right;
      }
    }
    // If the position is not in the range throw an error
    throw Error( "Mesh error: interpolation not in range " );
  }

  template <>
  Vector<std::complex<double> > OneD_node_mesh<std::complex<double>, double>::get_interpolated_vars( const double& x_pos ) const
  {
    for ( unsigned node = 0; node < NODES.size() - 1; ++node )
    {
      // find bracketing nodes - incl shameless hack for evaluations at the boundary
      if ( ( NODES[ node ] < x_pos  || std::abs( NODES[ node ] - x_pos ) < 1.e-7 ) &&
           ( NODES[ node + 1 ] > x_pos || std::abs( NODES[ node + 1 ] - x_pos ) < 1.e-7 ) )
      {
        // distance from left node
        double delta_x( x_pos - NODES[ node ] );
        // empty data to return
        Vector<std::complex<double> > left;
        Vector<std::complex<double> > right;
        Vector<std::complex<double> > deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        deriv = ( right - left ) / ( NODES[ node + 1 ] - NODES[ node ] );
        // overwrite right
        right = left + deriv * delta_x;
        return right;
      }
    }
    // If the position is not in the range throw an error
    throw Error( "Mesh error: interpolation not in range " );
  }

  template <>
  Vector<std::complex<double> > OneD_node_mesh<std::complex<double>, std::complex<double> >::get_interpolated_vars( const std::complex<double>& pos ) const
  {
    double x_pos( pos.real() );
    for ( unsigned node = 0; node < NODES.size() - 1; ++node )
    {
      // find bracketing nodes - incl shameless hack for evaluations at the boundary
      if ( ( NODES[ node ].real() < x_pos  || std::abs( NODES[ node ].real() - x_pos ) < 1.e-7 ) &&
           ( NODES[ node + 1 ].real() > x_pos || std::abs( NODES[ node + 1 ].real() - x_pos ) < 1.e-7 ) )
      {
        // distance from left node -- real coordinate is given. We also need to
        // interpolate between the two complex nodes 
        std::complex<double> delta_z = ( NODES[ node + 1 ] - NODES[ node ] ) * ( x_pos - NODES[ node ].real() ) / ( NODES[ node + 1 ].real() - NODES[ node ].real() );
        // empty data to return
        Vector<std::complex<double> > left;
        Vector<std::complex<double> > right;
        Vector<std::complex<double> > deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        // derivative of the data
        deriv = ( right - left ) / ( NODES[ node + 1 ] - NODES[ node ] );
        // overwrite right
        right = left + deriv * delta_z;
        return right;
      }
    }
    // If the position is not in the range throw an error
    throw Error( "Mesh error: interpolation not in range " );
  }

  // Output the mesh to a file 
  template < class T, class X >
  void OneD_node_mesh<T, X>::output( std::string filename, int precision ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( precision );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    for ( std::size_t i = 0; i < NODES.size(); ++i )
    {
      dump << NODES[ i ] << " ";
      for ( std::size_t var = 0; var < NV; ++var )
      {
        dump << VARS[ i * NV + var ] << " ";
      }
      dump << "\n";
    }    
  }

  template <>
  void OneD_node_mesh< std::complex<double>, double >::output( std::string filename, int precision ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( precision );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    for ( std::size_t i = 0; i < NODES.size(); ++i )
    {
      dump << NODES[ i ] << " ";
      for ( std::size_t var = 0; var < NV; ++var )
      {
        dump << real( VARS[ i * NV + var ] ) << " " << imag( VARS[ i * NV + var ] ) << " ";
      }
      dump << "\n";
    }
  }


  //the templated versions we require are:
  template class OneD_node_mesh<double>
  ;
  template class OneD_node_mesh< std::complex<double> >
  ;
  template class OneD_node_mesh< std::complex<double>, std::complex<double> >
  ;

} // End of namespace TSL
