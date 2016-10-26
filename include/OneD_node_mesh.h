/* OneD_node_mesh - A class specifying a one dimensional mesh object used
                    for storing and manipulating data.
*/

#ifndef ONED_NODE_MESH_H
#define ONED_NODE_MESH_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Vector.h"
#include "Error.h"


namespace TSL
{
  

	/// A templated class
	template <class T, class X = double>
	
	class OneD_node_mesh
	{

    protected:
      std::size_t NV;                         // Number of variables
      Vector<X> NODES;                        // Vector for storing nodal points
      Vector<T> VARS;                         // Vector for storing variables at each node

    public:

      /// Constructor (default)
      OneD_node_mesh()
      {
        NV = 0;
      }

      /// Constructor (given nodal distribution)
      OneD_node_mesh( const Vector<X>& nodes, const std::size_t& nvars ) : 
          NV( nvars ), NODES( nodes )
      {
        // set the contents to zero
        VARS = Vector<T>( NV * NODES.size(), T(0.0) );
      } 
    
      /// Destructor
      virtual ~OneD_node_mesh() {}

      /* ----- Operator overloading ----- */
    
      /// Read only access of a variable at a given node
      const T& operator()( const std::size_t i, const std::size_t var ) const
      {
        return VARS[ i * NV + var ];
      }

      /// Read and write access of a variable at a given node
      T& operator()( const std::size_t i, const std::size_t var )
      {
        return VARS[ i * NV + var ];
      }

      /* ----- Methods ----- */
      
      /// Read only access of the nodal value
      const X& coord( const std::size_t& node ) const
      {
        return NODES[ node ];
      }

      /// Read and write access of the nodal value
      X& coord( const std::size_t& node )
      {
        return NODES[ node ];
      }

      /// Set the variables stored at a specific node
      void set_nodes_vars( const std::size_t& node, const Vector<T>& u )
      {
        if ( u.size() != NV ) { throw Error( "Mesh error: set_nodes_vars " );}
        for ( std::size_t var = 0; var < u.size(); ++var )
        {
          VARS[ node * NV + var ] = u[ var ];
        }
      }
      
      /// Get the variables stored at a specific node
		  Vector<T> get_nodes_vars( const std::size_t& node ) const
      {
        if ( ( node >= NODES.size() ) || ( node < 0 ) )
        {
          throw Error( "Mesh error: get_nodes_vars range error." );
        }
        Vector<T> nodes_vars;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          nodes_vars.push_back( VARS[ node * NV + var ] );
        }
        return nodes_vars;
      }

      /// Get the variable data at an interpolated position using a first order scheme.
      Vector<T> get_interpolated_vars( const X& pos ) const;
      
      /// Return the number of nodal points in the mesh
      std::size_t get_nnodes() const { return NODES.size(); }

      /// Return the number of variables stored in the mesh
      std::size_t get_nvars() const { return NV; }

      /// Return the vector of nodal positions
      const Vector<X>& nodes() const { return NODES; }

      /// Output data to a file
      void output( std::string filename, int precision = 10 ) const;

      /// Return a vector of the variables stored in the mesh
      const Vector<T>& vars_as_vector() const { return VARS; }

      /// Set the variables of this mesh from a vector
      void set_vars_from_vector( const Vector<T>& vec )
      {
        if ( vec.size() != NV * NODES.size() ) 
        { throw Error( "Mesh error: set_vars_from_vector sizes do not agree." ); }
        VARS = vec;
      }

	}; // End of class OneD_node_mesh


} // End of namespace TSL

#endif
