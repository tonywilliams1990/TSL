/* TwoD_node_mesh - A class specifying a two dimensional mesh object used
                    for storing and manipulating data.
*/

#ifndef TWOD_NODE_MESH_H
#define TWOD_NODE_MESH_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Vector.h"
#include "Error.h"
#include "Matrix.h"
#include "OneD_node_mesh.h"


namespace TSL
{
  

	/// A templated class
	template <class T>
	
	class TwoD_node_mesh
	{

    protected:
      std::size_t NX, NY, NV;                 // Number of nodes and variables
      Vector<double> X_NODES;                 // Vector for storing nodal points (x)
      Vector<double> Y_NODES;                 // Vector for storing nodal points (y)
      Vector<T> VARS;                         // Vector for storing variables at each node

    public:

      /// Constructor (default)
      TwoD_node_mesh()
      {}

      /// Constructor (given nodal distribution)
      TwoD_node_mesh( const Vector<double>& x_nodes, const Vector<double>& y_nodes,
                      const std::size_t& nvars ) :
          NX( x_nodes.size() ),
          NY( y_nodes.size() ),
          NV( nvars ), 
          X_NODES( x_nodes ), 
          Y_NODES( y_nodes )     
      {
        // set the contents to zero
        VARS = Vector<T>( NV * NX * NY, T(0.0) );
      } 

      //TODO constructor from a file
    
      /// Destructor
      virtual ~TwoD_node_mesh() {}

      /* ----- Operator overloading ----- */
    
      /// Access operator for a nodal point/variable in the mesh
      T& operator()( const std::size_t nodex, const std::size_t nodey, 
                     const std::size_t var );

      /// Const access operator for a nodal point/variable in the mesh
      const T& operator()( const std::size_t nodex, const std::size_t nodey, 
                           const std::size_t var ) const;

      /* ----- Methods ----- */

      /// Return the spatial position of a given given node as a pair  
      std::pair<double, double> coord( const std::size_t nodex, 
                                       const std::size_t nodey ) const;

      /// Set the variables stored at a specified node
      void set_nodes_vars( const std::size_t nodex, const std::size_t nodey, 
                           const Vector<T>& U );

      /// Get the variables stored at a specified node
      Vector<T> get_nodes_vars( const std::size_t nodex, const std::size_t nodey ) const;
           
      /// Get a cross section of the 2D mesh at a specified (constant) x node
      OneD_node_mesh<T> get_xsection_at_xnode( const std::size_t nodex ) const;

      /// Get a cross section of the 2D mesh at a specified (constant) y node
      OneD_node_mesh<T> get_xsection_at_ynode( const std::size_t nodey ) const; 

      /// Assign an element to all entries in the mesh
      void assign( const T elt );

      /// Get the number of nodes in the two directions of the 2D mesh
      std::pair< std::size_t, std::size_t > get_nnodes() const;

      /// Get the number of variables that are stored at each node
      std::size_t get_nvars() const;

      /// Return A vector of the x-nodal positions for this mesh
      const Vector<double>& xnodes() const;

      /// Return A vector of the y-nodal positions for this mesh
      const Vector<double>& ynodes() const;

      /// Return a matrix corresponding to each nodal point in the mesh
      /// Each matrix element will contain a specified variable number
      Matrix<T> get_var_as_matrix( std::size_t var ) const;

      /// Interpolate this mesh data (bilinearly) into a new
      /// mesh with nodal points defined in the argument list.
      void remesh1( const Vector<double>& newX, const Vector<double>& newY );

      /// A simple method for dumping data to std::cout
      void dump() const;

      /// A simple method for dumping data to a file
      void dump( std::string filename ) const;
      
      /// A simple method for dumping a single variable to a file with no nodal information
      void dump_var( std::string filename, const std::size_t var ) const;

      /// A simple method for reading data from a file 
      void read( std::string filename, const bool reset = false );

      /// A simple method for dumping data to a file for gnuplot surface plotting
      void dump_gnu( std::string filename ) const;

      /// Get a bilinearly interpolated value at a specified point
      Vector<T> get_interpolated_vars( const double& x, const double& y)
      {
        // check start and end
        if ( ( x < X_NODES[0] ) || ( x>X_NODES[NX-1] ) )
        {
          std::string problem;
          problem = "The TwoD_node_mesh.get_interpolated_vars method has been called \n";
          problem += "with an x coordinate that lies outside the mesh. \n";
          throw Error( problem );
        }
        // check start and end
        if ( ( y < Y_NODES[0] ) || ( y>Y_NODES[NY-1] ) )
        {
          std::string problem;
          problem = "The TwoD_node_mesh.get_interpolated_vars method has been called \n";
          problem += "with a y coordinate that lies outside the mesh. \n";
          throw Error( problem );
        }
        int bottom_j(-2);
        for ( unsigned j = 0; j < NY-1; ++j )
	      {
	        if ( ( y > Y_NODES[j] ) && ( y < Y_NODES[j+1] ) )
	        {
	          bottom_j = j;
	        }
	        if ( ( abs(y-Y_NODES[j])<1.e-10 ) || ( abs(y-Y_NODES[j+1])<1.e-10 ) )
	        {
	          bottom_j = j;
	        }
	      }
        if ( bottom_j == -1 )
	      {
          std::string problem;
          problem = " The TwoD_node_mesh.get_interpolated_vars method is broken.\n";
          throw Error( problem );
	      }

        std::cout << y << " " << Y_NODES[bottom_j] << " " << Y_NODES[bottom_j+1] << "\n";
        //
        OneD_node_mesh<T> bottom_row = get_xsection_at_ynode( bottom_j );
        OneD_node_mesh<T> top_row = get_xsection_at_ynode( bottom_j+1 );
        const double y1 = Y_NODES[ bottom_j ]; 
        const double y2 = Y_NODES[ bottom_j+1 ];
        Vector<T> result = top_row.get_interpolated_vars(x)*( y2-y )/( y2-y1 )
          + bottom_row.get_interpolated_vars(x)*( y-y1 )/( y2-y1 );
        std::cout << "x,y,interp: " << x << " " << y << " " << result[0] << "\n"; 
        return result; 
      }





	}; // End of class TwoD_node_mesh

  template <class T>
  inline T& TwoD_node_mesh<T>::operator()( const std::size_t nodex, 
                                           const std::size_t nodey, const std::size_t var )
  {
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw Error( problem );
    }

    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <class T>
  inline const T& TwoD_node_mesh<T>::operator()( const std::size_t nodex, 
                                  const std::size_t nodey, const std::size_t var ) const
  {
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw Error( problem );
    }
    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <class T>
  inline std::pair<double, double> TwoD_node_mesh<T>::coord( const std::size_t nodex, 
                                                        const std::size_t nodey ) const
  {
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.coord method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    std::pair< double, double > pos;
    pos.first = X_NODES[ nodex ];
    pos.second = Y_NODES[ nodey ];
    return pos;
  }

} // End of namespace TSL

#endif
