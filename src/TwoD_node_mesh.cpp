/* TwoD_node_mesh - A class specifying a two dimensional mesh object used
                    for storing and manipulating data.
*/

#include "TwoD_node_mesh.h"

namespace TSL
{
  template <class T>
  void TwoD_node_mesh<T>::set_nodes_vars( const std::size_t nodex, const std::size_t nodey, 
                                          const Vector<T>& U )
  {
    if ( U.size() > NV )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.set_nodes_vars method is trying to use a \n";
      problem += " vector that has more entries than variables stored in the mesh. \n";
      throw Error( problem );
    }
    // assign contents of U to the member data
    std::size_t offset( ( nodex * NY + nodey ) * NV );
    for ( std::size_t var = 0; var < NV; ++var )
    {
      VARS[ offset++ ] = U[ var ];
    }
  }
  
	template <class T>
  Vector<T> TwoD_node_mesh<T>::get_nodes_vars( const std::size_t nodex, 
                                               const std::size_t nodey ) const
  {
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.get_nodes_vars method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    //Construct a vector with NV elements 
    Vector<T> nodes_vars;
    for ( std::size_t var = 0; var < NV; ++var )
    {
      nodes_vars.push_back( VARS[ ( nodex * NY + nodey ) * NV + var ] );
    }
    //TODO check this is working correctly
    return nodes_vars;
  }

  template <class T>
  void TwoD_node_mesh<T>::assign( const T elt )
  {
    VARS.assign( NX * NY * NV, elt );
  }

  template <class T>
  std::pair< std::size_t, std::size_t > TwoD_node_mesh<T>::get_nnodes() const
  {
    std::pair< std::size_t, std::size_t > nodes;
    nodes.first = NX;
    nodes.second = NY;
    return nodes;
  }

  template <class T>
  std::size_t TwoD_node_mesh<T>::get_nvars() const
  {
    return NV;
  }

  template <class T>
  const Vector<double>& TwoD_node_mesh<T>::xnodes() const
  {
    return X_NODES;
  }

  template <class T>
  const Vector<double>& TwoD_node_mesh<T>::ynodes() const
  {
    return Y_NODES;
  }

  template <class T>
  Matrix<T> TwoD_node_mesh<T>::get_var_as_matrix( std::size_t var ) const
  {
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.get_var_as_matrix method is trying to use a \n";
      problem += " variable index bigger than the number of variables in the mesh. \n";
      throw Error( problem );
    }
    Matrix<T> temp( NX, NY, 0.0 );
    for ( std::size_t i = 0; i < NX; ++i )
    {
      for  ( std::size_t j = 0; j < NY; ++j )
      {
        temp( i, j ) = VARS[ ( i * NY + j ) * NV + var ];
      }
    }
    return temp;
  }  
  
  template<class T>
  void TwoD_node_mesh<T>::remesh1( const Vector<double>& newX, const Vector<double>& newY )
  {
    // check start & end 
    if ( std::abs( X_NODES[ 0 ] - newX[ 0 ] ) > 1.e-10 ||
         std::abs( X_NODES[ X_NODES.size() - 1 ] - newX[ newX.size() - 1 ] ) > 1.e-10 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.remesh1 method has been called with \n";
      problem += " a passed X coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw Error( problem );
    }
    // check monotonic node positions
    for ( std::size_t i = 0; i < newX.size() - 1; ++i )
    {
      if ( newX[ i ] >= newX[ i + 1 ] )
      {
        std::string problem;
        problem = " The TwoD_node_mesh.remesh1 method has been passed \n";
        problem += " a non-monotonic X coordinate vector. \n";
        throw Error( problem );
      }
    }
    // check start and end
    if ( std::abs( Y_NODES[ 0 ] - newY[ 0 ] ) > 1.e-10 ||
         std::abs( Y_NODES[ Y_NODES.size() - 1 ] - newY[ newY.size() - 1 ] ) > 1.e-10 )
    {
      std::string problem;
      problem = " The TwoD_node_mesh.remesh1 method has been called with \n";
      problem += " a passed Y coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw Error( problem );
    }
    // check monotonic node positions
    for ( std::size_t i = 0; i < newY.size() - 1; ++i )
    {
      if ( newY[ i ] >= newY[ i + 1 ] )
      {
        std::string problem;
        problem = " The TwoD_node_mesh.remesh1 method has been passed \n";
        problem += " a non-monotonic Y coordinate vector. \n";
        throw Error( problem );
      }
    }

    // new variables storage
    Vector<T> newvars( newX.size() * newY.size() * NV, 0.0 );

    // left boundary
    {
      std::size_t xnode( 0 );
      // bottom left corner copy
      for ( unsigned var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + 0 ) * NV + var ] = get_nodes_vars( 0, 0 )[ var ];
      }
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaY( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
        {
          if ( ( Y_NODES[ j ] <= newY[ ynode ] ) && ( newY[ ynode ] < Y_NODES[ j + 1 ] ) )
          {
            below_j = j;
            deltaY = newY[ ynode ] - Y_NODES[ j ];
          }
        }
        Vector<T> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 ) 
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i, below_j + 1 ).second 
                              - coord( left_i, below_j ).second );
        Vector<T> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdY * deltaY;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
      // top left corner copy
      for ( std::size_t var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + newY.size() - 1 ) * NV + var ] = 
        get_nodes_vars( 0, NY - 1 )[ var ];
      }
    }
    // right boundary
    {
      std::size_t xnode( newX.size() - 1 );
      // bottom right corner copy
      for ( std::size_t var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + 0 ) * NV + var ] = 
        get_nodes_vars( NX - 1, 0 )[ var ];
      }
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( X_NODES.size() - 1 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaY( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
        {
          if ( ( Y_NODES[ j ] <= newY[ ynode ] ) && ( newY[ ynode ] < Y_NODES[ j + 1 ] ) )
          {
            below_j = j;
            deltaY = newY[ ynode ] - Y_NODES[ j ];
          }
        }
        Vector<T> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 ) 
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i, below_j + 1 ).second 
                              - coord( left_i, below_j ).second );
        Vector<T> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdY * deltaY;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
      // bottom right corner copy
      for ( std::size_t var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + newY.size() - 1 ) * NV + var ] = 
        get_nodes_vars( NX - 1, NY - 1 )[ var ];
      }
    }
    // bottom boundary
    {
      std::size_t ynode( 0 );
      for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaX( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X_NODES.size() - 1; ++i )
        {
          if ( ( X_NODES[ i ] <= newX[ xnode ] ) && ( newX[ xnode ] < X_NODES[ i + 1 ] ) )
          {
            left_i = i;
            deltaX = newX[ xnode ] - X_NODES[ i ];
          }
        }
        Vector<T> dvarsdX = ( get_nodes_vars( left_i + 1, below_j ) 
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i + 1, below_j ).first 
                              - coord( left_i, below_j ).first );
        Vector<T> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdX * deltaX;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
    }
    // top boundary
    {
      std::size_t ynode( newY.size() - 1 );
      for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( Y_NODES.size() - 1 ); // bracketing index
        double deltaX( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X_NODES.size() - 1; ++i )
        {
          if ( ( X_NODES[ i ] <= newX[ xnode ] ) && ( newX[ xnode ] < X_NODES[ i + 1 ] ) )
          {
            left_i = i;
            deltaX = newX[ xnode ] - X_NODES[ i ];
          }
        }
        Vector<T> dvarsdX = ( get_nodes_vars( left_i + 1, below_j ) 
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i + 1, below_j ).first 
                              - coord( left_i, below_j ).first );
        Vector<T> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdX * deltaX;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
    }
    // loop thru interior nodes of the destination mesh one node at a time
    for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
    {
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X_NODES.size() - 1; ++i )
        {
          if ( ( X_NODES[ i ] <= newX[ xnode ] ) && ( newX[ xnode ] < X_NODES[ i + 1 ] ) )
          {
            left_i = i;
          }
        }
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
        {
          if ( ( Y_NODES[ j ] <= newY[ ynode ] ) && ( newY[ ynode ] < Y_NODES[ j + 1 ] ) )
          {
            below_j = j;
          }
        }
        Vector<T> dvarsdX = ( get_nodes_vars( left_i + 1, below_j ) 
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i + 1, below_j ).first 
                              - coord( left_i, below_j ).first );
        Vector<T> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 ) 
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i, below_j + 1 ).second 
                              - coord( left_i, below_j ).second );

        Vector<T> interpolated_vars_bottom =
          ( get_nodes_vars( left_i, below_j ) * ( coord( left_i + 1, below_j ).first 
            - newX[ xnode ] ) + get_nodes_vars( left_i + 1, below_j ) * ( newX[ xnode ] -
           coord( left_i, below_j ).first ) ) / ( coord( left_i + 1, below_j ).first 
          - coord( left_i, below_j ).first );

        Vector<T> interpolated_vars_top =
          ( get_nodes_vars( left_i, below_j + 1 ) * ( coord( left_i + 1, below_j 
          + 1 ).first - newX[ xnode ] ) + get_nodes_vars( left_i + 1, below_j + 1 ) 
          * ( newX[ xnode ] - coord( left_i, below_j + 1 ).first ) ) /
          ( coord( left_i + 1, below_j + 1 ).first - coord( left_i, below_j + 1 ).first );

        Vector<T> interpolated_vars =
          (  interpolated_vars_bottom * ( coord( left_i, below_j + 1 ).second 
          - newY[ ynode ] ) +  interpolated_vars_top 
          * ( newY[ ynode ] - coord( left_i, below_j ).second ) ) /
          ( coord( left_i, below_j + 1 ).second - coord( left_i, below_j ).second );

        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
    }
    // finally replace the old nodes with the new ones
    X_NODES = newX;
    Y_NODES = newY;
    NX = newX.size();
    NY = newY.size();
    VARS = newvars;
  }

  template<class T>
  OneD_node_mesh<T> TwoD_node_mesh<T>::get_xsection_at_xnode( 
                                                          const std::size_t nodex ) const
  {
    OneD_node_mesh<T> xsection( Y_NODES, NV );
    for ( std::size_t nodey = 0; nodey < NY; ++nodey )
    {
      xsection.set_nodes_vars( nodey, this -> get_nodes_vars( nodex, nodey ) );
    }
    return xsection;
  }

  template<class T>
  OneD_node_mesh<T> TwoD_node_mesh<T>::get_xsection_at_ynode( 
                                                          const std::size_t nodey ) const
  {
    OneD_node_mesh<T> xsection( X_NODES, NV );
    for ( std::size_t nodex = 0; nodex < NX; ++nodex )
    {
      xsection.set_nodes_vars( nodex, this -> get_nodes_vars( nodex, nodey ) );
    }
    return xsection;
  }

  template <class T>
  void TwoD_node_mesh<T>::dump() const
  {
    for ( std::size_t var = 0; var < NV; ++var )
    {
      std::cout << "Variable : " << var << "\n";
      std::cout << " x = ";
      for ( std::size_t i = 0; i < NX; ++i )
      {
        std::cout << X_NODES[ i ] << ", ";
      }
      std::cout << "\n";
      for ( std::size_t j = 0; j < NY; ++j )
      {
        std::cout << " y = " << Y_NODES[ j ] << "\n";
        for ( std::size_t i = 0; i < NX; ++i )
        {
          std::cout << VARS[ ( i * NY + j ) * NV + var ] << ", ";
        }
        std::cout << "\n";
      }
    }
  }

  template<>
  void TwoD_node_mesh<double>::dump_gnu( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );

    for ( std::size_t i = 0; i < NX; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
      {
        dump << X_NODES[ i ] << " " << Y_NODES[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << VARS[ ( i * NY + j ) * NV + var ] << " ";
        }
        dump << "\n";
      }
      dump << "\n";
    }
    dump.close();
  }

  //TODO complex version of above 419 TwoD_Node_Mesh.cpp

  template<class T>
  void TwoD_node_mesh<T>::dump_var( std::string filename, const std::size_t var ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    dump.precision( 9 );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        dump << VARS[ ( i * NY + j ) * NV + var ] << "\n";
      }
    }
  }

  template<class T>
  void TwoD_node_mesh<T>::dump( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    //dump << NX << " " << NY << " " << NV << "\n";
    dump.precision( 9 );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        dump << X_NODES[ i ] << " " << Y_NODES[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << VARS[ ( i * NY + j ) * NV + var ] << " ";
        }
        dump << "\n";
      }
    }
  }

  template<class T>
  void TwoD_node_mesh<T>::read( std::string filename, bool reset )
  {
    std::ifstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        double x, y;
        dump >> x;
        dump >> y;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          double value;
          dump >> value;
          VARS[ ( i * NY + j ) * NV + var ] = value;
        }
        if ( reset != true )
        {
          // if not reseting the mesh we should check the node positions
          if ( ( std::abs( x - X_NODES[ i ] ) > 1.e-6 ) || 
               ( std::abs( y - Y_NODES[ j ] ) > 1.e-6 ) )
          {
            std::cout << " Read x = " << x << " Expected x = " << X_NODES[ i ] 
                      << "; Read y = " << y << " Expected y = " << Y_NODES[ j ] << " \n";
            std::cout << " Absolute differences are " << abs( x - X_NODES[i] ) << " and " 
                      << abs( y - Y_NODES[j] ) << "\n";              
            std::string problem;
            problem = " The TwoD_node_mesh.read method is trying to read a \n";
            problem += " file whose nodal points are in a different position. \n";
            throw Error( problem );
          }
        }
        else
        {
          X_NODES[ i ] = x;
          Y_NODES[ j ] = y;
        }
      }
    }
  }


  //the templated versions we require are:
  template class TwoD_node_mesh<double>
  ;
  template class TwoD_node_mesh< std::complex<double> >
  ;
  

} // End of namespace TSL
