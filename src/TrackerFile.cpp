/* TrackerFile class implementaion 
*/

#include "TrackerFile.h"

namespace TSL
{

  /// Constructor ( no filename )
  TrackerFile::TrackerFile( std::size_t prec )
  {
    this -> PREC = prec;
  }

  /// Constructor ( with filename )
  TrackerFile::TrackerFile( std::string filename, std::size_t prec )
  {
    file.open( filename.c_str() );
    this -> PREC = prec;
    file.precision( prec );
    file << std::scientific;
  }

  /// Destructor
  TrackerFile::~TrackerFile()
  {
    file.close();
  }

  /* ----- Methods ----- */

  /// Set the output precision
  void TrackerFile::precision( std::size_t prec )
  {
    this -> PREC = prec;
    file.precision( prec );
  }

  /// Set the filename
  void TrackerFile::set_filename( std::string filename )
  {
    file.close();
    file.open( filename.c_str() );
    file.precision( PREC );
  }

  /// Push pointer to a double variable for tracking
  void TrackerFile::push_ptr( double* scalar, std::string desc )
  {
    ptr_DOUBLES.push_back( scalar );
    DOUBLE_DESC.push_back( desc );
  }

  /// Push pointer to a std::complex<double> for tracking
  void TrackerFile::push_ptr( std::complex<double>* scalar, std::string desc )
  {
    double * r = &reinterpret_cast<double(&)[2]>( scalar[0] )[0];
    double * i = &reinterpret_cast<double(&)[2]>( scalar[0] )[1];
    ptr_DOUBLES.push_back( r );
    ptr_DOUBLES.push_back( i );
    // Old implementation worked with some compilers but not correct C++
    //ptr_DOUBLES.push_back( &( scalar -> real() ) );
    //ptr_DOUBLES.push_back( &( scalar -> imag() ) );
    DOUBLE_DESC.push_back( desc + " (real)" );
    DOUBLE_DESC.push_back( desc + " (imag)" );
  }

  /// Push pointer to Vector<double> for tracking
  void TrackerFile::push_ptr( Vector<double>* ptr_to_vector, std::string desc )
  {
    ptr_DVECTORS.push_back( ptr_to_vector );
    DVECTOR_DESC.push_back( desc );
  }

  /// Push pointer to a Vector<std::complex<double>> for tracking
  void TrackerFile::push_ptr( Vector< std::complex<double> >* ptr_to_vector,
                              std::string desc )
  {
    ptr_CVECTORS.push_back( ptr_to_vector );
    CVECTOR_DESC.push_back( desc );
  }

  /// Push pointer to a 1D mesh of doubles for tracking
  void TrackerFile::push_ptr( OneD_node_mesh<double>* ptr_to_mesh, std::string desc )
  {
    ptr_DMESH.push_back( ptr_to_mesh );
    DMESH_DESC.push_back( desc );
  }

  /// Push pointer to a 1D mesh of std::complex<double> for tracking
  void TrackerFile::push_ptr( OneD_node_mesh< std::complex<double> >* ptr_to_mesh,
                              std::string desc )
  {
    ptr_CMESH.push_back( ptr_to_mesh );
    CMESH_DESC.push_back( desc );
  }

  /// Push pointer to a complex 1D mesh of std::complex<double> for tracking
  void TrackerFile::push_ptr( OneD_node_mesh< std::complex<double>, 
                              std::complex<double> >* ptr_to_mesh, std::string desc )
  {
    ptr_CCMESH.push_back( ptr_to_mesh );
    CCMESH_DESC.push_back( desc );
  }

  /// Create a header in the file
  void TrackerFile::header()
  {
    // write the header
    file << " # Header : \n # ";
    for ( std::size_t i = 0; i < DOUBLE_DESC.size(); ++i )
    {
      file << DOUBLE_DESC[ i ] << " | ";
    }
    for ( std::size_t i = 0; i < DVECTOR_DESC.size(); ++i )
    {
      file << DVECTOR_DESC[ i ] << " | ";
    }
    for ( std::size_t i = 0; i < CVECTOR_DESC.size(); ++i )
    {
      file << CVECTOR_DESC[ i ] + " (Real)" << " | ";
      file << CVECTOR_DESC[ i ] + " (Imag)" << " | ";
    }
    
    for ( std::size_t i = 0; i < DMESH_DESC.size(); ++i )
    {
      file << DMESH_DESC[ i ] + " (nodes)" << " | ";
      for ( std::size_t var = 0; var < ( *ptr_DMESH[ i ] ).get_nvars(); ++var )
      {
        file << DMESH_DESC[ i ] + " var#" + make_string( var ) + " | ";
      }
    }
    
    for ( std::size_t i = 0; i < CMESH_DESC.size(); ++i )
    {
      file << CMESH_DESC[ i ] + "(nodes)" << " | ";
      for ( std::size_t var = 0; var < ( *ptr_CMESH[ i ] ).get_nvars(); ++var )
      {
        file << CMESH_DESC[ i ] + " var#" + make_string( var ) + " (Real) | ";
        file << CMESH_DESC[ i ] + " var#" + make_string( var ) + " (Imag) | ";
      }
    }
    
    for ( std::size_t i = 0; i < CCMESH_DESC.size(); ++i )
    {
      file << CCMESH_DESC[ i ] + "(nodes)_real" << " | ";
      file << CCMESH_DESC[ i ] + "(nodes)_imag" << " | ";
      for ( std::size_t var = 0; var < ( *ptr_CCMESH[ i ] ).get_nvars(); ++var )
      {
        file << CCMESH_DESC[ i ] + " var#" + make_string( var ) + " (Real) | ";
        file << CCMESH_DESC[ i ] + " var#" + make_string( var ) + " (Imag) | ";
      }
    } 
    file << std::endl;
  }

  /// Private method for writing scalar data to file
  void TrackerFile::write_scalar_data()
  {
    if ( !ptr_DOUBLES.empty() )
    {
      // simple flat data file
      for ( std::size_t i = 0; i < ptr_DOUBLES.size(); ++i )
      {
        file << *ptr_DOUBLES[ i ] << " ";
      }
    }
  }

  void TrackerFile::update()
  {
    std::size_t block_size ( 1 );
    if ( !ptr_DVECTORS.empty() )
    {
      block_size = ptr_DVECTORS[ 0 ] -> size();
    }
    if ( !ptr_CVECTORS.empty() )
    {
      block_size = ptr_CVECTORS[ 0 ] -> size();
    }
    if ( !ptr_DMESH.empty() )
    {
      block_size = ptr_DMESH[ 0 ] -> get_nnodes();
    }
    if ( !ptr_CMESH.empty() )
    {
      block_size = ptr_CMESH[ 0 ] -> get_nnodes();
    }
    if ( !ptr_CCMESH.empty() )
    {
      block_size = ptr_CCMESH[ 0 ] -> get_nnodes();
    }

    for ( std::size_t line = 0; line < block_size; ++line )
    {
      write_scalar_data();
      if ( !ptr_DVECTORS.empty() )
      {
        // for each vector ptr
        for ( std::size_t i = 0; i < ptr_DVECTORS.size(); ++i )
        {
          file << ( *ptr_DVECTORS[ i ] ) [ line ] << " ";
        }
      }
      if ( !ptr_CVECTORS.empty() )
      {
        // for each vector ptr
        for ( std::size_t i = 0; i < ptr_CVECTORS.size(); ++i )
        {
          file << ( *ptr_CVECTORS[ i ] ) [ line ].real() << " ";
          file << ( *ptr_CVECTORS[ i ] ) [ line ].imag() << " ";
        }
      }
      if ( !ptr_DMESH.empty() )
      {
        // for each mesh ptr
        for ( std::size_t i = 0; i < ptr_DMESH.size(); ++i )
        {
          file << ( *ptr_DMESH[ i ] ).coord( line ) << " ";
          for ( std::size_t var = 0; var < ptr_DMESH[ i ] -> get_nvars(); ++var )
          {
            file << ( *ptr_DMESH[ i ] )( line, var ) << " ";
          }
        }
      }
      if ( !ptr_CMESH.empty() )
      {
        // for each mesh ptr
        for ( std::size_t i = 0; i < ptr_CMESH.size(); ++i )
        {
          file << ( *ptr_CMESH[ i ] ).coord( line ) << " ";
          for ( std::size_t var = 0; var < ptr_CMESH[ i ] -> get_nvars(); ++var )
          {
            file << ( *ptr_CMESH[ i ] )( line, var ).real() << " ";
            file << ( *ptr_CMESH[ i ] )( line, var ).imag() << " ";
          }
        }
      }
      if ( !ptr_CCMESH.empty() )
      {
        // for each mesh ptr
        for ( std::size_t i = 0; i < ptr_CCMESH.size(); ++i )
        {
          file << ( *ptr_CCMESH[ i ] ).coord( line ).real() << " ";
          file << ( *ptr_CCMESH[ i ] ).coord( line ).imag() << " ";
          for ( std::size_t var = 0; var < ptr_CCMESH[ i ] -> get_nvars(); ++var )
          {
            file << ( *ptr_CCMESH[ i ] )( line, var ).real() << " ";
            file << ( *ptr_CCMESH[ i ] )( line, var ).imag() << " ";
          }
        }
      }
      file << std::endl;
    }
    // flush the buffer
    file.flush();
  }

} // End of namespace TSL
