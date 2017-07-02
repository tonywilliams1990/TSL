/* TrackerFile class - A class that can be passed pointers to scalars, vectors and mesh 
                       objects so that their contents can be written to a file each time the
                       update method is called. 
*/

#ifndef TRACKERFILE_H
#define TRACKERFILE_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "Error.h"
#include "Vector.h"
#include "OneD_node_mesh.h"

namespace TSL
{

  class TrackerFile
  {
    private:
      std::size_t PREC;                              // Output precision
      std::vector< double* > ptr_DOUBLES;            // Doubles
      std::vector< Vector<double>* > ptr_DVECTORS;   // Vector<double>
      
      // Vector<complex>
      std::vector< Vector< std::complex<double> >* > ptr_CVECTORS; 
      // 1D mesh<double>
      std::vector< OneD_node_mesh<double>* > ptr_DMESH;            
      // 1D mesh<complex>
      std::vector< OneD_node_mesh< std::complex<double> >* > ptr_CMESH;         
      /// A vector of pointers to complex-complex 1D meshes
      std::vector< OneD_node_mesh< std::complex<double>, 
                   std::complex<double> >* > ptr_CCMESH;

      /// a vector of pointers to descriptions
      std::vector< std::string > DOUBLE_DESC;        // Doubles
      std::vector< std::string > DVECTOR_DESC;       // Vector of doubles
      std::vector< std::string > CVECTOR_DESC;       // Vector of complex numbers
      std::vector< std::string > DMESH_DESC;         // Mesh of doubles
      std::vector< std::string > CMESH_DESC;         // Mesh of complex numbers
      std::vector< std::string > CCMESH_DESC; 

      /// Make the value into a string
      std::string make_string( const std::size_t& val )
      {
        std::stringstream temp;
        temp << val;
        return temp.str();
      }

      /// Private method for writing scalar data to file
      void write_scalar_data();

    protected:

      /// Output stream
      mutable std::ofstream file;

    public:

      /// Constructor ( no filename )
      TrackerFile( std::size_t prec = 4 );

      /// Constructor ( with filename )
      TrackerFile( std::string filename, std::size_t prec = 12 );

      /// Destructor
      ~TrackerFile();

      /* ----- Methods ----- */

      /// Set the output precision
      void precision( std::size_t prec );

      /// Set the filename 
      void set_filename( std::string filename );

      /// Push pointer to a double variable for tracking
      void push_ptr( double* scalar, std::string desc = "" );

      /// Push pointer to a std::complex<double> for tracking
      void push_ptr( std::complex<double>* scalar, std::string desc = "" );

      /// Push pointer to a Vector<double> for tracking
      void push_ptr( Vector<double>* ptr_to_vector, std::string desc = "" );

      /// Push pointer to a Vector<std::complex<double>> for tracking
      void push_ptr( Vector< std::complex<double> >* ptr_to_vector, 
                     std::string desc = "" );

      /// Push pointer to a 1D mesh of doubles for tracking
      void push_ptr( OneD_node_mesh<double>* ptr_to_mesh, std::string desc = "" );

      /// Push pointer to a 1D mesh of std::complex<double> for tracking
      void push_ptr( OneD_node_mesh< std::complex<double> >* ptr_to_mesh, 
                     std::string desc = "" );

      /// Push pointer to a complex 1D mesh of std::complex<double> for tracking
      void push_ptr( OneD_node_mesh< std::complex<double>, 
                     std::complex<double> >* ptr_to_mesh, std::string desc = "" );

      /// Create a new line in the file
      void newline(){ file << std::endl; }

      /// Create a header in the file
      void header();

      /// Update the file
      void update();


  }; // End of class TrackerFile

} // End of namespace TSL

#endif
