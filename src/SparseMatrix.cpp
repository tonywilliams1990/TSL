/// \file SparseVector.cpp
/// Implementation of the SparseVector class -- a sparse, variable size, vector object.
/// The sparse & dense vectors do not share a common base because the
/// sparse class encapsulates an STL map and operator[] assigns entries to
/// the map. Hence get/set methods are used here.

#include <complex>
#include <map>
#include <cmath>
#include <algorithm>

#include "SparseMatrix.h"
#include "Error.h"
//#include <Functors.h>

namespace TSL
{

#if defined(PETSC_Z)
    template <>
    void SparseMatrix<std::complex<double> >::get_row_petsc( PetscInt row, PetscScalar* storage, PetscInt* cols )
    {
      citer pos;
      std::size_t i(0);

      if ( MATRIX[ row ].nelts() > 0 )
      {
        // start at the begining of this row
        pos = MATRIX[ row ].begin();
        do
        {
          // for each non-zero elt in the row
          PetscScalar elt;
          elt = std::real(pos -> second) + PETSC_i * std::imag(pos -> second);
          int col( pos -> first );
          storage[ i ] = elt;
          // +1 to return FORTRAN indexing
          cols[ i ] = col;
          ++pos;
          ++i;
        }
        while ( pos != MATRIX[ row ].end() );
      }
    }

#endif

#if defined(PETSC_D)
   template <>
   void SparseMatrix<double>::get_row_petsc( PetscInt row, PetscScalar* storage, PetscInt* cols )
     {
     //std::string problem;
     //problem = "The SparseMatrix::get_row_petsc method was called for a SparseMatrix<double>\n";
     //problem += "even though PETSC_ARCH is currently pointing to a complex version of the library.\n";
     //throw ExceptionExternal( problem );
     // iterator to the maps that are used in SparseVector
     // this is bad form as it exposes the internals of the SparseVector storage
     citer pos;
     std::size_t i(0);
     //
     // matrix could be singular with an empty row for the mass matrix
     // of a generalised eigenvalue problem
     if ( MATRIX[row].nelts() > 0 )
     {
       // start at the begining of this row
       pos = MATRIX[ row ].begin();
       do
       {
         // for each non-zero elt in the row
         PetscScalar elt;
         elt = pos -> second;
         int col( pos -> first );
         storage[ i ] = elt;
         // +1 to return FORTRAN indexing
         cols[ i ] = col;
         ++pos;
         ++i;
       }
       while ( pos != MATRIX[ row ].end() );
     }
   }
 #endif

 #if defined(PETSC_D)
   template <>
   void SparseMatrix<std::complex<double> >::get_row_petsc( PetscInt row, PetscScalar* storage, PetscInt* cols )
   {
     std::string problem;
     problem = "The SparseMatrix::get_row_petsc method was called for a SparseMatrix<D_complex>\n";
     problem += "even though PETSC_ARCH is currently pointing to a double version of the library.\n";
     throw Error( problem );
   }
 #endif


  // the templated versions we require are:
  template class SparseMatrix<double>;
  template class SparseMatrix< std::complex<double> >;

} // end namespace
