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
      row_iter ii;
      row_iter iter;
      // Find the required row iterator
      for(iter=this->CONTAINER.begin(); iter!=this->CONTAINER.end(); iter++)
      {
        if( (*iter).first == row )
        {
          ii = iter;
        }
      }

      col_iter jj;
      std::size_t i(0);
      if ( this->numel_row( row ) > 0 )
      {
        for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
        {
          // for each non-zero elt in the row
          PetscScalar elt;
          elt = std::real( (*jj).second ) + PETSC_i * std::imag( (*jj).second );
          int col( (*jj).first );
          storage[ i ] = elt;
          // +1 to return FORTRAN indexing
          cols[ i ] = col;
          ++i;
        }
      }
    }
#endif

  // the templated versions we require are:
  template class SparseMatrix<double>;
  template class SparseMatrix< std::complex<double> >;

} // end namespace
