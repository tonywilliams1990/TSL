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


  // the templated versions we require are:
  template class SparseMatrix<double>;
  template class SparseMatrix< std::complex<double> >;

} // end namespace
