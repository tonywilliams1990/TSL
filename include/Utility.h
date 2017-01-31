#include <sstream>
#include <string>

#include "Vector.h"

namespace TSL
{
  namespace Utility
  {

    /// Return the dot product of two vectors
    template <typename T>
    T dot( const Vector<T>& X, const Vector<T>& Y )
    {
      if ( X.size() != Y.size() )	{ throw Error( "Vector dot product: size error" );}
      T dp;
      Vector<T> temp( X );
      dp = temp.dot( Y );
      return dp;
    }

    /// Return an integer value as a string
    std::string stringify( const int &val );

    /// Return a double value as a string
    std::string stringify( const double &val, int p, std::string str="" );

  } // End of namespace Utility
} // End of namespace TSL
