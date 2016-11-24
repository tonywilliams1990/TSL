
#include "Vector.h"

namespace TSL
{
  namespace Utility
  {

    template <typename T>
    T dot( const Vector<T>& X, const Vector<T>& Y )
    {
      if ( X.size() != Y.size() )	{ throw Error( "Vector dot product: size error" );}
      T dp;
      Vector<T> temp( X );
      dp = temp.dot( Y );
      return dp;
    }

  } // End of namespace Utility
} // End of namespace TSL
