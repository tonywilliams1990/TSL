#include <string>
#include <sstream>

#include "Utility.h"

namespace TSL
{
  namespace Utility
  {

    /// Return an integer value as a string
    std::string stringify( const int &val )
    {
      std::stringstream temp;
      temp << val;
      return temp.str();
    }

    /// Return a double value as a string
    std::string stringify( const double &val, int p )
    {
      std::stringstream temp;
      temp.precision( p );
      temp << val;
      return temp.str();
    }

  } // End of namespace Utility
} // End of namespace TSL
