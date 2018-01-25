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
    std::string stringify( const double &val, int p, std::string str )
    {
      std::stringstream temp;
      temp.precision( p );
      if ( str == "fixed" )
      {
        temp << std::fixed << val;
      } else {
        temp << val;
      }
      return temp.str();
    }

    /// Check if a file exists
    bool file_exists( const std::string& name )
    {
      std::ifstream f( name.c_str() );
      return f.good();
    }

  } // End of namespace Utility
} // End of namespace TSL
