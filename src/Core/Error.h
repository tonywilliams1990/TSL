/* Error.h - Defines a class for simple error handling.
   Created - 06/08/2016
*/

#ifndef ERROR_H
#define ERROR_H

#include <string>
#include <iostream>
#include <stdexcept>

namespace TSL
{

  /// A generic runtime exception
  class Error : public std::runtime_error
  {
  public:
    Error( const std::string &problem ) : std::runtime_error( problem )
    {
      error_header();
      std::cout << problem << std::endl;
    }

  private:

    void error_header()
    {
      std::cout << "----------------------------------------------------" << std::endl;
      std::cout << " Error: A TSL routine has had a RUNTIME problem! " << std::endl;
      std::cout << "----------------------------------------------------" << std::endl;
    }

  };	// End of class Error

}  // End of namespace TSL

#endif
