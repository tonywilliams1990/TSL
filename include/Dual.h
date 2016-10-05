/* Dual - A class of generalised dual numbers 
*/

#ifndef DUAL_H
#define DUAL_H

#include <iostream>

#include "Vector.h"
#include "Error.h"

namespace TSL
{

  template<class T>
	
	class Dual 
	{
		
      protected:
        T REAL;                          // Real part of the dual number
        Vector<T> EPSILONS;              // Vector of epsilons of the dual number

      public:
      
        /// Default constructor
        Dual() : REAL( 0.0 ) { }

        /// Initialised constructor (Generalised dual number)
        Dual( const T& real, const Vector<T>& eps ) : REAL( real ), 
                                                                EPSILONS( eps ) { }

        /// Destructor
        ~Dual() { }

        /* ----- Operator overloading ----- */

        /// Output operator <<
        template <class Type>
			  friend std::ostream& operator<<( std::ostream& os, const Dual<Type>& dual );

        /// Assignment
        Dual<T>& operator=( const Dual<T>& original )
        {
          REAL = original.REAL;
          EPSILONS.resize( original.EPSILONS.size() );
          EPSILONS = original.EPSILONS;
          return *this;
        }

        /// Indexing operator ( read only )
			  const T& operator[] ( const std::size_t& i ) const
			  {
					  // Range check 
					  if ( i<0 || EPSILONS.size()<i )	{ throw Error( "Dual range error" );}
					  if ( i==0 ) { return REAL; }
            else { return EPSILONS[ i - 1 ]; }
			  }

        /// Indexing operator ( read/write )
			  T& operator[] ( const std::size_t& i )
			  {
					  // Range check 
					  if ( i<0 || EPSILONS.size()<i )	{ throw Error( "Dual range error" );}
					  if ( i==0 ) { return REAL; }
            else { return EPSILONS[ i - 1 ]; }
			  }

        /// Unary +
        Dual<T> operator+() const
			  {
				  return *this;
			  }  
      
        /// Unary -
        Dual<T> operator-() const
        {
          Dual<T> temp( *this );
          temp.REAL = -REAL;
		      temp.EPSILONS = -EPSILONS;
		      return temp;
        }

        /// Binary +
        Dual<T> operator+( const Dual<T>& d_plus ) const
			  {
				  if ( d_plus.EPSILONS.size() != EPSILONS.size() )
          { throw Error( "Dual epsilon + size error" );}
				  Dual<T> temp( *this );
		      temp.EPSILONS += d_plus.EPSILONS;
          temp.REAL += d_plus.REAL;
		      return temp;
			  }

        /// Binary -
        Dual<T> operator-( const Dual<T>& d_minus ) const
        {
          if ( d_minus.EPSILONS.size() != EPSILONS.size() )
          { throw Error( "Dual epsilon - size error" );}
          Dual<T> temp( *this );
          temp.EPSILONS -= d_minus.EPSILONS;
          temp.REAL -= d_minus.REAL;
          return temp;
        }

        /// Multiplication
        Dual<T> operator*( const Dual<T>& d_times ) const
        {
          if ( d_times.EPSILONS.size() != EPSILONS.size() )
          { throw Error( "Dual epsilon * size error" );}
          Dual<T> result;
          result.EPSILONS.resize( EPSILONS.size() );
          result.REAL = REAL * d_times.REAL;
          for ( std::size_t i=0; i<EPSILONS.size(); ++i)
          {
            result.EPSILONS[ i ] = REAL * d_times.EPSILONS[ i ] 
                                 + d_times.REAL * EPSILONS[ i ];
          }
          return result;
        }

        /// Scalar addition (both sides)
        Dual<T> operator+( const T& plus ) const
			  {
				  Dual<T> temp( *this );
          temp.REAL += plus;
		      return temp;
			  }
        // Friend function so the + operator can be on either side
        friend Dual<T> operator+( const T& plus, const Dual<T>& dual )
			  {
				  return dual + plus;
			  }	

        /// Scalar subtraction
        Dual<T> operator-( const T& minus ) const
        {
          Dual<T> temp( *this );
          temp.REAL -= minus;
          return temp;
        }

        /// Scalar multiplication
        Dual<T> operator*( const T& times ) const
        {
          Dual<T> temp( *this );
          temp.REAL *= times;
          temp.EPSILONS = temp.EPSILONS * times;
          return temp;
        }
        // Friend function so the * operator can be on either side
        friend Dual<T> operator*( const T& times, const Dual<T>& dual )
        {
          return dual * times;
        }

        /// Scalar division
        Dual<T> operator/( const T& div ) const
        {
          Dual<T> temp( *this );
          temp.REAL /= div;
          temp.EPSILONS = temp.EPSILONS / div;
          return temp;
        } 

        /* ----- Methods ----- */

        /// Return the real part of the dual number
        double real() const { return REAL; }

        /// Return a pointer to the real part of the dual number
        double& real() { return REAL; }

        /// Return the Vector of epsilons
        Vector<double> epsilons() const { return EPSILONS; }

        /// Return a pointer to the Vector of epsilons
        Vector<double>& epsilons() { return EPSILONS; } 

	}; // End of class Dual

  /// Output operator <<
  template <class Type>
  inline std::ostream& operator<<( std::ostream& os, const Dual<Type>& dual )
  {
    os << dual.REAL << "\t[";
    for ( std::size_t i=0; i<dual.EPSILONS.size() - 1; ++i )
    {
      os << dual.EPSILONS[ i ] << ", ";
    }
    os << dual.EPSILONS[ dual.EPSILONS.size() - 1 ] << "]";
    return os;
  } 

} // End of namespace TSL

#endif
