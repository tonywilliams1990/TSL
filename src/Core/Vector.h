/* Vector.h - Encapsulates the std::vector container
   adding operator overloading and a few useful functions.
	 Created  - 13/08/2016
*/

#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <functional>

#include "Error.h"
#include "Matrix.h"

namespace TSL
{
	
	/// A templated vector class
	template<class T> class Vector
	{
		protected:
			std::vector<T> CONTAINER;		// Container for storing the elements

		public:

			typedef typename std::vector<T>::iterator iter;

			/* ----- Constructors and destructor ----- */ 
			
			/// Constructor for an empty vector of unspecified size
			Vector() { }

			/// Constructor for a vector of specfied size
			Vector( const std::size_t& size )
			{
				CONTAINER.resize( size );
			}

			/// Constructor for an initialised vector of specified size
			Vector( const std::size_t& size, const T& elem )
			{
				CONTAINER.assign( size, elem );
			}

			/// Copy constructor
			Vector( const Matrix<T>& source )
			{
				CONTAINER = source.CONTAINER;
			}

			/// Destructor
			~Vector() { }

			/* ----- Operator overloading ----- */
			
			/// Indexing operator ( read only )
			const T& operator() ( const std::size_t& i ) const
			{
					// Range check 
					if ( i<0 || CONTAINER.size()<=i )	{ throw Error( "Vector range error." );}
					return this->CONTAINER[ i ];
			}

			/// Indexing operator ( read/write )
			T& operator() ( const std::size_t& i )
			{
					// Range check 
					if ( i<0 || CONTAINER.size()<=i )	{ throw Error( "Vector range error." );}
					return this->CONTAINER[ i ];
			}

			/// Output operator <<
			template <class Type>
			friend std::ostream& operator<<( std::ostream& os, const Vector<Type>& v );	

			/// Binary +
			Vector<T> operator+( const Vector<T>& v_plus ) const
			{
				if ( v_plus.size() != size() ) { throw Error( "(+): Vectors must be of equal length." );}
				Vector<T> result( v_plus.size() );
				std::transform (CONTAINER.begin(), CONTAINER.end(), v_plus.CONTAINER.begin(), 
						result.CONTAINER.begin(), std::plus<T>());
				return result;
			}

			/// Binary -
			Vector<T> operator-( const Vector<T>& v_minus ) const
			{
				if ( v_minus.size() != size() ) { throw Error( "(-): Vectors must be of equal length." );}
				Vector<T> result( v_minus.size() );
				std::transform (CONTAINER.begin(), CONTAINER.end(), v_minus.CONTAINER.begin(), 
						result.CONTAINER.begin(), std::minus<T>());
				return result;
			}

			/// Addition assignment +=
			Vector<T>& operator+=( const Vector<T>& v_plus )
			{
				if ( v_plus.size() != size() ) { throw Error( "(+=): Vectors must be of equal length." );}
				
				std::transform (CONTAINER.begin(), CONTAINER.end(), v_plus.CONTAINER.begin(), 
						CONTAINER.begin(), std::plus<T>());
				return *this;
			}

			/// Subtraction assignment -=
			Vector<T>& operator-=( const Vector<T>& v_minus )
			{
				if ( v_minus.size() != size() ) { throw Error( "(-=): Vectors must be of equal length." );}
				std::transform (CONTAINER.begin(), CONTAINER.end(), v_minus.CONTAINER.begin(), 
						CONTAINER.begin(), std::minus<T>());
				return *this;
			}

			/*----- Element-wise operations -----*/

			/// Apply functions
			template<class F>	void apply(F f)
			{
				for (std::size_t i=0; i<CONTAINER.size(); ++ i ) { f( CONTAINER[i] ); }
			}	

			template<class F>	void apply(F f, const T& c)
			{
				for (std::size_t i=0; i<CONTAINER.size(); ++ i ) { f( CONTAINER[i], c ); }
			}	

			/// Element-wise operations:

			Vector& operator=(const T& c)  { this->apply(Assign<T>(),c);   			return *this; }
			Vector& operator+=(const T& c) { this->apply(Add_assign<T>(),c);   	return *this; }
			Vector& operator*=(const T& c) { this->apply(Mul_assign<T>(),c);   	return *this; }
			Vector& operator-=(const T& c) { this->apply(Minus_assign<T>(),c); 	return *this; }
			Vector& operator/=(const T& c) { this->apply(Div_assign<T>(),c);   	return *this; }

			/// Scalar multiplication
			Vector<T> operator*( const T& m ) const
			{
				Vector<T> result( *this );
				result *= m;
				return result;
			}
			// Friend function so the operator can be on either side
			friend Vector<T> operator*( const T& m, Vector<T>& v )
			{
				return v * m;
			}	

			

			/* ----- Methods ----- */

			/// Size of the vector
			std::size_t size() const
			{
				return CONTAINER.size();
			}	

			/// Push_back a new element into the vector
			void push_back( const T& new_elem )
			{
				CONTAINER.push_back( new_elem );
			}

			/// Resize the vector so that it contains n elements
			void resize( const std::size_t& n )
			{
				CONTAINER.resize( n );
			}	

			// Clear all the elements from the vector
			void clear()
			{
				CONTAINER.clear();
			}	

			/// Assign new contents to the vector
			void assign( const std::size_t& n, const T& elem )
			{
				CONTAINER.assign( n, elem );
			} 

			/// Create a linearly spaced vector (of doubles) with n elements
			void linspace( const double& a, const double& b, const std::size_t& n );
	

			
			/* ----- Iterators ----- */

			/// Iterator pointing to the first element
			iter begin()
			{
				return CONTAINER.begin();
			}

			/// Iterator pointing to the last element
			iter end()
			{
				return CONTAINER.end();
			}

	};	// End of class Vector

	/// Output operator
	template<class T> inline std::ostream& operator<<( std::ostream& os, const Vector<T>& v)
	{
		// Scientific output with 3 decimal points precision
		os << std::scientific << std::setprecision(3);

		for (std::size_t i = 0; i<v.CONTAINER.size(); ++i)
		{
			os << v(i) << "  ";
		}
		return os;
	}

	/// Create a linearly spaced vector (of doubles) with n elements (specialised function)
	template<>
	void Vector<double>::linspace( const double& a, const double& b, const std::size_t& n )
	{
		CONTAINER.resize( n );
		const double h = ( b - a ) / (n - 1)  ;
		for ( std::size_t i=0; i < n; ++i ) 
		{		
			CONTAINER[ i ] = a + h * i;
		}
	}	 

}  // End of namespace TSL

#endif
