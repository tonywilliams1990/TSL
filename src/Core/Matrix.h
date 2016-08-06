/* Matrix.h - Defines a matrix class using an STL vector as a container 
   Created  - 06/08/2016
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Error.h"

namespace TSL
{

	//-----------------------------------------------------------------------------------------------

	/// Function objects for various apply() operations:

	template<class T> struct Assign { void operator()(T& a, const T& c) { a = c; } };

	template<class T> struct Add_assign { void operator()(T& a, const T& c) { a += c; } };

	template<class T> struct Mul_assign { void operator()(T& a, const T& c) { a *= c; } };
	
	template<class T> struct Minus_assign { void operator()(T& a, const T& c) { a -= c; } };

	template<class T> struct Div_assign { void operator()(T& a, const T& c) { a /= c; } };

	//-----------------------------------------------------------------------------------------------
	
	/// A templated matrix class
	template<class T> class Matrix
	{
		protected:
			std::size_t ROWS;						// Number of rows
			std::size_t COLS;						// Number of columns
			std::vector<T> CONTAINER;		// Container for storing the elements

		public:
			
			/// Constructor for an empty matrix of unspecified size
			Matrix() : ROWS( 0 ), COLS( 0 ) { }

			/// Constructor for a matrix of specified size and elements
			Matrix( const std::size_t& rows, const std::size_t& cols, const T& elem ) : ROWS( rows ), COLS( cols )
			{
				CONTAINER.assign( ROWS * COLS, elem );
			}

			/// Copy constructor
			Matrix( const Matrix<T>& source )
			{
				ROWS = source.ROWS;
				COLS = source.COLS;
				CONTAINER = source.CONTAINER;
			}

			/// Destructor
			~Matrix(){ }

			/* ----- Operator overloading ----- */
			
			/// Indexing operator ( read only )
			const T& operator() ( const std::size_t& i, const std::size_t& j ) const
			{
					// Range check 
					if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
					if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}
					
					return this->CONTAINER[ i*COLS + j ];
			}
		
			/// Indexing operator ( read/write )
			T& operator() ( const std::size_t& i, const std::size_t& j )
			{
					// Range check 
					if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
					if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}

					return this->CONTAINER[ i*COLS + j ];
			}

			/// Output operator <<
			template <class Type>
			friend std::ostream& operator<<( std::ostream& os, const Matrix<Type>& m );	

			/// Element-wise operations:

			/// Apply functions
			template<class F>	void apply(F f)
			{
				for (std::size_t i=0; i<CONTAINER.size(); ++ i ) { f( CONTAINER[i] ); }
			}	

			template<class F>	void apply(F f, const T& c)
			{
				for (std::size_t i=0; i<CONTAINER.size(); ++ i ) { f( CONTAINER[i], c ); }
			}	

			/*----- Element-wise operations -----*/

			Matrix& operator=(const T& c)  { this->apply(Assign<T>(),c);   			return *this; }
			Matrix& operator+=(const T& c) { this->apply(Add_assign<T>(),c);   	return *this; }
			Matrix& operator*=(const T& c) { this->apply(Mul_assign<T>(),c);   	return *this; }
			Matrix& operator-=(const T& c) { this->apply(Minus_assign<T>(),c); 	return *this; }
			Matrix& operator/=(const T& c) { this->apply(Div_assign<T>(),c);   	return *this; }
			

	};	// End of class Matrix

	/// Output operator
	template<class T> inline std::ostream& operator<<( std::ostream& os, const Matrix<T>& m)
	{
		os << std::endl;
		// Scientific output with 3 decimal points precision
		os << std::scientific << std::setprecision(3);

		for (std::size_t i = 0; i<m.ROWS; ++i)
		{
			for (std::size_t j=0; j<m.COLS; ++j)
			{
				os << m(i,j) << "\t";
			}
			std::cout << std::endl;
		}

		return os;
	}

}  // End of namespace TSL

#endif
