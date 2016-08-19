/* Matrix.h - Defines a matrix class using an STL vector as a container 
   Created  - 06/08/2016
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <iterator>

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

			typedef typename std::vector<T>::iterator iter;
			
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

			/// Binary +
			Matrix<T> operator+( const Matrix<T>& m_plus ) const
			{
				if ( m_plus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_plus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
				Matrix<T> result( *this );
				std::transform( CONTAINER.begin(), CONTAINER.end(), m_plus.CONTAINER.begin(),
								result.CONTAINER.begin(), std::plus<T>() );
				return result;
			}

			/// Binary -
			Matrix<T> operator-( const Matrix<T>& m_minus ) const
			{
				if ( m_minus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_minus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
				Matrix<T> result( *this );
				std::transform( CONTAINER.begin(), CONTAINER.end(), m_minus.CONTAINER.begin(),
								result.CONTAINER.begin(), std::minus<T>() );
				return result;
			}

			/// Unary +
			Matrix<T> operator+() const
			{
				return *this;
			}

			/// Unary -
			Matrix<T> operator-() const
			{
				Matrix<T> result( *this );
				std::transform (CONTAINER.begin(), CONTAINER.end(),result.CONTAINER.begin(),
						 std::negate<T>());
				return result;
			}

			/// Matrix multiplication A*B
			Matrix<T> operator*( const Matrix<T>& B ) const
			{
				if ( COLS != B.ROWS ) { throw Error( "Matrix error (*): dimensions do not agree." );}
				Matrix<T> result;
				result.CONTAINER.resize( ROWS * B.COLS );
				result.ROWS = ROWS;
				result.COLS = B.COLS;
				std::size_t n( ROWS ), p( B.COLS ), m( COLS );
				Matrix<T> temp( B );
				temp.transpose_in_place();
				for ( std::size_t i=0; i<n; ++i )
				{
					for ( std::size_t j=0; j<p; ++j )
					{
						for ( std::size_t k=0; k<m; ++k )
						{
							result.CONTAINER[ i*p + j ] += CONTAINER[ i*m + k ] * temp.CONTAINER[ j*m + k ];
						}
					}
				} 
				return result;
			}

			/// Addition assignment +=
			Matrix<T>& operator+=( const Matrix<T>& m_plus )
			{
				if ( m_plus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_plus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
				std::transform (CONTAINER.begin(), CONTAINER.end(), m_plus.CONTAINER.begin(), 
						CONTAINER.begin(), std::plus<T>());
				return *this;
			}

			/// Subtraction assignment -=
			Matrix<T>& operator-=( const Matrix<T>& m_minus )
			{
				if ( m_minus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_minus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
				std::transform (CONTAINER.begin(), CONTAINER.end(), m_minus.CONTAINER.begin(), 
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

			Matrix& operator=(const T& c)  { this->apply(Assign<T>(),c);   			return *this; }
			Matrix& operator+=(const T& c) { this->apply(Add_assign<T>(),c);   	return *this; }
			Matrix& operator*=(const T& c) { this->apply(Mul_assign<T>(),c);   	return *this; }
			Matrix& operator-=(const T& c) { this->apply(Minus_assign<T>(),c); 	return *this; }
			Matrix& operator/=(const T& c) { this->apply(Div_assign<T>(),c);   	return *this; }

			/*----- Methods -----*/

			/// Return the number of rows in the matrix
			std::size_t rows() const
			{
				return ROWS;
			}

			/// Return the number of columns in the matrix
			std::size_t cols() const
			{
				return COLS;
			}

			/// Return the number of elements in the matrix
			std::size_t numel() const
			{
				return CONTAINER.size();
			}	

			/// Fill the matrix with specified elements (same as assign operator)
			void fill( const T& elem )
			{
				this->apply( Assign<T>(), elem);
			}

			/// Fill the leading diagonal with specified elements
			void fill_diag( const T& elem )
			{
				std::size_t N( ROWS );
				if ( COLS < ROWS )
				{
					N = COLS;
				}
				for (std::size_t i=0; i < N; ++i)
				{
					CONTAINER[ i*COLS + i ] = elem;
				}
			}

			/// Fill a diagonal band of the matrix
			void fill_band( const std::size_t& offset, const T& value )
			{
				for ( std::size_t row = 0; row < ROWS; ++row )
				{
            if ( ( row + offset < COLS ) && ( row + offset >= 0 ) )
            {
                CONTAINER[ row*COLS + row + offset ] = value;
            }
				}
			}

			/// Fill the main three diagonals of the matrix
      void fill_tridiag( const T& L, const T& D, const T& U )
			{
				// Fill the lower band
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
            if ( ( row - 1 < COLS ) && ( row - 1 >= 0 ) )
            {
                CONTAINER[ row*COLS + row - 1 ] = L;
            }
        }
        // Fill the main diagonal
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
            if ( ( row < COLS ) && ( row >= 0 ) )
            {
                CONTAINER[ row*COLS + row ] = D;
            }
        }
        // Fill the upper band
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
            if ( ( row + 1 < COLS ) && ( row + 1 >= 0 ) )
            {
                CONTAINER[ row*COLS + row + 1 ] = U;
            }
        }
			}

			/// Transpose the matrix in place
			void transpose_in_place()
			{
				const std::size_t m = COLS;
				iter first = CONTAINER.begin();
				iter last  = CONTAINER.end();
				const std::size_t mn1 = last - first - 1;
				const std::size_t n = (last - first) / m;
				std::vector<bool> visited( last - first );
				iter cycle = first;
				
				while (++cycle != last )
				{
					if	(visited[cycle - first])
						continue;
					std::size_t a = cycle -first;
					do {
							a = a == mn1 ? mn1 : (n * a) % mn1;
							std::swap(*(first + a), *cycle);
							visited[a] = true;
					} while ((first+a) != cycle);
				}
				COLS = ROWS;
				ROWS = m;				
			}

			/// Return the transpose of a matrix
			Matrix<T> transpose() const
			{
				Matrix<T> temp( *this );
				temp.transpose_in_place();
				return temp;
			}

			/// Resize the matrix ( clears all elements )
      void resize( const std::size_t& rows, const std::size_t& cols )
			{
				CONTAINER.clear();
				CONTAINER.resize( rows * cols );
				ROWS = rows;
				COLS = cols;
			}

			/// Swap the selected rows of the matrix
			void swap_rows( const std::size_t& row_1, const std::size_t& row_2 )
			{
				if ( row_1<0 || ROWS<=row_1 )	{ throw Error( "Matrix range error: swap rows" );}
				if ( row_2<0 || ROWS<=row_2 )	{ throw Error( "Matrix range error: swap rows" );}
				T temp;
				for (std::size_t k=0; k<COLS; ++k )
				{
					temp = CONTAINER[ row_1 * COLS + k ];
				 	CONTAINER[ row_1 * COLS + k ] = CONTAINER[ row_2 * COLS + k ];
					CONTAINER[ row_2 * COLS + k ] = temp;
				}
			}			

			/* ----- Norms ----- */

      // Return the maximum absolute column sum of the matrix
      double norm_1() const
			{
        double max( 0.0 );
        for ( std::size_t j = 0; j < COLS; ++j )
        {
            double sum( 0.0 );
            for ( std::size_t i = 0; i < ROWS; ++i )
            {
                sum += std::abs( CONTAINER[ i*COLS + j ] );
            }
            max = std::max( max, sum );
        }
        return max;
    	}

      // Return the maximum absolute row sum of the matrix
      double norm_inf() const
			{
        double max( 0.0 );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
            double sum( 0.0 );
            for ( std::size_t j = 0; j < COLS; ++j )
            {
                sum += std::abs( CONTAINER[ i*COLS + j ] );
            }
            max = std::max( max, sum );
        }
        return max;
    	}

      // Return the Frobenius norm of the matrix
      double norm_frob() const
			{
				return this->norm_p( 2.0 );
			}	

      // Return the entrywise max-norm of the matrix
      double norm_max() const
			{
        double max( 0.0 );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
            for ( std::size_t j = 0; j < COLS; ++j )
            {
                max = std::max( max, std::abs( CONTAINER[ i*COLS + j ] ) );
            } 
        }
        return max;
    	}

      // Return the entrywise p-norm of the matrix (p=2 is Frobenius, p=inf is max norm)
      double norm_p( const double& p ) const
			{
        double sum( 0.0 );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
            for ( std::size_t j = 0; j < COLS; ++j )
            {
                sum += std::pow( std::abs( CONTAINER[ i*COLS + j ] ), p );
            }
        }
        return std::pow( sum , 1.0/p );
    	}

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
