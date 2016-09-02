/* Matrix class - Encapsulates the Eigen matrix container
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <Eigen/Dense>

#include "Error.h"

namespace TSL
{
	// Templated matrix class
	template <class T>

	
	class Matrix 
	{
		
		protected:
				std::size_t ROWS;					            // Number of rows 
				std::size_t COLS;					            // Number of columns
				Eigen::Matrix<T, -1, -1> MATRIX;	    // Dynamic matrix object

		public:
		  /// Constructor for an empty matrix of unspecified size
			Matrix() : ROWS( 0 ), COLS( 0 ) { }

      /// Copy constructor
			Matrix( const Matrix<T>& source )
			{
				ROWS = source.ROWS;
				COLS = source.COLS;
				MATRIX = source.MATRIX;
			}

      // Constructor for a matrix with specfied initial elements
			Matrix( const std::size_t& rows, const std::size_t& cols, const T& elem )
      {
        ROWS = rows;
        COLS = cols;
        MATRIX.resize( ROWS, COLS );
        MATRIX.fill( elem );
      }

			/// Destructor
	   	~Matrix(){ }

      /* ----- Operator overloading ----- */

      /// Output operator <<
			template <class Type>
			friend std::ostream& operator<<( std::ostream& os, const Matrix<Type>& mat );

      /// Indexing operator ( read only )
			const T& operator() ( const std::size_t& i, const std::size_t& j ) const
			{
					// Range check 
					if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
					if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}
					return this->MATRIX( i, j );
			}

      /// Indexing operator ( read/write )
			T& operator() ( const std::size_t& i, const std::size_t& j )
			{
					// Range check 
					if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
					if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}
					return this->MATRIX( i, j );
			}

      /// Assignment
    	Matrix& operator=( const Matrix& original )
      {		
		    MATRIX = original.MATRIX;
		    ROWS = original.ROWS;
		    COLS = original.COLS;
		    return *this;
	    }

      /// Unary +
      Matrix<T> operator+() const
			{
				return *this;
			}      

      /// Unary -
      Matrix<T> operator-() const
      {
        Matrix<T> temp( *this );
		    temp.MATRIX = -temp.MATRIX;
		    return temp;
      }

      /// Binary +
      Matrix<T> operator+( const Matrix<T>& m_plus ) const
			{
				if ( m_plus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_plus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
				Matrix<T> temp( *this );
		    temp.MATRIX += m_plus.MATRIX;
		    return temp;
			}

      /// Binary -
      Matrix<T> operator-( const Matrix<T>& m_minus ) const
			{
				if ( m_minus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_minus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
				Matrix<T> temp( *this );
		    temp.MATRIX -= m_minus.MATRIX;
		    return temp;
			}

      /// Scalar multiplication (both sides)
      Matrix<T> operator*( const T& m_times )const
      {
        Matrix<T> temp( *this );
		    for (size_t j = 0; j < temp.COLS; ++j)
		    {
			    for (size_t i = 0; i < temp.ROWS; ++i)
			    {
				    temp.MATRIX( i, j ) *= m_times;
			    }
		    }
		    return temp;
      }
      // Friend function so the * operator can be on either side
			friend Matrix<T> operator*( const T& m, Matrix<T>& mat )
			{
				return mat * m;
			}	

      /// Scalar division 
      Matrix<T> operator/( const T& m_div )const
      {
        Matrix<T> temp( *this );
		    for (size_t j = 0; j < temp.COLS; ++j)
		    {
			    for (size_t i = 0; i < temp.ROWS; ++i)
			    {
				    temp( i, j ) /= m_div;
			    }
		    }
		    return temp;
      }

      /// Matrix multiplication
      Matrix<T> operator*( const Matrix<T>& B ) const
      {
        if ( COLS != B.ROWS ) { throw Error( "Matrix error: dimensions do not agree." );}
        Matrix<T> result;
				result.MATRIX.resize( ROWS , B.COLS );
				result.ROWS = ROWS;
				result.COLS = B.COLS;
        result.MATRIX = MATRIX * B.MATRIX;
        return result;
      }

      /// Addition assignment
      Matrix<T> operator+=( const Matrix<T>& m_plus )
      {
        if ( m_plus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_plus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
        MATRIX = MATRIX + m_plus.MATRIX;
        return *this;
      }

      /// Subtraction assignment
      Matrix<T> operator-=( const Matrix<T>& m_minus )
      {
        if ( m_minus.ROWS != ROWS ) { throw Error( "Matrix error: dimension 1 " );}
				if ( m_minus.COLS != COLS ) { throw Error( "Matrix error: dimension 2 " );}
        MATRIX = MATRIX - m_minus.MATRIX;
        return *this; 
      }

      /// Scalar multiplication assignment TODO
      
      /// Scalar division assignment TODO

      /* ----- Methods ----- */

      /// Return the number of rows in the matrix
      std::size_t rows() const { return ROWS; }

      /// Return the number of columns in the matrix
      std::size_t cols() const { return COLS; }

      /// Return the number of elements in the matrix
      std::size_t numel() const { return ROWS * COLS; }

      /// Return the transpose of the matrix
      Matrix<T> transpose() const
      {
        Matrix<T> temp( *this );
        Matrix<T> transp;
        transp.MATRIX = temp.MATRIX.transpose();
        transp.ROWS = temp.COLS;
        transp.COLS = temp.ROWS;
        return transp;
      }

      /// Transpose the matrix in place
      void transpose_in_place()
      {
        const std::size_t m = COLS;
        MATRIX.transposeInPlace();
        COLS = ROWS;
        ROWS = m;
      }

      /// Return the adjoint of the matrix
      Matrix<T> adjoint() const
      {
        Matrix<T> temp( *this );
        Matrix<T> adj;
        adj.MATRIX = temp.MATRIX.adjoint();
        adj.ROWS = temp.COLS;
        adj.COLS = temp.ROWS;
        return adj;
      }

      /// Return the conjugate of the matrix
      Matrix<T> conjugate() const
      {
         Matrix<T> temp( *this );
         Matrix<T> conj;
         conj.MATRIX = temp.MATRIX.conjugate();
         conj.ROWS = temp.ROWS;
         conj.COLS = temp.COLS;
         return conj;
      }

      /// Resize the matrix
      void resize( const std::size_t& rows, const std::size_t& cols )
      {
        MATRIX.conservativeResize( rows, cols );
        ROWS = rows;
        COLS = cols;
      }

      /// Swap the selected rows of the matrix TODO

      /// Fill the matrix with specified elements
      void fill( const T& elem )
      {
        MATRIX.fill( elem );
      }

      /// Fill the leading diagonal with specified elements
      void fill_diag( const T& elem )
      {
        std::size_t N( ROWS );
        if ( COLS < ROWS ) { N = COLS; }
        for ( std::size_t i=0; i<N; ++i )
        {
          MATRIX( i, i ) = elem;
        }        
      }

      /// Fill a diagonal band of the matrix
      void fill_band( const std::size_t& offset, const T& value )
      {
        for ( std::size_t row = 0; row < ROWS; ++row )
				{
            if ( ( row + offset < COLS ) && ( row + offset >= 0 ) )
            {
                MATRIX( row, row + offset ) = value;
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
                MATRIX( row, row - 1 ) = L;
            }
        }
        // Fill the main diagonal
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
            if ( ( row < COLS ) && ( row >= 0 ) )
            {
                MATRIX( row, row ) = D;
            }
        }
        // Fill the upper band
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
            if ( ( row + 1 < COLS ) && ( row + 1 >= 0 ) )
            {
                MATRIX( row, row + 1 ) = U;
            }
        }
      }

      /* ----- Norms -----  */

      /// Return the maximum absolute column sum of the matrix
      double norm_1() const
			{
        double max( 0.0 );
        for ( std::size_t j = 0; j < COLS; ++j )
        {
            double sum( 0.0 );
            for ( std::size_t i = 0; i < ROWS; ++i )
            {
                sum += std::abs( MATRIX( i, j ) );
            }
            max = std::max( max, sum );
        }
        return max;
    	}

      /// Return the maximum absolute row sum of the matrix
      double norm_inf() const
			{
        double max( 0.0 );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
            double sum( 0.0 );
            for ( std::size_t j = 0; j < COLS; ++j )
            {
                sum += std::abs( MATRIX( i, j ) );
            }
            max = std::max( max, sum );
        }
        return max;
    	}

      /// Return the entrywise p-norm of the matrix (p=2 is Frobenius, p=inf is max norm)
      double norm_p( const double& p ) const
			{
        double sum( 0.0 );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
            for ( std::size_t j = 0; j < COLS; ++j )
            {
                sum += std::pow( std::abs( MATRIX( i, j ) ), p );
            }
        }
        return std::pow( sum , 1.0/p );
    	}

      /// Return the Frobenius norm of the matrix
      double norm_frob() const
			{
				return this->norm_p( 2.0 );
			}	
				
      /// Return the entrywise max-norm of the matrix
      double norm_max() const
			{
        double max( 0.0 );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
            for ( std::size_t j = 0; j < COLS; ++j )
            {
                max = std::max( max, std::abs( MATRIX( i, j ) ) );
            } 
        }
        return max;
    	}


	}; // End of class Matrix

  /// Output operator <<
  template <class Type>
  inline std::ostream& operator<<( std::ostream& os, const Matrix<Type>& mat )
  {
    os << mat.MATRIX;
    return os;
  }
} // End of namespace TSL

#endif
