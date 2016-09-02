/* SparseMatrix class - Encapsulates the Eigen sparse matrix container
*/

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "Error.h"
#include "Vector.h"

namespace TSL
{
	// Templated sparse matrix class
	template <class T>

	
	class SparseMatrix 
	{
		
		protected:
      Eigen::SparseMatrix<T> S_MATRIX;      // Eigen sparse matrix container
      

		public:
		  /// Constructor for an empty matrix of unspecified size
			SparseMatrix(){ }

      /// Copy constructor
			SparseMatrix( const SparseMatrix<T>& source )
			{
				S_MATRIX = source.S_MATRIX;
			}

      /// Constructor for a matrix of specified size
			SparseMatrix( const std::size_t& rows, const std::size_t& cols )
      {
        S_MATRIX.resize( rows, cols );
      }

			/// Destructor
	   	~SparseMatrix(){ }

      /* ----- Operator overloading ----- */

      /// Indexing operator (read only)
      const T operator() ( const std::size_t& i, const std::size_t& j ) const
      {
					// Range check 
        if ( i<0 || this->rows()<=i )	{ throw Error( "Matrix range error: dimension 1" );}
        if ( j<0 || this->cols()<=j )	{ throw Error( "Matrix range error: dimension 2" );}
        return this->get_elem( i, j );
      }

      /// Indexing operator (read/write)
      T& operator() ( const std::size_t& i, const std::size_t& j )
      {
        // Range check 
        if ( i<0 || this->rows()<=i )	{ throw Error( "Matrix range error: dimension 1" );}
        if ( j<0 || this->cols()<=j )	{ throw Error( "Matrix range error: dimension 2" );}
        return this->value_ref( i, j );
      } 
 
      /// Assignment
    	SparseMatrix& operator=( const SparseMatrix& original )
      {		
		    S_MATRIX = original.S_MATRIX;
		    return *this;
	    } 

      /// Binary +
      SparseMatrix<T> operator+( const SparseMatrix<T>& m_plus ) const
			{
				if ( m_plus.rows() != this->rows() ) {throw Error( "Matrix error: dimension 1 " );}
				if ( m_plus.cols() != this->cols() ) {throw Error( "Matrix error: dimension 2 " );}
				SparseMatrix<T> temp( *this );
		    temp.S_MATRIX += m_plus.S_MATRIX;
		    return temp;
			}

      /// Binary +
      SparseMatrix<T> operator-( const SparseMatrix<T>& m_minus ) const
			{
				if ( m_minus.rows() != this->rows() ){throw Error( "Matrix error: dimension 1 " );}
				if ( m_minus.cols() != this->cols() ){throw Error( "Matrix error: dimension 2 " );}
				SparseMatrix<T> temp( *this );
		    temp.S_MATRIX -= m_minus.S_MATRIX;
		    return temp;
			}

      /* ----- Methods ----- */
      
      /// Return the number of rows in the matrix
      std::size_t rows() const { return S_MATRIX.rows(); }

      /// Return the number of columns in the matrix
      std::size_t cols() const { return S_MATRIX.cols(); }

      /// Return the size of the matrix
      std::size_t size() const { return S_MATRIX.size(); }

      /// Return the number of non zero elements
      std::size_t numel() const { return S_MATRIX.nonZeros(); }

      /// Set the matrix to the identity matrix
      void eye() { S_MATRIX.setIdentity(); } 

      /// Set element in the matrix
      void set_elem( const std::size_t i, const std::size_t j, const T& elem )
      {
        S_MATRIX.coeffRef( i, j ) = elem;
      }

      /// Get element in the matrix
      const T get_elem( const std::size_t i, const std::size_t j ) const
      {
        return S_MATRIX.coeff( i, j );
      }

      /// Reference to the value of the matrix at position (i,j)
      T& value_ref( const std::size_t i, const std::size_t j ) 
      {
        return S_MATRIX.coeffRef( i, j );
      }

      /// Scale the matrix by m
      void scale( const T& m ) { S_MATRIX = S_MATRIX * m; }

      /// Scalar multiplication (both sides)
      SparseMatrix<T> operator*( const T& m_times )const
      {
        SparseMatrix<T> temp( *this );
		    temp.scale( m_times );
		    return temp;
      }
      // Friend function so the * operator can be on either side
			friend SparseMatrix<T> operator*( const T& m, SparseMatrix<T>& mat )
			{
				return mat * m;
			}

      /* ----- Solve sparse linear systems ----- */
      
      /// Solve system of equations Ax=b where x and b are vectors
      Vector<T> solve( const Vector<T>& b ) const
      {
        if ( this->rows() != b.SIZE ) {throw Error( "Sparse solver error: dimension 1 " );}
        Eigen::SparseLU< Eigen::SparseMatrix<T> > solver;
        Eigen::SparseMatrix<T> S_mat( S_MATRIX );
        S_mat.makeCompressed(); 
        solver.compute( S_mat );
        Eigen::Matrix<T, -1, 1> X;
	      X.resize( b.SIZE, 1 );
        X = solver.solve( b.VECTOR );
        Vector<T> x( b.SIZE );
        x.VECTOR = X;
        return x;
      }
       
	}; // End of class SparseMatrix

} // End of namespace TSL

#endif
