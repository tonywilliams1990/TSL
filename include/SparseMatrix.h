/* Sparse Matrix class
*/

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>
#include <map>
#include <vector>
#include <cstdlib>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>


#include "Error.h"
#include "Vector.h"

namespace TSL
{

  template <class T>

  class SparseMatrix
  {
    typedef std::map<std::size_t, std::map<std::size_t , T> > s_mat;
    typedef typename s_mat::const_iterator row_iter;
    typedef std::map<std::size_t, T> s_col;
    typedef typename s_col::const_iterator col_iter;
    typedef Eigen::SparseMatrix<T, Eigen::ColMajor, long long> SpMat;
    typedef Eigen::Triplet<T> Triplet;

    private:
      s_mat CONTAINER;                              // Sparse container
      std::size_t ROWS;                             // Number of rows
      std::size_t COLS;                             // Number of columns

    public:

      /// Constructor for an empty matrix of unspecified size
			SparseMatrix()
      {
        ROWS = 0;
        COLS = 0;
      }

      /// Copy constructor
			SparseMatrix( const SparseMatrix<T>& source )
			{
				CONTAINER = source.CONTAINER;
        ROWS = source.ROWS;
        COLS = source.COLS;
			}

      /// Constructor for a matrix of specified size
			SparseMatrix( const std::size_t& rows, const std::size_t& cols )
      {
        ROWS = rows;
        COLS = cols;
      }

			/// Destructor
	   	~SparseMatrix(){ }

      /* ----- Operator overloading ----- */

      /// Indexing operator (read only)
      const T operator() ( const std::size_t& i, const std::size_t& j ) const;

      /// Indexing operator (read/write)
      T& operator() ( const std::size_t& i, const std::size_t& j );

      /// Assignment
    	SparseMatrix<T>& operator=( const SparseMatrix<T>& original );

      /* ----- Methods ----- */

      /// Return the number of rows in the matrix
      std::size_t rows() const { return ROWS; }

      /// Return the number of columns in the matrix
      std::size_t cols() const { return COLS; }

      /// Return the size of the matrix
      std::size_t size() const { return ROWS * COLS; }

      /// Return the number of non zero elements
      std::size_t numel() const
      {
        std::size_t sum( 0 );
        row_iter ii;
        // Loop over the rows and add up the size of the columns
        for(ii=this->CONTAINER.begin(); ii!=this->CONTAINER.end(); ii++)
        {
          sum += (*ii).second.size();
        }
        return sum;
      }

      /// Clear all the elements
      void clear()
      {
        CONTAINER.clear();
      }

      /// Set the matrix to the identity matrix
      void eye()
      {
        CONTAINER.clear();
        for ( std::size_t i=0; i < ROWS; ++i )
        {
          CONTAINER[ i ][ i ] = 1.0;
        }
      }

      /// Scale the matrix by m
      void scale( const T& m )
      {
        row_iter ii;
        col_iter jj;
        for(ii=this->CONTAINER.begin(); ii!=this->CONTAINER.end(); ii++)
        {
          for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
          {
            std::size_t i = (*ii).first;
            std::size_t j = (*jj).first;
            T elem = (*jj).second;
            CONTAINER[ i ][ j ] = elem * m ;
          }
        }
      }

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

      /// Print the matrix
      void print()
      {
        row_iter ii;
        col_iter jj;
        for(ii=this->CONTAINER.begin(); ii!=this->CONTAINER.end(); ii++)
        {
          for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
          {
            std::cout << "( " << (*ii).first << ", "
                      << (*jj).first << ", "
                      << (*jj).second << " ), ";
          }
        }
        std::cout << std::endl;

      }

      /* ----- Solve sparse linear systems ----- */

      /// Solve system of equations Ax=b where x and b are vectors (SparseLU)
      Vector<T> solve( const Vector<T>& b ) const
      {
        if ( ROWS != b.SIZE ) {throw Error( "Sparse solver error: dimension 1 " );}
        Eigen::SparseLU< SpMat > solver;
        SpMat S_mat( ROWS, COLS );
        // Fill Eigen sparse matrix
        S_mat = this->convert_to_Eigen();
        /*std::vector<Triplet> tripletList;             // Vector of triplets
        row_iter ii;
        col_iter jj;
        for(ii=this->CONTAINER.begin(); ii!=this->CONTAINER.end(); ii++)
        {
          for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
          {
            std::size_t i = (*ii).first;
            std::size_t j = (*jj).first;
            T elem = (*jj).second;
            tripletList.push_back( Triplet(i, j, elem) );
          }
        }
        S_mat.setFromTriplets( tripletList.begin(), tripletList.end() );
        S_mat.makeCompressed();*/
        solver.compute( S_mat );

        Eigen::Matrix<T, -1, 1> X;
	      X.resize( b.SIZE, 1 );
        X = solver.solve( b.VECTOR );
        Vector<T> x( b.SIZE );
        x.VECTOR = X;
        return x;
      }

      /// Return the Eigen::SparseMatrix of the TSL::SparseMatrix
      SpMat convert_to_Eigen() const
      {
        SpMat S_mat( ROWS, COLS );
        // Fill Eigen sparse matrix
        std::vector<Triplet> tripletList;             // Vector of triplets
        row_iter ii;
        col_iter jj;
        for(ii=this->CONTAINER.begin(); ii!=this->CONTAINER.end(); ii++)
        {
          for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
          {
            std::size_t i = (*ii).first;
            std::size_t j = (*jj).first;
            T elem = (*jj).second;
            tripletList.push_back( Triplet(i, j, elem) );
          }
        }
        S_mat.setFromTriplets( tripletList.begin(), tripletList.end() );
        S_mat.makeCompressed();

        return S_mat;
      }

  }; // End of class SparseMatrix

  /* Indexing operator (read only)
  template <class T>
  inline const T SparseMatrix<T>::operator() ( const std::size_t& i,
                                               const std::size_t& j ) const
  {
    // Range check
    if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
    if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}
    return CONTAINER[ i ][ j ];
  }*/

  /// Indexing operator (read/write)
  template <class T>
  inline T& SparseMatrix<T>::operator() ( const std::size_t& i, const std::size_t& j )
  {
    // Range check
    if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
    if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}
    return CONTAINER[ i ][ j ];
  }

  /// Assignment
  template <class T>
  inline SparseMatrix<T>& SparseMatrix<T>::operator=( const SparseMatrix<T>& original )
  {
    CONTAINER = original.CONTAINER;
    ROWS = original.ROWS;
    COLS = original.COLS;
    return *this;
  }


} // End of namespace TSL

#endif
