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
#include "SparseVector.h"

#if defined(PETSC_D) || defined(PETSC_Z)
#include "petscsys.h"
#include "petsc.h"
#include "petscksp.h"
#include "mpi.h"
#endif

#ifdef SLEPC
#include "slepceps.h"
#endif

namespace TSL
{

  template <class T>

  class SparseMatrix
  {
    //typedef std::map<std::size_t, std::map<std::size_t , T> > s_mat;
    //typedef typename s_mat::const_iterator row_iter;
    //typedef std::map<std::size_t, T> s_col;
    //typedef typename s_col::const_iterator col_iter;
    typedef Eigen::SparseMatrix<T, Eigen::ColMajor, long long> SpMat;
    typedef Eigen::Triplet<T> Triplet;

    typedef typename std::map< std::size_t, T >::const_iterator citer;
    typedef typename std::map< std::size_t, T >::iterator iter;

    private:

      std::size_t ROWS;                             // Number of rows
      std::size_t COLS;                             // Number of columns
      std::vector< SparseVector<T> > MATRIX;        // STL vector of SparseVectors

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
				//CONTAINER = source.CONTAINER;
        ROWS = source.ROWS;
        COLS = source.COLS;
        MATRIX.reserve( ROWS );
        MATRIX = source.MATRIX;
			}

      /// Constructor for a matrix of specified size
			SparseMatrix( const std::size_t& rows, const std::size_t& cols )
      {
        ROWS = rows;
        COLS = cols;
        MATRIX.reserve( ROWS );
        SparseVector<T> sparse_row( COLS );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
          MATRIX.push_back( sparse_row );
        }
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

      /// Return a SparseVector containing the specified row
      SparseVector<T> get_row( const std::size_t& row ) const
      {
        return MATRIX[ row ];
      }

      /// Set a specified row of the SparseMatrix using a SparseVector
      void set_row( const std::size_t& row, const SparseVector<T>& row_vector )
      {
        MATRIX[row] = row_vector;
      }

      /// Return the number of rows in the matrix
      std::size_t rows() const { return ROWS; }

      /// Return the number of columns in the matrix
      std::size_t cols() const { return COLS; }

      /// Return the size of the matrix
      std::size_t size() const { return ROWS * COLS; }

      /// Return the number of non zero elements
      std::size_t numel() const
      {
        std::size_t nelts( 0 );
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
          nelts += MATRIX[ row ].nelts();
        }
        return nelts;
      }

      /// Return the number of non zero elements in a specified row
      std::size_t numel_row( std::size_t row ) const
      {
        return MATRIX[ row ].nelts();
      }

      /// Clear all the elements
      void clear()
      {
        //CONTAINER.clear();
        MATRIX.clear();
        MATRIX.reserve( ROWS );
        SparseVector<T> sparse_row( COLS );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
          MATRIX.push_back( sparse_row );
        }
      }

      /// Set the matrix to the identity matrix
      void eye()
      {
        this->clear();
        for ( std::size_t i=0; i < ROWS; ++i )
        {
          MATRIX[ i ][ i ] = 1.0;
        }
      }

      /// Scale the matrix by m
      void scale( const T& m )
      {
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
          for ( citer pos = MATRIX[row].begin(); pos != MATRIX[row].end(); ++pos )
          {
            std::size_t col = pos -> first;
            T elem = pos -> second;
            MATRIX[ row ][ col ] = elem * m;
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

      /// Print the matrix to screen
      void print()
      {
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
          for ( citer pos = MATRIX[row].begin(); pos != MATRIX[row].end(); ++pos )
          {
            //std::size_t col = pos -> first;
            //T elem = pos -> second;

            std::cout << "( " << row << ", "
                      << pos -> first << ", "
                      << pos -> second << " ), ";
          }
        }
        std::cout << std::endl;
      }

      /// Output the matrix to a file
      void output( std::string filename, int precision = 10 ) const //TODO
      {
        std::ofstream dump;
        dump.open( filename.c_str() );
        dump.precision( precision );
        dump.setf( std::ios::showpoint );
        dump.setf( std::ios::showpos );
        //dump.setf( std::ios::scientific );
        for ( std::size_t row = 0; row < ROWS; ++row )
        {
          dump << "row " << row << " : ";
          for ( citer pos = MATRIX[row].begin(); pos != MATRIX[row].end(); ++pos )
          {
            dump << "[ " << pos -> first << " ] = "
                 << pos -> second << " ";
          }
          dump << std::endl;
        }
        dump << std::endl;
        dump.close();
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

        for ( std::size_t row = 0; row < ROWS; ++row )
        {
          for ( citer pos = MATRIX[row].begin(); pos != MATRIX[row].end(); ++pos )
          {
            std::size_t col = pos -> first;
            T elem = pos -> second;
            tripletList.push_back( Triplet( row, col, elem  ) );
          }
        }


        S_mat.setFromTriplets( tripletList.begin(), tripletList.end() );
        S_mat.makeCompressed();

        return S_mat;
      }

#if defined (PETSC_D) || defined (PETSC_Z)
    /// Converts the SparseMatrix to a compressed format for a specified row.
    void get_row_petsc( PetscInt row_number, PetscScalar* storage, PetscInt* cols );

    /// Extracts the number of non-zero elements in each row and returns a PetscInt array.
    void nelts_all_rows( PetscInt* row_nnz)
    {
      for ( std::size_t i = 0; i < ROWS; ++i )
      {
        row_nnz[i] = this->numel_row( i );
      }
    }

#endif

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
    return MATRIX[ i ][ j ];
  }

  /// Assignment
  template <class T>
  inline SparseMatrix<T>& SparseMatrix<T>::operator=( const SparseMatrix<T>& original )
  {
    //CONTAINER = original.CONTAINER;
    MATRIX = original.MATRIX;
    ROWS = original.ROWS;
    COLS = original.COLS;
    return *this;
  }
} // End of namespace TSL

#endif
