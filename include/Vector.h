/* Vector class - Encapsulates the Eigen vector container
*/

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "Error.h"

namespace TSL
{
  // Forward declare the Matrix class that will be friends with Vector
  template <class T>
  class Matrix;

  // Forward declare the SparseMatrix class that will be friends with Vector
  template <class T>
  class SparseMatrix;

	// Templated vector class
	template <class T>
	
	class Vector 
	{
		
      protected:
      friend class Matrix<T>;			        // Friend class that may access VECTOR directly
      friend class SparseMatrix<T>;			  // Friend class that may access VECTOR directly
      std::size_t SIZE;                   // Number of elements in the vector
      Eigen::Matrix<T, -1, 1> VECTOR;	    // Dynamic column vector

      public:
      /// Constructor for an empty matrix of unspecified size
      Vector() : SIZE( 0 ) { }

      /// Constructor for a matrix of specified size
      Vector( const std::size_t& size )
      {
        SIZE = size;
        VECTOR.resize( SIZE, 1 );
      }

      /// Constructor for a matrix of specified size with specified elements
      Vector( const std::size_t& size, const T& elem )
      {
        SIZE = size;
        VECTOR.resize( SIZE, 1 );
        VECTOR.fill( elem );
      }

      /// Copy constructor
      Vector( const Vector<T>& source )
      {
        SIZE = source.SIZE;
        VECTOR = source.VECTOR;
      }

      /// Destructor
      ~Vector() { }

      /* ----- Operator overloading ----- */

      /// Output operator <<
			template <class Type>
			friend std::ostream& operator<<( std::ostream& os, const Vector<Type>& vec );

      /// Indexing operator ( read only )
			const T& operator[] ( const std::size_t& i ) const
			{
					// Range check 
					if ( i<0 || SIZE<=i )	{ throw Error( "Vector range error" );}
					return this->VECTOR( i, 0 );
			}

      /// Indexing operator ( read/write )
			T& operator[] ( const std::size_t& i )
			{
					// Range check 
					if ( i<0 || SIZE<=i )	{ throw Error( "Vector range error" );}
					return this->VECTOR( i, 0 );
			}

      /// Assignment
    	Vector& operator=( const Vector& original )
      {	
        SIZE = original.SIZE;	
		    VECTOR = original.VECTOR;
		    return *this;
	    }

      /// Unary +
      Vector<T> operator+() const
			{
				return *this;
			}  
      
      /// Unary -
      Vector<T> operator-() const
      {
        Vector<T> temp( *this );
		    temp.VECTOR = -temp.VECTOR;
		    return temp;
      }

      /// Binary +
      Vector<T> operator+( const Vector<T>& v_plus ) const
			{
				if ( v_plus.SIZE != SIZE ) { throw Error( "Vector dimension error" );}
				Vector<T> temp( *this );
		    temp.VECTOR += v_plus.VECTOR;
		    return temp;
			}

      /// Binary -
      Vector<T> operator-( const Vector<T>& v_minus ) const
			{
				if ( v_minus.SIZE != SIZE ) { throw Error( "Vector dimension error" );}
				Vector<T> temp( *this );
		    temp.VECTOR -= v_minus.VECTOR;
		    return temp;
			}

      /// Scalar multiplication (both sides)
      Vector<T> operator*( const T& v_times )const
      {
        Vector<T> temp( *this );
		    for (size_t i = 0; i < temp.SIZE; ++i)
		    {
				  temp.VECTOR( i, 0 ) *= v_times;
		    }
		    return temp;
      } 
      // Friend function so the * operator can be on either side
      friend Vector<T> operator*( const T& v, Vector<T>& vec )
			{
				return vec * v;
			}				

      /// Scalar division 
      Vector<T> operator/( const T& v_div )const
      {
        Vector<T> temp( *this );
		    for (size_t i = 0; i < temp.SIZE; ++i)
		    {
				  temp.VECTOR( i, 0 ) /= v_div;
		    }
		    return temp;
      } 
      
      /// Addition assignment  
      Vector<T> operator+=( const Vector<T>& v_plus )
      {
        if ( v_plus.SIZE != SIZE ) { throw Error( "Vector dimension error" );}
        VECTOR = VECTOR + v_plus.VECTOR;
        return *this;
      }    

      /// Subtraction assignment
      Vector<T> operator-=( const Vector<T>& v_minus )
      {
        if ( v_minus.SIZE != SIZE ) { throw Error( "Vector dimension error" );}
        VECTOR = VECTOR - v_minus.VECTOR;
        return *this;
      } 
    
      /* ----- Methods ----- */
      
      /// Return the size of the vector
      std::size_t size() const { return SIZE; }

      /// Push_back a new element into the end of the vector
      void push_back( const T& new_elem )
      {
        VECTOR.conservativeResize( SIZE + 1, 1 );
        VECTOR( SIZE, 0 ) = new_elem;
        ++SIZE;
      }

      /// Resize the vector
      void resize( const std::size_t size )
      {
        VECTOR.conservativeResize( size, 1 );
        SIZE = size;
      }

      /// Clear all the elements from the vector
      void clear()
      {
        VECTOR.resize( 0, 1 );
        SIZE = 0;
      }

      /// Return a vector of absolute values of the elements
      Vector<double> abs() const
      {
        Vector<double> abs_vals( SIZE );		
        for (size_t i=0; i < SIZE; ++i)
        {
          abs_vals.VECTOR( i, 0 ) = std::abs( VECTOR( i, 0 ) ) ;
        }
        return abs_vals; 
      }

      /// Swap elements i and j
      void swap( const std::size_t& i, const std::size_t& j )
      {
        if ( i<0 || SIZE<=i )	{ throw Error( "Vector range error" );}
        if ( j<0 || SIZE<=j )	{ throw Error( "Vector range error" );}
        if ( i == j ) { return; }
        std::swap<T>( VECTOR( i, 0 ), VECTOR( j, 0 ) );
      }

		  /// Assign new contents to the vector
      void assign( const std::size_t& size, const T& elem )
      {
        SIZE = size;
        VECTOR.resize( SIZE, 1 );
        VECTOR.fill( elem );
      }

      /// Create a linearly spaced vector (of doubles) with n elements
      void linspace( const double& a, const double& b, const std::size_t& n );

      /// Create a nonuniform vector using a power law with exponent p (p=1 ->linear)
      void power( const double& a, const double& b, const std::size_t& n, const double& p);

      /// Product of the elements in the vector (from index start to end)
      T product( const std::size_t& start, const std::size_t& end )
      {
        if ( start > end )	{ throw Error( "Vector product: start > end" );}
        if ( start<0 || SIZE<=start )	{ throw Error( "Vector range error" );}
        if ( end<0 || SIZE<=end )	{ throw Error( "Vector range error" );}
		    T prod( VECTOR( start, 0 ) ); // First value
		    for ( std::size_t i=start+1; i<=end; ++i )
		    {
			    prod *= VECTOR( i, 0 );
		    }
        return prod;
      }

      T product( )
      {
        return this->product( 0, SIZE - 1);
      }

      /// Sum of the elements in the vector (from index start to end)
      T sum( const std::size_t& start, const std::size_t& end )
      {
        if ( start > end )	{ throw Error( "Vector sum: start > end" );}
        if ( start<0 || SIZE<=start )	{ throw Error( "Vector range error" );}
        if ( end<0 || SIZE<=end )	{ throw Error( "Vector range error" );}
		    T S( VECTOR( start, 0 ) ); // First value
		    for ( std::size_t i=start+1; i<=end; ++i )
		    {
			    S += VECTOR( i, 0 );
		    }
        return S;
      }	

      T sum( )
      {
        return this->sum( 0, SIZE - 1);
      }

      /// Scale each element of the vector, equivalent to *=
      void scale( const T& m )
      {
        for( std::size_t i=0; i<SIZE; ++i)
        {
          VECTOR( i, 0 ) *= m;
        }
      }

      /// Return the dot product of two vectors v.dot(w)
      T dot( const Vector<T>& w )
      {
        if ( SIZE != w.SIZE )	{ throw Error( "Vector dot product: size error" );}
        T dp = VECTOR.adjoint()*w.VECTOR;
        return dp;
      }
      
      /* ----- Norms ----- */

      /// L1 norm: sum of absolute values 
      double norm_1() const
      {
        double sum( 0.0 );
        for (size_t i=0; i < SIZE; ++i)
        {
			    sum += std::abs( VECTOR( i, 0 ) );
		    }
		    return sum;  
	    }

      /// L2 norm: square root of the sum of the squares 
      double norm_2() const
      {
        double sum( 0.0 );
        for (size_t i=0; i < SIZE; ++i)
        {
			    sum += std::pow( std::abs( VECTOR( i, 0 ) ), 2.0 );
		    }
		    return std::sqrt( sum );  
	    }

      /// Lp norm: p-th root of the sum of the absolute values raised to the power p
      double norm_p( const double& p ) const
      {
        double sum( 0.0 );
        for (size_t i=0; i < SIZE; ++i)
        {
			    sum += std::pow( std::abs( VECTOR( i, 0 ) ), p );
		    }
		    return std::pow( sum , 1.0/p ); 
	    }

      /// Inf norm: largest absolute value element (p -> infinity)
      double norm_inf() const
	    {
		    std::vector<double> abs_vals;		
		    for (size_t i=0; i < SIZE; ++i)
		    {
			    abs_vals.push_back( std::abs( VECTOR( i, 0 ) ) );
		    }
		    return *std::max_element(abs_vals.begin(), abs_vals.end());
	    }


	}; // End of class Vector

  /// Output operator <<
  template <class Type>
  inline std::ostream& operator<<( std::ostream& os, const Vector<Type>& vec )
  {
    os << vec.VECTOR;
    return os;
  }  

  /// Create a linearly spaced vector (of doubles) with n elements
  template <>
	inline void Vector<double>::linspace( const double& a, const double& b,
                                        const std::size_t& n )
	{
		VECTOR.resize( n, 1 );
    SIZE = n;
		const double h = ( b - a ) / (n - 1)  ;
		for ( std::size_t i=0; i < n; ++i ) 
		{		
			VECTOR( i, 0 ) = a + h * i;
		}
	}

  /// Create a nonuniform vector using a power law with exponent p (p=1 ->linear)
  template <>
  inline void Vector<double>::power( const double& a, const double& b,
                                     const std::size_t& n, const double& p)
  {
    VECTOR.resize( n, 1 );
    SIZE = n;
    for ( std::size_t i = 0; i < n; ++i )
    {
      VECTOR( i , 0 ) = a + ( b - a ) * std::pow( ( double )i / ( n - 1 ), p );
    }
  }

} // End of namespace TSL

#endif
