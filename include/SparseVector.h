/* A templated SparseVector class sparse, variable size, vector object.
*/

#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H

#include <map>

#include "Error.h"

namespace TSL
{

  /// An SparseVector class -- a sparse vector object.
  /// This is templated but intended ONLY for double or
  /// std::complex<double>.
  template <typename T>
  class SparseVector
  {
    typedef typename std::map< std::size_t, T >::const_iterator citer;
    typedef typename std::map< std::size_t, T >::iterator iter;

  public:

    /// Constructor for a non-filled vector, to be filled by the user.
    explicit SparseVector( const std::size_t& max );

    /// Copy constructor.
    SparseVector( const SparseVector& source );

    /// Assignment operator.
    SparseVector& operator=( const SparseVector& source );

    /// Return an iterator pointing to the beginning of the STL map storage
    iter begin();

    /// Return a constant iterator pointing to the beginning of the STL map storage
    citer begin() const;

    /// Return an iterator pointing to the end of the STL map storage
    iter end();

    /// Return a constant iterator pointing to the end of the STL map storage
    citer end() const;

    /// Resize the maximum length of the sparse vector.
    void resize( const std::size_t& length );

    /// Remove all elements from the sparse vector
    void clear();

    /// Erase an element from the vector
    void erase( const iter& pos );

    /// Erase an element from the vector
    void erase( const std::size_t& index );

    /// Swap elements i and j.
    void swap( const std::size_t& i, const std::size_t& j );

    /// Swap ALL elements with those of another SparseVector
    void swap( SparseVector<T>& X );

    /// Find the (max) size of the vector.
    std::size_t size() const;

    /// Find the number of non-zero elements in the vector.
    std::size_t nelts() const;

    /// Set an element of the vector
    T& set( const std::size_t& i );

    /// Get an element of the vector (zero if the element has not been set)
    const T& get( const std::size_t& i ) const;

    /// Equivalent to the 'get' method
    const T& operator[] ( const std::size_t& i ) const;

    /// Equivalent to the 'set' method
    T& operator[] ( const std::size_t& i );

    /// Operator overloading for sparse vector addition
    SparseVector<T> operator+( const SparseVector<T>& X ) const;

    /// Overloading for +
    SparseVector<T> operator+() const;

    /// Operator overloading for sparse vector subtraction
    SparseVector<T> operator-( const SparseVector<T>& X ) const;

    /// Overloading for -
    SparseVector<T> operator-() const;

    /// Overloading multiplication for a scalar
    SparseVector<T> operator*( const T& m ) const;

    /// Overloading *= for scalar multiplication
    SparseVector<T>& operator*=( const T& m );

    /// Overloading -= for sparse vectors
    SparseVector<T>& operator-=( const SparseVector<T>& X );

    /// Overloading += for sparse vectors
    SparseVector<T>& operator+=( const SparseVector<T>& X );

    /// Look for an element index
    citer find( std::size_t i ) const
    {
      citer pos;
      pos = VEC.find( i );
      return pos;
    }

    /// A dot product.
    const T dot( const SparseVector<T>& x ) const;

    /// l1-norm
    double one_norm() const;

    /// l2-norm.
    double two_norm() const;

    /// Infinity norm.
    double inf_norm() const;

    /// Scale each element of the vector.
    void scale( const T& scale );

    /// Add a vector
    void add( const SparseVector<T>& X );

    /// Subtract a vector, element wise.
    void sub( const SparseVector<T>& X );

    /// Find the index of the element nearest in value to that specified
    std::size_t nearest_index( const T& value ) const;

    /// Find the index of the maximum element in the vector
    std::size_t maxabs_index() const;

    /// Find the index of the maximum element in the vector
    std::size_t minabs_index() const;

    /// Output the sparse vector's contents.
    void dump() const;

  private:
    // private storage of the data - encapsulated in std::vector
    std::map< std::size_t, T > VEC;
    T ZERO;
    std::size_t MAX_SIZE;

  }
  ; // end class


  // INLINE METHODS BELOW

  template <typename T>
  inline typename std::map< std::size_t, T >::iterator SparseVector<T>::begin()
  {
    return VEC.begin();
  }

  template <typename T>
  inline typename std::map< std::size_t, T >::const_iterator SparseVector<T>::begin() const
  {
    return VEC.begin();
  }

  template <typename T>
  inline typename std::map< std::size_t, T >::iterator SparseVector<T>::end()
  {
    return VEC.end();
  }

  template <typename T>
  inline typename std::map< std::size_t, T >::const_iterator SparseVector<T>::end() const
  {
    return VEC.end();
  }

  template <typename T>
  inline T& SparseVector<T>::operator[] ( const std::size_t& i )
  {

    if ( i >= size() )
    {
      std::string problem;
      problem = " The SparseVector.operator[] method is trying to access \n";
      problem += " outside the container. \n";
      throw Error( problem );
    }

    return VEC[ i ];
  }

  template <typename T>
  inline const T& SparseVector<T>::operator[] ( const std::size_t& i ) const
  {

    if ( i >= size() )
    {
      std::string problem;
      problem = " The SparseVector.operator[] method is trying to access \n";
      problem += " outside the container. \n";
      throw Error( problem );
    }

    return get( i );
  }

  template <typename T>
  inline T& SparseVector<T>::set( const std::size_t& i )
  {

    if ( i >= size() )
    {
      std::string problem;
      problem = " The SparseVector.set method is trying to access \n";
      problem += " outside the container. \n";
      throw Error( problem );
    }

    return VEC[ i ];
  }

  template <typename T>
  inline const T& SparseVector<T>::get( const std::size_t& i ) const
  {

    if ( i >= size() )
    {
      std::string problem;
      problem = " The SparseVector.get method is trying to access \n";
      problem += " outside the container. \n";
      throw Error( problem );
    }

    citer pos;
    pos = VEC.find( i );
    if ( pos != VEC.end() )
    {
      return pos -> second;
    }
    else
    {
      return ZERO;
    }
  }

  template <typename T>
  inline std::size_t SparseVector<T>::size() const
  {
    return MAX_SIZE;
  }

  template <typename T>
  inline std::size_t SparseVector<T>::nelts() const
  {
    return VEC.size();
  }

  template <typename T>
  inline void SparseVector<T>::erase( const std::size_t& index )
  {
    VEC.erase( index );
  }

  template <typename T>
  inline void SparseVector<T>::erase( const iter& pos )
  {
    VEC.erase( pos );
  }

  template <typename T>
  inline void SparseVector<T>::swap( SparseVector<T>& X )
  {
#ifdef PARANOID
    if ( X.size() != size() )
    {
      std::string problem;
      problem = " The SparseVector.swap method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw Error( problem );
    }
#endif
    VEC.swap( X.VEC );
  }

  template <typename T>
  inline SparseVector<T> SparseVector<T>::operator+() const
  {
    return * this;
  }


} // end namespace

#endif
