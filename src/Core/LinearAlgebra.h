/* LinearAlgebra.h - Defines a class for various linear algebra methods.
   Created  		   - 06/08/2016
*/

#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include "Error.h"
#include "Matrix.h"
#include "Vector.h"

namespace TSL
{
	
	/// A templated linear algebra class
	template<class T> class LinearAlgebra
	{
		protected:
			

		public:
			
			/// Constructor 
			LinearAlgebra(){ }

			/// Destructor
			~LinearAlgebra(){ }

			/* ----- Methods ----- */

			// Calculate the permutation for LU decomposition
			std::vector<std::size_t> LU_permutation( Matrix<T>& A_mat )
			{
				if ( A_mat.rows() != A_mat.cols() ) { throw Error( "LU error: matrix not square." );}
				std::size_t n = A_mat.rows();
				std::vector<std::size_t> permutation( n );
				for ( std::size_t k=0; k < n; ++k )
				{
					permutation[k] = k;   // Initially set to be identity
				}
				
				for ( std::size_t j=0; j < n; ++j )
				{
					std::size_t max_i = j;
					for ( std::size_t i=j; i < n; ++ i )
					{
						if ( std::abs( A_mat(i,j) ) > std::abs( A_mat(max_i,j) ) )
						{
							max_i = i;
						}
					}
					if ( max_i != j ) // No need to swap if already in the right place
					{
						std::iter_swap( permutation.begin() + j, permutation.begin() + max_i );
					}
				}

				return permutation;
			}

			// LU decomposition PA=LU of a matrix ( returns a vector for the permutation )
			std::vector<std::size_t> LU_decompose( Matrix<T>& A_mat )
			{
				if ( A_mat.rows() != A_mat.cols() ) { throw Error( "LU error: matrix not square." );}
				std::size_t n = A_mat.rows();
				std::vector<std::size_t> permutation( n );
				permutation = this->LU_permutation( A_mat );

				// Create the permutation matrix
				TSL::Matrix<T> P;
				P.resize( n, n );
				for ( std::size_t k=0; k < n; ++k )
				{
					P(k,permutation[k]) = 1;
				}

				// Apply the permutation
				A_mat = P * A_mat;

				// Now decompose the permuted matrix 
				for ( std::size_t i=1; i<n; ++i ) 			// The first row remains unchanged
				{
					for ( std::size_t j=0; j<n; ++j )
					{
						if ( i > j ) // Lower entries
						{
							for ( std::size_t k=0; k<j; ++k )
							{
								A_mat(i,j) -= A_mat(k,j) * A_mat(i,k); 
							}
							A_mat(i,j) /=  A_mat(j,j); 		
						}
						if ( i <= j ) // Upper entries
						{
							for ( std::size_t k=0; k<i; ++k )
							{
								A_mat(i,j) -= A_mat(k,j) * A_mat(i,k); 
							}
						}		
					}
				}	
				return permutation;
			}

	};	// End of class LinearAlgebra

}  // End of namespace TSL

#endif
