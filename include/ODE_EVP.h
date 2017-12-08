/* ODE_EVP - Here we define the ODE_EVP class used for solving ODE linear
             eigenvalue problems.
*/

#ifndef ODE_EVP_H
#define ODE_EVP_H

#include <memory>
#include <algorithm>

#include "Vector.h"
#include "Matrix.h"
#include "Equation_2matrix.h"
#include "Error.h"
#include "OneD_node_mesh.h"
#include "Eigensystem.h"
#include "Residual.h"

namespace TSL
{

  template <typename T>
  class ODE_EVP
  {
    private:

      /// A method that constructs the banded matrix problem
      void assemble_dense_problem();

      /// The function associated with this instance.
      Equation_2matrix<T > *p_EQUATION;
      /// Pointer to the residual defining the LHS BC
      Residual<T > *p_LEFT_RESIDUAL;
      /// Pointer to the residual defining the RHS BC
      Residual<T > *p_RIGHT_RESIDUAL;
      /// The linear eigensystem
      Eigensystem<T> SYSTEM;
      /// Matrices
      Matrix<T> A_DENSE;
      Matrix<T> B_DENSE;
      /// the nodal distribution
      Vector<double> NODES;
      /// A vector of uniform meshes that store the eigenfunctions
      /// and eigenvalue (as the last dof at each nodal point)
      /// for use in later local refinement via the ODE_BVP class
      std::vector< OneD_node_mesh< std::complex<double> > > MESHES;
      /// has the eigensystem been constructed/solved?
      bool CONSTRUCTED;
      bool EIGENVALUES_COMPUTED;
      bool EIGENVECTORS_COMPUTED;

  public:

    /// Constructor
    ODE_EVP( Equation_2matrix<T > *equation_ptr,
             const Vector<double> &nodes,
             Residual<T>* ptr_to_left_residual,
             Residual<T>* ptr_to_right_residual );

    /// Destructor
    ~ODE_EVP();

    /// Formulate and solve the global eigenvalue problem
    void eigensolve( bool compute_evecs = false );

    /// Allow access to the underlying linear eigensystem
    Eigensystem<T>& eigensystem();

    /// Return the computed eigenvalues in a vector
    Vector< std::complex<double> > eigenvalues() const
    {
      if ( EIGENVALUES_COMPUTED ) { return SYSTEM.eigenvalues(); }
      else { throw Error( "Eigensystem: eigenvalues not computed." ); }
    }

    /// Return the matrix of eigenvectors ( each column is an eigenvector )
    Matrix< std::complex<double> > eigenvector_matrix() const
    {
      if ( EIGENVECTORS_COMPUTED ) { return SYSTEM.eigenvector_matrix(); }
      else { throw Error( "Eigensystem: eigenvectors not computed." ); }
    }

    /// Return an std::vector of eigenvectors
    std::vector< Vector< std::complex<double> > > eigenvectors() const
    {
      if ( EIGENVECTORS_COMPUTED ){ return SYSTEM.eigenvectors(); }
      else { throw Error( "Eigensystem: eigenvectors not computed." ); }
    }

    /*void add_tagged_to_mesh()
    {
      // clear the existing data if any
      MESHES.clear();
      // order of the equation
      unsigned order = p_EQUATION -> get_order();
      // get the eigenvalues
      Vector< std::complex<double> > vals( p_SYSTEM -> get_tagged_eigenvalues() );
      // get the eigenvectors
      Matrix< std::complex<double> > vecs( p_SYSTEM -> get_tagged_eigenvectors() );
      // loop through the eigenvectors
      for ( unsigned ivec = 0; ivec < vals.size(); ++ivec )
      {
        // make a mesh with the right node distribution
        // we'll increase the order by 1 to allow the eigenvalue to be
        // stored at each nodal point -- this is wasteful but very useful
        // for feeding this into a BVP for local refinement.
        OneD_node_mesh< std::complex<double> > eigfn( NODES, order + 1 );
        // loop through all nodes
        for ( unsigned node = 0; node < NODES.size(); ++node )
        {
          // complex vector of the dof at this node ( + 1 for the eigenvalue)
          Vector< std::complex<double> > vars_at_node( order + 1, 0.0 );
          // get the dof from the eigenvector
          for ( unsigned var = 0; var < order; ++var )
          {
            vars_at_node[ var ] = vecs[ ivec ][ node * order + var ];
          }
          // the last variable at each node is the corresponding eigenvalue
          vars_at_node[ order ] = vals[ ivec ];
          // set the first 'order' dof to be those from the eigenvector
          eigfn.set_nodes_vars( node, vars_at_node );
        }
        //// store the eigenvalue in the mesh at each node ... wasteful, but useful
        //// a complex vector filled with the same value N times
        //DenseVector<D_complex> c( NODES.size(), vals[ ivec ] );
        //// add it to the mesh -- for use in nonlinear local refinement via ODE_BVP
        //eigfn.push_var( c );
        // add the eigenfunction to the vector of meshes
        MESHES.push_back( eigfn );
      }
    }*/

    OneD_node_mesh< std::complex<double> > get_mesh( const unsigned& i ) const
    {
      if ( i > MESHES.size() )
      {
        std::string problem;
        problem = "You have tried to extract an eigenfunction from the ODE_EVP class\n";
        problem += "whose index is outside the range of stored meshes.\n";
        throw Error( problem );
      }
      return MESHES[ i ];
    }

  }; // End of class ODE_EVP
} // End of namespace TSL

#endif
