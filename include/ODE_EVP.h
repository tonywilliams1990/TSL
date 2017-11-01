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

  template <typename _Type>
  class ODE_EVP
  {
    private:

      /// A method that constructs the banded matrix problem
      void assemble_dense_problem();

      /// The function associated with this instance.
      Equation_2matrix<_Type > *p_EQUATION;
      /// Pointer to the residual defining the LHS BC
      Residual<_Type > *p_LEFT_RESIDUAL;
      /// Pointer to the residual defining the RHS BC
      Residual<_Type > *p_RIGHT_RESIDUAL;
      /// The linear eigensystem
      Eigensystem<_Type> SYSTEM;
      /// Matrices
      Matrix<_Type> A_DENSE;
      Matrix<_Type> B_DENSE;
      /// the nodal distribution
      Vector<double> NODES;
      /// A vector of uniform meshes that store the eigenfunctions
      /// and eigenvalue (as the last dof at each nodal point)
      /// for use in later local refinement via the ODE_BVP class
      std::vector< OneD_node_mesh< std::complex<double> > > MESHES;
      /// has the eigensystem been constructed/solved?
      bool CONSTRUCTED;
      bool EIGENVALUES_COMPUTED;

  public:

    /// The class is defined by a vector function for the system.
    /// \param equation_ptr A pointer to an equation with 2 associated matrices; matrix1 will define the eigenvalue problem.
    /// \param nodes A vector of nodal points.
    /// \param ptr_to_left_residual A pointer to a residual object that defines the LHS boundary conditions.
    /// \param ptr_to_right_residual A pointer to a residual object that defines the RHS boundary conditions.
    ODE_EVP( Equation_2matrix<_Type > *equation_ptr,
             const Vector<double> &nodes,
             Residual<_Type>* ptr_to_left_residual,
             Residual<_Type>* ptr_to_right_residual );

    /// Destructor
    ~ODE_EVP();

    /// Formulate and solve the global eigenvalue problem
    /// for a linear system.
    void eigensolve();

    /// Allow access to the underlying dense linear eigensystem
    /// through a pointer to the private member data.
    Eigensystem<_Type>& eigensystem();

    /// Return the computed eigenvalues in a vector
    Vector< std::complex<double> > eigenvalues() const
    {
      if ( EIGENVALUES_COMPUTED ) { return SYSTEM.eigenvalues(); }
      else { throw Error( "Eigensystem: eigenvalues not computed." ); }
    }

    //TODO function for returning eigenvectors + decide if we want to compute evecs

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
#ifdef PARANOID
      if ( i > MESHES.size() )
      {
        std::string problem;
        problem = "You have tried to extract an eigenfunction from the ODE_EVP class\n";
        problem += "whose index is outside the range of stored meshes.\n";
        throw ExceptionRange( problem, MESHES.size(), i );
      }
#endif
      return MESHES[ i ];
    }

  }; // End of class ODE_EVP

} // End of namespace TSL

#endif
