#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>

#include "PETScSparseLinearSystem.h"
#include "Error.h"
#include "Timer.h"

#if defined(PETSC_D) || defined(PETSC_Z)
   #include "petscksp.h"
   #include "petsc.h"
   #include "mpi.h"
#endif

namespace TSL
{

  template <typename T>
  PETScSparseLinearSystem<T>::PETScSparseLinearSystem( SparseMatrix<T >* Aptr, Vector<T >* Bptr )
  {
    p_A = Aptr;
    p_B = Bptr;
    factorised_ = false;

#if defined(PETSC_D) || defined(PETSC_Z)
    int flag(0);
    MPI_Initialized( &flag );
    if ( flag != 1 )
    {
      std::string problem;
      problem = "The PETScSparseLinearSystem has been instantiated for a petsc solver.\n";
      problem += "You must call PetscInitialize before calling the petsc solver.\n";
      throw Error( problem );
    }
    MPI_Comm_size(MPI_COMM_WORLD,&size_);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
    if ( size_ > 1 )
    {
       std::string problem;
       problem = " The SparseLinearSystem object links to PETSc which makes\n";
       problem += " use of MPI, but you probably won't gain anything yet by\n";
       problem += " using >1 processes. The SparseMatrix object is still too\n";
       problem += " dumb, and will be stored in each process.\n";
       throw Error( problem );
    }

#endif

  }

  template<typename T>
   PETScSparseLinearSystem<T>::~PETScSparseLinearSystem()
   {
     cleanup();
   }

   template<typename T>
   void PETScSparseLinearSystem<T>::cleanup()
   {
     std::cout << "in cleanup.\n";
     #if defined(PETSC_D) || defined(PETSC_Z)
       // delete objects used in the factorisation
       if (factorised_)
       {
         VecDestroy(&x_);
         VecDestroy(&B_);
         KSPDestroy(&ksp_);
         factorised_ = false;
       }
     #endif
   }

  template <typename T>
  void PETScSparseLinearSystem<T>::solve()
  {
    factorise();
    solve_using_factorisation();
  }

  template<>
   void PETScSparseLinearSystem<double>::factorise()
   {
   #if !defined(PETSC_D)
     std::string problem;
     problem = "CppNoddy is linked against the COMPLEX version of PETSc\n";
     problem += "but you are trying to factorise a DOUBLE matrix. Either\n";
     problem += "redefine your matrix as complex, or recompile with $PETSC_ARCH\n";
     problem += "pointing to a DOUBLE version of the PETSc code.";
     throw Error( problem );
   #endif
   #if defined(PETSC_D)
     if (factorised_)
     {
       // already factorised -- so delete and re-create below
       cleanup();
     }

     // store a boolean to indicate that we
     factorised_ = true;
     PetscInt Istart,Iend,n;
     Mat A;

     // size of the (assumed square) matrix
     n = p_A -> rows();

     // configure the A matrix
     MatCreate(PETSC_COMM_WORLD,&A);
     MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
     MatSetFromOptions(A);

     // get: all_rows_nnz[i] is the number of nonzero elts in row i
     PetscInt* all_rows_nnz = new PetscInt[ n ];
     p_A -> nelts_all_rows( all_rows_nnz );

     // pre-allocate memory using the number of non-zero elts in each row (the 0 is ignored here)
     MatSeqAIJSetPreallocation(A, 0, all_rows_nnz );
     // need to allocate for MPI codes too
     // \todo if we every get MPI running, we need to sort out preallocation
     // MatMPIAIJSetPreallocation(A, 800, NULL, 800, NULL);
     //
     // finish the A definition
     MatSetUp(A);

     /*
        Currently, all PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors.  Determine which
        rows of the matrix are locally owned.
     */
     MatGetOwnershipRange(A,&Istart,&Iend);
     // populate the A matrix using the CppNoddy sparse matrix data
     for ( PetscInt i = Istart; i<Iend; ++i )
     {
       // move the matrix data into PETSc format 1 row at a time
       std::size_t nelts_in_row = all_rows_nnz[i]; //p_A -> nelts_in_row(i);
       // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
       PetscInt* cols = new PetscInt[all_rows_nnz[i]];
       // store the non-zero elts in this row
       PetscScalar* storage = new PetscScalar[all_rows_nnz[i]];
       // get the data from the CppNoddy sparse matrix structure
       p_A -> get_row_petsc( i, storage, cols );
       MatSetValues(A,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
       // delete temp storage made in the conversion
       delete[] cols; delete[] storage;
     }

     /*
        Assemble matrix, using the 2-step process:
          MatAssemblyBegin(), MatAssemblyEnd()
        Computations can be done while messages are in transition
        by placing code between these two statements.
     */
     MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
     MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

     /*
        Create parallel vectors.
     */
     VecCreate(PETSC_COMM_WORLD,&B_);
     VecSetSizes(B_,PETSC_DECIDE,p_A->rows());
     VecSetFromOptions(B_);
     VecDuplicate(B_,&x_);

     PetscInt low, high;
     VecGetOwnershipRange(x,&low,&high);
     VecAssemblyBegin(B_);
     VecAssemblyEnd(B_);

     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   Create the linear solver and set various options
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
     KSPCreate(PETSC_COMM_WORLD,&ksp_);
     KSPSetOperators(ksp_,A,A);
     KSPSetType(ksp_,KSPPREONLY);
     PetscInt  ival,icntl;
     PetscReal val;
     KSPGetPC(ksp_,&pc_);
     // hardwire a DIRECT SOLVER via MUMPS
     PCSetType(pc_,PCLU);
     //PCFactorSetMatSolverPackage(pc_,MATSOLVERMUMPS);
     //PCFactorSetUpMatSolverPackage(pc_);
     PCFactorSetMatSolverType(pc_,MATSOLVERMUMPS); // Updated in PETSC 3.9.1
     /* call MatGetFactor() to create F */
     //PCFactorGetMatrix(pc_,&F_);

     /* sequential ordering */
     //icntl = 7; ival = 2;
     //MatMumpsSetIcntl(F_,icntl,ival);

     /* threshhold for row pivot detection */
     //MatMumpsSetIcntl(F_,24,1);
     //icntl = 3; val = 1.e-6;
     //MatMumpsSetCntl(F_,icntl,val);

     /* Get info from matrix factors */
     KSPSetUp(ksp_);
     PCFactorGetMatrix(pc_,&F_);

     /* sequential ordering */
     icntl = 7; ival = 2;
     MatMumpsSetIcntl(F_,icntl,ival);

     /* threshhold for row pivot detection */
     MatMumpsSetIcntl(F_,24,1);
     icntl = 3; val = 1.e-6;
     MatMumpsSetCntl(F_,icntl,val);

     MatDestroy(&A);
     delete[] all_rows_nnz;
   #endif // check for PETSC_D/Z
   }

   template <>
   void PETScSparseLinearSystem<double>::solve_using_factorisation()
   {
   #if !defined(PETSC_D)
     std::string problem;
     problem = "CppNoddy is linked against the COMPLEX version of PETSc\n";
     problem += "but you are trying to solve a DOUBLE matrix. Either\n";
     problem += "redefine your matrix as complex, or recompile with $PETSC_ARCH\n";
     problem += "pointing to a DOUBLE version of the PETSc code.";
     throw Error( problem );
   #endif
   #if defined(PETSC_D)
     // size of the (assumed square) matrix
     /*PetscInt n = p_A -> rows();

     // populate the RHS vector using the CppNoddy DenseVector content
     for ( PetscInt i = 0; i < n; ++i )
     {
       VecSetValue(B_,i,p_B->operator[](i),INSERT_VALUES);
     }

     KSPSolve(ksp_,B_,x_);

     // We can now gather the parallel result back to ALL processes
    //   This is temporary as the SparseMatrix is stored on each processes
    //   and is too dumb for "proper" parallelization
     Vec y;
     // a scatter context
     VecScatter ctx = 0;
     // map all elts of the parallel vector to a sequential copy
     VecScatterCreateToAll(x_,&ctx,&y);
     // scatter it
     VecScatterBegin(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
     VecScatterEnd(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
     // clean up
     VecScatterDestroy(&ctx);
     // this array is a pointer not a copy
     PetscScalar* array;
     VecGetArray(y,&array);
     // now copy to the CppNoddy densevctor
     for (PetscInt i=0; i<n; i++)
     {
       p_B -> operator[](i) = array[i];
     }
     // follow the docs and Restore after get
     VecRestoreArray(x_,&array);
     VecDestroy(&y);*/
   #endif
   }



  template<>
  void PETScSparseLinearSystem<std::complex<double> >::factorise()
  {
  #if !defined(PETSC_Z)
    std::string problem;
    problem = "TSL is linked against the DOUBLE version of PETSc\n";
    problem += "but you are trying to factorise a COMPLEX matrix.\n";
    problem += "Recompile with $PETSC_ARCH\n";
    problem += "pointing to a COMPLEX version of the PETSc code.";
    throw Error( problem );
  #endif

  #if defined(PETSC_Z)
    if (factorised_)
    {
      // already factorised -- so delete and re-create below
      cleanup();
    }

    // store a boolean to indicate that we
    factorised_ = true;
    PetscInt Istart,Iend,n;
    Mat A;

    // size of the (assumed square) matrix
    n = p_A -> rows();
    /*
       Create parallel vectors.
    */
    VecCreate( PETSC_COMM_WORLD, &B_ );
    VecSetSizes( B_, PETSC_DECIDE, p_A->rows() );
    VecSetFromOptions( B_ );
    VecDuplicate( B_, &x_ );

    // configure the A matrix
    MatCreate(PETSC_COMM_WORLD,&A);
    // set A to be an nxn matrix
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    MatSetFromOptions(A);

    // get: all_rows_nnz[i] is the number of nonzero elts in row i
    PetscInt* all_rows_nnz = new PetscInt[ n ];
    p_A -> nelts_all_rows( all_rows_nnz );

    // pre-allocate memory using the number of non-zero elts
    // in each row (the 0 is ignored here)
    MatSeqAIJSetPreallocation(A, 0, all_rows_nnz );
    // need to allocate for MPI codes too
    // \todo if we every get MPI running, we need to sort out preallocation
    // MatMPIAIJSetPreallocation(A, 800, NULL, 800, NULL);
    //
    // finish the A definition
    MatSetUp(A);

    /*
       Currently, all PETSc parallel matrix formats are partitioned by
       contiguous chunks of rows across the processors.  Determine which
       rows of the matrix are locally owned.
    */
    MatGetOwnershipRange(A,&Istart,&Iend);
    // populate the A matrix using the CppNoddy sparse matrix data
    for ( PetscInt i = Istart; i<Iend; ++i )
    {
      // move the matrix data into PETSc format 1 row at a time
      std::size_t nelts_in_row = all_rows_nnz[i];
      // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
      PetscInt* cols = new PetscInt[nelts_in_row];
      // store the non-zero elts in this row
      PetscScalar* storage = new PetscScalar[nelts_in_row];
      // get the data from the CppNoddy sparse matrix structure
      p_A -> get_row_petsc( i, storage, cols );
      MatSetValues(A,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
      // delete temp storage made in the conversion
      delete[] cols; delete[] storage;
    }

    /*
       Assemble matrix, using the 2-step process:
         MatAssemblyBegin(), MatAssemblyEnd()
       Computations can be done while messages are in transition
       by placing code between these two statements.
    */
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  Create the linear solver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    */
    KSPCreate(PETSC_COMM_WORLD,&ksp_);
    KSPSetOperators(ksp_,A,A);
    KSPSetType(ksp_,KSPPREONLY);
    PetscInt  ival,icntl;
    PetscReal val;
    KSPGetPC(ksp_,&pc_);
    // hardwire a DIRECT SOLVER via MUMPS

    PCSetType(pc_,PCLU);
    //PCFactorSetMatSolverPackage(pc_,MATSOLVERMUMPS);
    PCFactorSetMatSolverType(pc_,MATSOLVERMUMPS); // Updated in PETSC 3.9.1
    //PCFactorSetUpMatSolverPackage(pc_);
    /* call MatGetFactor() to create F */
    //PCFactorGetMatrix(pc_,&F_); // moved after KSPSetUp

    /* sequential ordering */
    //icntl = 7; ival = 2;
    //MatMumpsSetIcntl(F_,icntl,ival);

    /* threshhold for row pivot detection */
    //MatMumpsSetIcntl(F_,24,1);
    //icntl = 3; val = 1.e-6;
    //MatMumpsSetCntl(F_,icntl,val);

    /* Get info from matrix factors */
    KSPSetUp(ksp_);


    PCFactorGetMatrix(pc_,&F_);

    /* sequential ordering */
    icntl = 7; ival = 2;
    MatMumpsSetIcntl(F_,icntl,ival);

    /* threshhold for row pivot detection */
    MatMumpsSetIcntl(F_,24,1);
    icntl = 3; val = 1.e-6;
    MatMumpsSetCntl(F_,icntl,val);

    MatDestroy(&A);
    delete[] all_rows_nnz;

  #endif // check for PETSC_D/Z
  }

   template <>
   void PETScSparseLinearSystem<std::complex<double> >::solve_using_factorisation()
   {
   #if !defined(PETSC_Z)
     std::string problem;
     problem = "CppNoddy is linked against the DOUBLE version of PETSc\n";
     problem += "but you are trying to solve e a COMPLEX matrix.\n";
     problem += "Recompile with $PETSC_ARCH\n";
     problem += "pointing to a COMPLEX version of the PETSc code.";
     throw Error( problem );
   #endif
   #if defined(PETSC_Z)
     // size of the (assumed square) matrix
     PetscInt n = p_A -> rows();

     // populate the RHS vector using the CppNoddy DenseVector content
     for ( PetscInt i = 0; i < n; ++i )
     {
       VecSetValue(B_,i,p_B->operator[](i),INSERT_VALUES);
     }

     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         Solve the linear system
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
     KSPSolve(ksp_,B_,x_);

     /* We can now gather the parallel result back to ALL processes
       This is temporary as the SparseMatrix is stored on each processes
       and is too dumb for "proper" parallelization */
     Vec y;
     // a scatter context
     VecScatter ctx = 0;
     // map all elts of the parallel vector to a sequential copy
     VecScatterCreateToAll(x_,&ctx,&y);
     // scatter it
     VecScatterBegin(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
     VecScatterEnd(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
     // clean up
     VecScatterDestroy(&ctx);
     // this array is a pointer not a copy
     PetscScalar* array;
     VecGetArray(y,&array);
     // now copy to the CppNoddy densevctor
     for (PetscInt i=0; i<n; i++)
     {
       p_B -> operator[](i) = array[i];
     }
     // follow the docs and Restore after get
     VecRestoreArray(x_,&array);
     VecDestroy(&y);
     #endif
   }

#ifdef PETSC_Z
  template class PETScSparseLinearSystem< std::complex<double> >;
#endif
#ifdef PETSC_D
  template class PETScSparseLinearSystem<double>;
#endif



} // end namespace
