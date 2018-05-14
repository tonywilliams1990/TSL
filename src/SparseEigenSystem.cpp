#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>

#include "SparseEigenSystem.h"
#include "Error.h"
#include "Timer.h"

#ifdef SLEPC

namespace TSL
{

  template <typename T>
  SparseEigenSystem<T>::SparseEigenSystem( SparseMatrix<T >* Aptr, SparseMatrix<T >* Bptr )
  {
    REGION_DEFINED = false;
    GUESS_DEFINED = false;
    //
    p_A = Aptr;
    p_B = Bptr;
    NEV = 8; NCONV = 0;
    ORDER = (EPSWhich)7; // target magnitude is the default
    //CALC_EIGENVECTORS = true; // SLEPc methods *always* obtain eigenvecs ???
    CALC_EIGENVECTORS = false;
  }

  template <typename T>
  SparseEigenSystem<T>::~SparseEigenSystem()
  {
  }

  template <typename T>
  unsigned SparseEigenSystem<T>::get_nconv() const
  {
    return NCONV;
  }

  template <typename T>
  void SparseEigenSystem<T>::set_nev( unsigned n )
  {
    NEV = n;
  }

  template <typename T>
  void SparseEigenSystem<T>::set_target( std::complex<double> target )
  {
    // defaults to (0,0)
    SHIFT = target;
  }

  template <typename T>
  void SparseEigenSystem<T>::set_order( std::string order_string )
  {
    int flag(0);
    if ( order_string == "EPS_LARGEST_MAGNITUDE" ) { ORDER=(EPSWhich)1; flag=1; }
    if ( order_string == "EPS_SMALLEST_MAGNITUDE" ) { ORDER=(EPSWhich)2; flag=1; }
    if ( order_string == "EPS_LARGEST_REAL" ) { ORDER=(EPSWhich)3; flag=1; }
    if ( order_string == "EPS_SMALLEST_REAL" ) { ORDER=(EPSWhich)4; flag=1; }
    if ( order_string == "EPS_LARGEST_IMAGINARY" ) { ORDER=(EPSWhich)5; flag=1; }
    if ( order_string == "EPS_SMALLEST_IMAGINARY" ) { ORDER=(EPSWhich)6; flag=1; }
    if ( order_string == "EPS_TARGET_MAGNITUDE" ) { ORDER=(EPSWhich)7; flag=1; }
    if ( order_string == "EPS_TARGET_REAL" ) { ORDER=(EPSWhich)8; flag=1; }
    if ( order_string == "EPS_TARGET_IMAGINARY" ) { ORDER=(EPSWhich)9; flag=1; }
    if ( order_string == "EPS_ALL" ) { ORDER=(EPSWhich)10; flag=1; }
    if (flag==0)
    {
      std::string problem;
      problem = "The SparseEigenSystem::set_order method has been called\n";
      problem += "with an ordering_string that is not recognised.\n";
      throw Error( problem );
    }
  }

  template <typename T>
  bool& SparseEigenSystem<T>::region_defined()
  {
    return REGION_DEFINED;
  }

  template <typename T>
  bool& SparseEigenSystem<T>::calc_eigenvectors()
  {
    return CALC_EIGENVECTORS;
  }

  template <typename T>
  void SparseEigenSystem<T>::set_region( const double& left, const double& right, const double& bottom, const double& top )
  {
    REGION_DEFINED = true;
    REAL_L = left;
    REAL_R = right;
    IMAG_B = bottom;
    IMAG_T = top;
  }

  template <typename T>
  bool& SparseEigenSystem<T>::guess_defined()
  {
    return GUESS_DEFINED;
  }

  template <typename T>
  void SparseEigenSystem<T>::set_initial_guess( const Vector<T>& guess )
  {
    GUESS_DEFINED = true;
    INITIAL_GUESS = guess;
  }

  template <typename T>
  void SparseEigenSystem<T >::eigensolve()
  {
    // only one method available: SLEPc's generalized shift-invert solver
    eigensolve_slepc();
  }


  template <typename T>
  void SparseEigenSystem<T>::eigensolve_slepc()
  {
#ifndef SLEPC
    std::string problem;
    problem = "The SparseEigenSystem::eigensolve method has been called\n";
    problem += "but the compiler option -DSLEPC was not provided when\n";
    problem += "the library was built.";
    throw Error( problem );
#else
    // the SLEPc and PETSc data types needed
    Mat A,B;
    PetscInt n;
#ifdef PETSC_Z
    // for complex PETSC we only need one (complex) vector to store eigenvec.
    Vec x;
#endif
#ifdef PETSC_D
    // for double PETSC we need separate real and imag parts for eigenvec.
    Vec xr,xi;
#endif
    ST st;
    EPS eps;


    // assuming A & B are square
    n = p_A -> rows();

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Define the matrices that define the eigensystem, Ax=lambdaBx
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    // we need to convert from the native sparse format to that required by SLEPc/PETSc

    std::cout << "  * Starting matrix assembly\n";
    //Timer timer;
    //timer.start();

    // define A using PETSc structures
    //MatCreate(p_LIBRARY -> get_Comm(),&A);
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,n,n,n,n);
    MatSetFromOptions(A);
    // we ABSOLUTELY MUST pre-allocate, otherwise the performance really is AWFUL!
    // get the number of non-zero elts in each row as a vector
    PetscInt* all_rows_nnz = new PetscInt[ n ];
    p_A -> nelts_all_rows( all_rows_nnz );
    // allocate memory using the number of non-zero elts in each row (the 0 is ignored here)
    MatSeqAIJSetPreallocation(A, 0, all_rows_nnz );
    MatSetUp(A);
    for ( PetscInt i = 0; i<n; ++i )
    {
      // move the matrix data into PETSc format 1 row at a time
      std::size_t nelts_in_row = p_A -> numel_row(i);
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
    // delete the temp storage
    delete[] all_rows_nnz;

    // MatSetValue inserted values are generally cached
    // so we need to explicitly do final assembly
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    std::cout << "A assembled \n";
    //timer.print();
    //timer.stop();
    //timer.print();
    //timer.reset();
    //timer.start();

    // define B using PETSc structures
    //MatCreate(p_LIBRARY -> get_Comm(),&B);
    MatCreate(PETSC_COMM_WORLD,&B);
    MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n);
    MatSetFromOptions(B);
    // we ABSOLUTELY MUST pre-allocate, otherwise the performance really is AWFUL!
    // get the number of non-zero elts in each row as a vector
    all_rows_nnz = new PetscInt[ n ];
    p_B -> nelts_all_rows( all_rows_nnz );
    // allocate memory using the number of non-zero elts in each row (the 0 is ignored here)
    MatSeqAIJSetPreallocation(B, 0, all_rows_nnz );
    MatSetUp(B);

    for ( PetscInt i = 0; i<n; ++i )
    {
      // move the matrix data into PETSc format 1 row at a time
      std::size_t nelts_in_row = p_B -> numel_row(i);
      // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
      PetscInt* cols = new PetscInt[all_rows_nnz[i]];
      // store the non-zero elts in this row
      PetscScalar* storage = new PetscScalar[all_rows_nnz[i]];
      // get the data from the CppNoddy sparse matrix structure
      p_B -> get_row_petsc( i, storage, cols );
      MatSetValues(B,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
      // delete temp storage made in the conversion
      delete[] cols; delete[] storage;
    }
    // delete the temp storage
    delete[] all_rows_nnz;
    // MatSetValue inserted values are generally cached
    // so we need to explicitly do final assembly
    MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

    std::cout << "B assembled \n";

    //timer.print();
    //timer.stop();

    // PETSc storage for the eigenvector, using A to define the size
#ifdef PETSC_D
    MatCreateVecs(A,NULL,&xr);
    MatCreateVecs(A,NULL,&xi);
#endif
#ifdef PETSC_Z
    MatCreateVecs(A,NULL,&x);
#endif

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  Create the eigensolver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    // create the eigensolver environment -- the MPI_COMM_WORLD is stored in the
    // singleton instance to avoid creating/deleting it, which causes issues.
    //EPSCreate(p_LIBRARY -> get_Comm(),&eps);
    EPSCreate(PETSC_COMM_WORLD,&eps);
    // define a generalised problem with a B -- default case is non-Hermitian
    EPSSetOperators(eps,A,B);

    // Default is generalizzed non-Hermitian
    EPSSetProblemType(eps,EPS_GNHEP);
    // Method is Krylov Schur
    EPSSetType(eps, EPSKRYLOVSCHUR);
    // target spectrum shift - defaults to (0,0)
#ifdef PETSC_D
    EPSSetTarget(eps, SHIFT.real());
#endif
#ifdef PETSC_Z
    EPSSetTarget(eps, SHIFT);
#endif
    // set the order of the returned ev's, as set by the get_order method.
    EPSSetWhichEigenpairs(eps, ORDER);
    // set the number of requested ev's. Not sure if this is relevant if REGION_DEFINED
    EPSSetDimensions(eps,NEV,2*NEV+10,PETSC_DEFAULT);

    /*
       Define the region containing the eigenvalues of interest
    */
    if ( REGION_DEFINED )
    {
      RG rg;
      EPSGetRG(eps, &rg);
      RGSetType(rg, RGINTERVAL);
      RGIntervalSetEndpoints(rg,REAL_L,REAL_R,IMAG_B,IMAG_T);
      // it is possible to "invert" (take the complement of) the region
      //RGSetComplement(rg,PETSC_TRUE);
    }

    // get access to the spectral transformation
    EPSGetST(eps, &st);
    // we have to use "STSINVERT" instead of the default, because B is
    // typically singular for all problems I'm interested in.
    STSetType(st, STSINVERT);

    /*
       Defaults to using the MUMPS solver
    */
    // KSP is the linear solver object of the PETSc library
    KSP ksp;
    STGetKSP(st, &ksp);
    // set to precondition only
    KSPSetType(ksp, KSPPREONLY);
    // get a preconditioner object
    PC pc;
    KSPGetPC(ksp,&pc);
    // set it to LU factorization is precondition
    PCSetType(pc,PCLU);
    // solve using the SUPERLU_DIST library
    //PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);
    //PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
    PCFactorSetMatSolverType(pc,MATSOLVERMUMPS); // Updated in PETSC 3.9.1


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve the eigensystem
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    EPSSolve(eps);
    // EPSGetDimensions(eps,&nev,NULL,NULL);
    // update the NEV private data with the number of returned eigenvalues
    // NEV = (unsigned)nev; // is this always the same as the input nev?

    // Optional: Get some information from the solver and display it
/*  // DEBUG
    PetscInt its, lits, maxit;
    EPSType type;
    PetscReal tol;
    //KSP ksp;
    //
    std::cout << "[DEBUG] Target location for eigenvalue  = " << SHIFT << "\n";
    std::cout << "[DEBUG] Target ordering of returned eigenvalues (see EPSWhich enum) = " << ORDER << "\n";
    //
    EPSGetIterationNumber(eps,&its);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Number of iterations of the method: %D\n",its);
    EPSGetST(eps,&st); STGetKSP(st,&ksp);
    KSPGetTotalIterations(ksp,&lits);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Number of linear iterations of the method: %D\n",lits);
    EPSGetType(eps,&type);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Solution method: %s\n\n",type);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Number of requested eigenvalues: %D\n",NEV);
    EPSGetTolerances(eps,&tol,&maxit);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);
*/
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
    EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);
    EPSErrorView(eps,EPS_ERROR_ABSOLUTE,PETSC_VIEWER_STDOUT_WORLD);
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
//#endif

    /*
       Save eigenvectors, if requested
    */
    PetscInt nconv;
    EPSGetConverged(eps,&nconv);
    // store it in the class
    NCONV = (unsigned)nconv;
    // create a complex eigenvalue vector
    ALL_EIGENVALUES = Vector< std::complex<double> >( NCONV, 0.0 );
    // complex eigenvector matrix
    ALL_EIGENVECTORS = Matrix< std::complex<double> >( NCONV, n, 0.0 );

    //
    for ( unsigned i=0; i<NCONV; i++ )
    {
#ifdef PETSC_Z
      EPSGetEigenvalue(eps,i,&ALL_EIGENVALUES[i],NULL);
      //std::cout << ALL_EIGENVALUES[i] << "\n";
#endif
#ifdef PETSC_D
      double lambda_r,lambda_i;
      EPSGetEigenvalue(eps,i,&lambda_r,&lambda_i);
      ALL_EIGENVALUES[i] = std::complex<double>( lambda_r, lambda_i );
#endif

      if ( CALC_EIGENVECTORS )
      {
#ifdef PETSC_D
        // get the i-th eigenvector from SLEPc
        EPSGetEigenvector(eps,i,xr,xi);
        // convert to a more accessible data structure
        PetscScalar* arrayr;
        VecGetArray1d( xr, n, 0, &arrayr );
        PetscScalar* arrayi;
        VecGetArray1d( xi, n, 0, &arrayi );
        for ( int j=0; j<n; ++j )
        {
          //ALL_EIGENVECTORS[i][j]=std::complex<double>( arrayr[j], arrayi[j] );
          ALL_EIGENVECTORS( i, j ) = std::complex<double>( arrayr[j], arrayi[j] );
        }
        // documentation says to "restore", though it might not matter as we're done with it now
        VecRestoreArray1d( xr, n, 0, &arrayr );
        VecRestoreArray1d( xi, n, 0, &arrayi );
#endif
#ifdef PETSC_Z
        // get the i-th eigenvector from SLEPc
        EPSGetEigenvector(eps,i,x,NULL);
        // convert to a more accessible data structure
        PetscScalar* array;
        VecGetArray1d( x, n, 0, &array );
        for ( int j=0; j<n; ++j )
        {
          //ALL_EIGENVECTORS[i][j]=array[j];
          ALL_EIGENVECTORS( i, j ) = array[j];
        }

        // documentation says to "restore", though it might not matter as we're done with it now
        VecRestoreArray1d( x, n, 0, &array );
#endif
      }
    }
    /*
       Free work space
    */

    EPSDestroy(&eps);
    MatDestroy(&A); MatDestroy(&B);
#ifdef PETSC_D
    VecDestroy(&xr); VecDestroy(&xi);
#endif
#ifdef PETSC_Z
    VecDestroy(&x);
#endif

#endif
  }


#ifdef PETSC_Z
  template class SparseEigenSystem< std::complex<double> >;
#endif
#ifdef PETSC_D
  template class SparseEigenSystem<double>;
#endif


} // end namespace

#endif
