#ifndef DUNE_FEM_FMATRIXEIGENVALUES_HH 
#define DUNE_FEM_FMATRIXEIGENVALUES_HH 

#include <iostream>
#include <cmath>
#include <cassert>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/fem/io/parameter.hh>

#if HAVE_ALGLIB
#include <alglib/evd.h>
#endif

#if HAVE_LAPACK 
// dsyev declaration (in liblapack)
extern "C" {

/*
    *
    **  purpose
    **  =======
    **
    **  xsyev computes all eigenvalues and, optionally, eigenvectors of a
    **  BASE DATA TYPE symmetric matrix a.
    **
    **  arguments
    **  =========
    **
    **  jobz    (input) char
    **          = 'n':  compute eigenvalues only;
    **          = 'v':  compute eigenvalues and eigenvectors.
    **
    **  uplo    (input) char
    **          = 'u':  upper triangle of a is stored;
    **          = 'l':  lower triangle of a is stored.
    **
    **  n       (input) long int
    **          the order of the matrix a.  n >= 0.
    **
    **  a       (input/output) BASE DATA TYPE array, dimension (lda, n)
    **          on entry, the symmetric matrix a.  if uplo = 'u', the
    **          leading n-by-n upper triangular part of a contains the
    **          upper triangular part of the matrix a.  if uplo = 'l',
    **          the leading n-by-n lower triangular part of a contains
    **          the lower triangular part of the matrix a.
    **          on exit, if jobz = 'v', then if info = 0, a contains the
    **          orthonormal eigenvectors of the matrix a.
    **          if jobz = 'n', then on exit the lower triangle (if uplo='l')
    **          or the upper triangle (if uplo='u') of a, including the
    **          diagonal, is destroyed.
    **
    **  lda     (input) long int
    **          the leading dimension of the array a.  lda >= max(1,n).
    **
    **  w       (output) BASE DATA TYPE array, dimension (n)
    **          if info = 0, the eigenvalues in ascending order.
    **
    **
    **
    **  info    (output) long int
    **          = 0:  successful exit
    **          < 0:  if info = -i, the i-th argument had an illegal value
    **          > 0:  if info = i, the algorithm failed to converge; i
    **                off-diagonal elements of an intermediate tridiagonal
    **                form did not converge to zero.
    **
**/
extern void dsyev_(const char* jobz, const char* uplo, const long
                   int* n, double* a, const long int* lda, double* w,
                   double* work, const long int* lwork, long int* info);
/*
    *
    **  purpose
    **  =======
    **
    **  xgeev computes for an n-by-n DATA TYPE nonsymmetric matrix a, the
    **  eigenvalues and, optionally, the left and/or right eigenvectors.
    **
    **  the right eigenvector v(j) of a satisfies
    **                   a * v(j) = lambda(j) * v(j)
    **  where lambda(j) is its eigenvalue.
    **  the left eigenvector u(j) of a satisfies
    **                u(j)**h * a = lambda(j) * u(j)**h
    **  where u(j)**h denotes the conjugate transpose of u(j).
    **
    **  the computed eigenvectors are normalized to have euclidean norm
    **  equal to 1 and largest component BASE DATA TYPE.
    **
    **  arguments
    **  =========
    **
    **  jobvl   (input) char
    **          = 'n': left eigenvectors of a are not computed;
    **          = 'v': left eigenvectors of are computed.
    **
    **  jobvr   (input) char
    **          = 'n': right eigenvectors of a are not computed;
    **          = 'v': right eigenvectors of a are computed.
    **
    **  n       (input) long int
    **          the order of the matrix a. n >= 0.
    **
    **  a       (input/output) DATA TYPE array, dimension (lda,n)
    **          on entry, the n-by-n matrix a.
    **          on exit, a has been overwritten.
    **
    **  lda     (input) long int
    **          the leading dimension of the array a.  lda >= max(1,n).
    *  WR      (output) DOUBLE PRECISION array, dimension (N)
    *  WI      (output) DOUBLE PRECISION array, dimension (N)
    *          WR and WI contain the real and imaginary parts,
    *          respectively, of the computed eigenvalues.  Complex
    *          conjugate pairs of eigenvalues appear consecutively
    *          with the eigenvalue having the positive imaginary part
    *          first.

    **
    **  vl      (output) DATA TYPE array, dimension (ldvl,n)
    **          if jobvl = 'v', the left eigenvectors u(j) are stored one
    **          after another in the columns of vl, in the same order
    **          as their eigenvalues.
    **          if jobvl = 'n', vl is not referenced.
    **          u(j) = vl(:,j), the j-th column of vl.
    **
    **  ldvl    (input) long int
    **          the leading dimension of the array vl.  ldvl >= 1; if
    **          jobvl = 'v', ldvl >= n.
    **
    **  vr      (output) DATA TYPE array, dimension (ldvr,n)
    **          if jobvr = 'v', the right eigenvectors v(j) are stored one
    **          after another in the columns of vr, in the same order
    **          as their eigenvalues.
    **          if jobvr = 'n', vr is not referenced.
    **          v(j) = vr(:,j), the j-th column of vr.
    **
    **  ldvr    (input) long int
    **          the leading dimension of the array vr.  ldvr >= 1; if
    **          jobvr = 'v', ldvr >= n.
    **
    **
    **
    **
    **  info    (output) long int
    **          = 0:  successful exit
    **          < 0:  if info = -i, the i-th argument had an illegal value.
    **          > 0:  if info = i, the qr algorithm failed to compute all the
    **                eigenvalues, and no eigenvectors have been computed;
    **                elements and i+1:n of w contain eigenvalues which have
    **                converged.
    **
**/
   extern void dgeev_(
        const char* jobvl,
        const char* jobvr,
        const long int* n,
        double* a,
        const long int* lda,
        double* wr,
        double* wi,
        const double* vl,
        const long int* ldvl,
        const double* vr,
        const long int* ldvr,
        double* work, const long int* lwork, 
        long int* info);
/////////////////////////////////////////////////////////////////////////////////
} // end extern C 
#endif

namespace Dune {

namespace Fem {

/** \brief calculates the eigen values of a field matrix 
    \param[in]  matrix matrix eigen values are calculated for 
    \param[out] eigenvalues FieldVector that contains eigen values in 
                ascending order 
*/
template <class MatrixType, class VectorType> 
static void eigenValues(const MatrixType& matrix,
                        VectorType& rew,
                        double &ev)  
{
  const long int N = matrix.rows();
  assert( matrix.rows() == matrix.cols() );
#if HAVE_ALGLIB
  {
    const unsigned int precision = 256;
    typedef amp::ampf<precision> Field;
    //typedef ap::template_1d_array<Field> VectorType;
    typedef ap::template_2d_array<Field> DenseMatrixType;
    VectorType lr,li;
    lr.setbounds( 0, N-1 );
    li.setbounds( 0, N-1 );
    DenseMatrixType a,vl,vr;
    a.setbounds( 0, N-1, 0, N-1 );
    vl.setbounds( 0, N-1, 0, N-1 );
    vr.setbounds( 0, N-1, 0, N-1 );
  
    // copy matrix  
    bool symetric = true ;
    for(int i=0; i<N; ++i) 
    {
      for(int j=0; j<N; ++j) 
      {
        // make sure matrix is symetric 
        if( std::abs( matrix( i , j ) -  matrix( j , i ) ) > 
            1e-12 * ( std::abs( matrix( i , j ) ) +  
                      std::abs( matrix( j , i ) ) )
          && std::abs( matrix( i , j ) -  matrix( j , i ) ) > 1e-12 
          )
        {
          std::cout << i << " " << j << " " << " : "
                    << matrix( i , j ) << " " << matrix( j , i ) << " -> "
                    << std::abs( matrix( i , j ) -  matrix( j , i ) )
                    << std::endl;
          symetric = false ; 
        }
        a(i,j) = matrix( i , j );
      }
    }

    
    evd::rmatrixevd(a,N,0,lr,li,vl,vr);

    std::cout << "Eigenvalues: " << std::endl;
    bool negative = false ; 
    bool imagEV = false;
    double minE=1e10;
    double maxE=-1e10;
    for(int i=0; i<N; ++i) 
    {
      if( lr( i ) < -1e-20 ) negative = true ;
      if ( std::abs( li(i).toDouble() ) > 1e-20 ) imagEV = true;
      minE = std::min( minE, lr( i ).toDouble() );
      maxE = std::max( maxE, lr( i ).toDouble() );
      if (Parameter::verbose())
        std::cout << std::scientific << lr( i ).toDouble() << " " 
                  << li(i).toDouble() << "i" << std::endl;
    }
    std::cout << "Symmetric: " << symetric << " "
              << "negative: " << negative << " "
              << "imaginary: " << imagEV << std::endl;
    if (!symetric || negative || imagEV)
      for (int i=0;i<N;++i)
        std::cout << std::scientific << lr( i ).toDouble() << " " 
                  << li(i).toDouble() << "i" << std::endl;
    std::cout << "Range: [" << minE << "," << maxE << "]" << std::endl;
  }
#elif HAVE_LAPACK 
  {
    const long int N = matrix.rows();
    assert( matrix.rows() == matrix.cols() );
    const char jobz = 'n'; // only calculate eigen values  
    const char uplo = 'u'; // use upper triangular matrix 
    const int dim = N;

    // length of matrix vector 
    const long int w = N * N ;

    // matrix to put into dsyev 
    double* matrixVector = new double [dim * dim]; 
    double* eigenvalues = new double [dim];

    // copy matrix  
    int row = 0;
    bool symetric = true ;
    for(int i=0; i<dim; ++i) 
    {
      for(int j=0; j<dim; ++j, ++row) 
      {
        // make sure matrix is symetric 
        if( ! ( std::abs( matrix( i , j ) -  matrix( j , i ) ) < 1e-12 ) )
        {
          symetric = false ; 
        }

        matrixVector[ row ] = matrix( i , j );
      }
    }

    symetric = false;

    // working memory 
    double* workSpace = new double [dim * dim]; 

    // return value information 
    long int info = 0;

    if( symetric ) 
    {
      // call LAPACK dsyev 
      dsyev_(&jobz, &uplo, &N, &matrixVector[0], &N, 
             &eigenvalues[0], &workSpace[0], &w, &info);
    }
    else 
    { 
      const char jobsvl = 'n';
      const char jobsvr = 'n';

      double *VL = new double [ N ];
      double *VR = new double [ N ];
      double *imEv = new double [ N ];

      dgeev_( &jobsvl, &jobsvr, &N, &matrixVector[0], &N, &eigenvalues[0], &imEv[0],
              &VL[0], &N, &VR[0], &N, &workSpace[0], &w, &info );

      delete [] VL;
      delete [] VR;
      delete [] imEv;
    }

    delete [] workSpace;
    delete [] matrixVector;

    std::cout << "Eigenvalues: " << std::endl;
    bool negative = false ; 
    double minE=1e10;
    double maxE=-1e10;
    for(int i=0; i<dim; ++i) 
    {
      if( eigenvalues[ i ] < 0 ) negative = true ;
      minE = std::min( minE, eigenvalues[ i ] );
      maxE = std::max( maxE, eigenvalues[ i ] );
      if (Parameter::verbose())
        std::cout << std::scientific << eigenvalues[ i ] << std::endl;
    }

    std::cout << "Range: [" << minE << "," << maxE << "]" << std::endl;

    // assert( ! negative );

    delete [] eigenvalues;

    if( info != 0 ) 
    {
      DUNE_THROW(InvalidStateException,"eigenValues: Eigenvalue calculation failed!");
    }
  }
#else 
#error
  DUNE_THROW(NotImplemented,"LAPACK is not availble, therefore no eigen value calculation");
#endif
}

} // end namespace FMatrixHelp 
} // end namespace Dune 
#endif
