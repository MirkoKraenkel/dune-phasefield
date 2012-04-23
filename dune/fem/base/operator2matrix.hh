#ifndef DUNE_FEM_OPERATOR2MATRIX_HH
#define DUNE_FEM_OPERATOR2MATRIX_HH

#include <dune/fem/operator/matrix/spmatrix.hh>

#if HAVE_ALGLIB
#include <alglib/evd.h>
#endif

namespace Dune {

template <class OperatorType>   
void analyseOperator(const OperatorType& op)
{
#if HAVE_ALGLIB
  typedef typename OperatorType :: DestinationType  DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType  SpaceType;

  const SpaceType& space = op.space();

  DestinationType vector("vec", space );
  DestinationType matrixCol("col", space );

  const int size = space.size();

  const unsigned int precision = 256;
  typedef amp::ampf<precision> Field;
  typedef ap::template_1d_array<Field> VectorType;
  typedef ap::template_2d_array<Field> DenseMatrixType;
  VectorType lr,li;
  lr.setbounds( 0, size-1 );
  li.setbounds( 0, size-1 );
  DenseMatrixType matrix,vl,vr;
  matrix.setbounds( 0, size-1, 0, size-1 );
  vl.setbounds( 0, size-1, 0, size-1 );
  vr.setbounds( 0, size-1, 0, size-1 );

  // test on linearity
  bool linearity = true;
  {
    DestinationType A("A", space );
    DestinationType B("B", space );
    DestinationType resA("A", space );
    DestinationType resB("B", space );
    DestinationType res3A("sqrt(3)A", space );
    DestinationType res3ApiB("sqrt(3)A+pi*B", space );
    for (int i=0;i<size;++i)
      A.leakPointer()[i] = (i%3)?sin(sqrt(2.)*double(i)):0;
    op(A,resA);
    for (int i=0;i<size;++i)
      B.leakPointer()[i] = cos(sqrt(M_PI)*sqrt(double(i)));
    op(B,resB);
    for (int i=0;i<size;++i)
      vector.leakPointer()[i] = sqrt(3.)*A.leakPointer()[i];
    op(vector,res3A);
    for (int i=0;i<size;++i)
      vector.leakPointer()[i] += M_PI*B.leakPointer()[i];
    op(vector,res3ApiB);
    for (int k=0;k<size;++k)
    {
      double res1 = res3A.leakPointer()[k] - sqrt(3.)*resA.leakPointer()[k];
      if (  std::abs(res1)>1e-10 &&
            std::abs(res1)>1e-10*std::abs(sqrt(3.)*resA.leakPointer()[k]) 
          )
        linearity = false;
      double res2 = res3ApiB.leakPointer()[k] - 
                    ( sqrt(3.)*resA.leakPointer()[k] + M_PI*resB.leakPointer()[k] );
      if ( std::abs(res2)>1e-10 &&
           std::abs(res2)>1e-10*std::abs(sqrt(3.)*resA.leakPointer()[k] + M_PI*resB.leakPointer()[k]) 
         )
        linearity = false;
    }
  }

  bool symetric = true ;

  vector.clear();  
  for(int i=0; i<size; ++i)
  {
    vector.leakPointer()[ i ] = 1;
    matrixCol.clear();
    op( vector, matrixCol );
    vector.leakPointer()[ i ] = 0;

    for( int j =0; j< size; ++j) 
    {
      const double value = matrixCol.leakPointer()[ j ];
      matrix(j,i) = -value;
      if (j<i) 
      {
        // test if matrix is symetric 
        double aij = matrix( i , j ).toDouble();
        double aji = matrix( j , i ).toDouble();
        if( std::abs(aij-aji) > 1e-10 * ( std::abs(aij)+std::abs(aji) )
            && std::abs(aij-aji) > 1e-10
            && std::abs(aij)+std::abs(aji) > 1e-10 )
        {
          std::cout << i << " " << j << " " << aij << " " << aji << std::endl;
          symetric = false ; 
        }
      }
    }
  }

  if (Parameter::verbose())
  {
    std::cout << "Matrix: " << std::endl;
    for(int i=0; i<size; ++i)
    {
      for( int j =0; j< size; ++j) 
        std::cout << -matrix(i,j).toDouble() << " ";
      std::cout << std::endl;
    }
  }

  bool converged = evd::rmatrixevd(matrix,size,1,lr,li,vl,vr);
  if (!converged)
  {
    std::cout << "Eigenvalue computation failed!" << std::endl;
    return;
  }

  std::cout << "Eigenvalues: " << std::endl;
  bool negative = false ; 
  bool imagEV = false;
  bool isEV[size];
  double minE=1e10;
  double maxE=-1e10;
  for(int i=0; i<size; ++i) 
  {
    if( lr( i ) < -1e-20 ) negative = true ;
    if ( std::abs( li(i).toDouble() ) > 1e-20 ) imagEV = true;
    minE = std::min( minE, lr( i ).toDouble() );
    maxE = std::max( maxE, lr( i ).toDouble() );
    if (Parameter::verbose())
      std::cout << std::scientific << lr( i ).toDouble() << " " 
                << li(i).toDouble() << "i" << std::endl;
    if ( std::abs(li(i).toDouble())<1e-12)
    {
      for(int j=0; j<size; ++j)
        vector.leakPointer()[ j ] = vr(j,i).toDouble();
      matrixCol.clear();
      op( vector, matrixCol );
      isEV[i] = true;
      for(int j=0; j<size; ++j)
      {
        if ( std::abs( vector.leakPointer()[j]*(-lr(i).toDouble())-
                       matrixCol.leakPointer()[j] ) > 1e-10 &&
             std::abs( vector.leakPointer()[j]*(-lr(i).toDouble())-
                       matrixCol.leakPointer()[j] ) > 1e-10*std::abs(vector.leakPointer()[j]*lr(i).toDouble())
            )
        {
          isEV[i] = false;
          break;
        }
      }
      if ( !isEV[i] )
      {
        std::cout << "error in EV " << i << std::endl;
        for(int j=0; j<size; ++j)
        {
          Field val = 0;
          for (int k=0;k<size; ++k) 
            val += (-matrix(j,k)) * vr(k,i);
          std::cout << vector.leakPointer()[j]*(-lr(i).toDouble()) << " "
                    << matrixCol.leakPointer()[j] << " " 
                    << val.toDouble() 
                    << std::endl;
        }
      }
    }
  }
  std::cout << "Symmetric: " << symetric << " "
            << "negative: " << negative << " "
            << "imaginary: " << imagEV << " "
            << "linearity: " << linearity 
            << std::endl;
  if ( (!symetric || negative || imagEV) && Parameter::verbose() )
    for (int i=0;i<size;++i) 
    {
      std::cout << std::scientific << lr( i ).toDouble() << " " 
                << li(i).toDouble() << "i" << "    ";
      for (int k=0;k<size;++k)
        std::cout << vr(k,i).toDouble() << " ";
      std::cout << std::endl;
    }
  std::cout << "Range: [" << minE << "," << maxE << "]" << std::endl;
#else
  typedef typename OperatorType :: DestinationType  DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType  SpaceType;

  const SpaceType& space = op.space();

  DestinationType vector("vec", space );
  DestinationType matrixRow("row", space );
  vector.clear();  

  const int size = space.size();

  Dune::SparseRowMatrix< double > matrix( size, size, 20 );
  matrix.clear();

  for(int i=0; i<size; ++i)
  {
    vector.leakPointer()[ i ] = 1;
    matrixRow.clear();
    op( vector, matrixRow );
    vector.leakPointer()[ i ] = 0;

    for( int j=0; j< size; ++j) 
    {
      const double value = matrixRow.leakPointer()[ j ];
      if( std::abs( value ) > 0 )
        matrix.add( i, j, -value );
      if( j<i ) 
      {
        if( std::abs( matrix( j, i ) ) - std::abs( value ) > 1e-14 )
        {
          std::cout << "Matrix non-symetric \n";
          abort();
        }
        else if( std::abs( matrix( j, i ) - value ) > 1e-14 )
        {
          std::cout << "Sign error \n";
          abort();
        }
      }
    }
    std::cout << std::endl;
  }

  if (Parameter::verbose())
    matrix.print( std::cout );

  //Fem :: eigenValues( matrix );
#endif
} // end of analyseOperator()

} // end of namespace Dune
#endif
