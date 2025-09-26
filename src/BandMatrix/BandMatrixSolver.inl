/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once
#include <BandMatrix/BandMatrixSolver.h>
#include <sofa/linearalgebra/FullMatrix.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/linearalgebra/SparseMatrix.h>
#include <sofa/core/ObjectFactory.h>
#include <iostream>
#include <sofa/helper/system/thread/CTime.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <cmath>

namespace sofa::component::linearsolver::direct
{

template<class TMatrix, class TVector>
BandMatrixSolver<TMatrix,TVector>::BandMatrixSolver()
{
}

/// This function will call every methods from LAPACK (dbgtf2, idamax, dger etc...) to compute LU factorization
template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::invert(Matrix& M)
{   
    sofa::helper::AdvancedTimer::stepBegin("factorization");
    LUFactorizedMatrix = M;
    
    
    Index nbRows = LUFactorizedMatrix.rowSize();
    Index nbCols = LUFactorizedMatrix.colSize();
    Index nbSubDiag = nbCols / 2;
    Index nbSupDiag = nbCols / 2;
    Index dimLUMatrix = nbCols ;
    indexPivot.resize(std::min(nbRows,nbCols));
    computeLUBandMatrix(nbRows, nbCols, nbSubDiag, nbSupDiag, LUFactorizedMatrix, dimLUMatrix-1, indexPivot);
    sofa::helper::AdvancedTimer::stepEnd("factorization");
    
}

/// This function will call every methods from LAPACK to solve the system
template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::solve (Matrix& /*M*/, Vector& x, Vector& b)
{
    
    Index nbSubDiag = LUFactorizedMatrix.colSize() / 2;
    Index nbSupDiag = LUFactorizedMatrix.colSize() / 2;
    Index matrixOrder = LUFactorizedMatrix.colSize();

    Index nbColumnB = 1;
    Index dimB = b.size();

    Matrix B;
    B.resize(dimB,1);
    for(Index i = 0; i<dimB; i++)
    {
        B.set(i,0,b[i]);
        
    }
    solveAxB(matrixOrder, nbSubDiag,nbSupDiag, nbColumnB, LUFactorizedMatrix, indexPivot, B, dimB);
    for(Index i = 0; i<dimB;i++)
    {
        x[i] = B.element(i,0); 
        
    }
}

template<class TMatrix, class TVector>
bool BandMatrixSolver<TMatrix,TVector>::addJMInvJt(linearalgebra::BaseMatrix* result, linearalgebra::BaseMatrix* J, SReal fact)
{
    return false;
}



template<class TMatrix, class TVector>
typename BandMatrixSolver<TMatrix,TVector>::Index BandMatrixSolver<TMatrix,TVector>::indexMax(Index nbElem, Real * dx, Index increment)
{
    Index idmax = 0;
    if ((nbElem < 1) || (increment <= 0))
    {
        return idmax;
    }
    idmax = 0;
    if(nbElem == 1) return idmax;
    if(increment == 1) /// code for increment equal to 1
    {
        Real dmax = std::abs(dx[0]);
        for (Index i = 1; i < nbElem; i++)
        {
            if(std::abs(dx[i]) > dmax)
            {
                idmax = i;
                dmax = std::abs(dx[i]);
            }
        }
    }
    else /// code for increment not equal to 1
    {
        Index ix = 0;
        Real dmax = std::abs(dx[0]);
        ix += increment;
        for (Index i = 1; i < nbElem; i++)
        {
            if(std::abs(dx[ix]) > dmax)
            {
                idmax = i;
                dmax= std::abs(dx[ix]);
            }
            ix += increment;
        } 
    }

    return idmax;

}


template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::rank1MatrixUpdate(
    Index nbRows, Index nbCols, Real alpha, Real * x, Index incrementX, Real * y, Index incrementY, Real * updatedMatrix)
{
    Index ix, jy, kx;
    Real temp;
    
    if(incrementY > 0 ) jy=0;
    else jy = 1 - (nbCols - 1)*incrementY - 1;

    if(incrementX == 1)
    {
        for(Index j=0; j<nbCols; j++)
        {
            if(y[jy] != 0 )
            {
                temp = alpha*y[jy];
                for(Index i = 0; i<nbRows; i++)
                {
                    updatedMatrix[i*nbCols + j] += x[i]*temp;
                } 
            }
            jy += incrementY;
        }
    }
    else
    {
        if(incrementX > 0) kx = 0;
        else kx = 1 - (nbRows - 1)*incrementX - 1;

        for(Index j = 0; j<nbCols; j++)
        {
            if(y[jy] != 0 )
            {
                temp = alpha*y[jy];
                ix = kx;
                for(Index i=0; i<nbRows; i++)
                {
                    updatedMatrix[i*nbCols + j] += x[ix]*temp;
                    std::cout << "x(ix) = " << x[ix] <<" ,and temp = "<<temp<<std::endl;
                    ix += incrementX;
                }
            }
            jy += incrementY;
        }
    }


}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::swapVector(Index n, Real * x, Index incx, Real * y, Index incy) 
{
    
    if(n <= 0 ) return;
    // code for both increments equal to 1
    if(incx == 1 && incy == 1)
    {
        Index m = n%3;
        if(m != 0)
        {
            for(Index i = 0 ; i<m ; i++)
            {
                Real dtemp = x[i];
                x[i] = y[i];
                y[i] = dtemp;
            }
            if(n<3) return;
        }
        Index mp1 = m + 1;
        for(Index i = mp1 - 1; i<n; i+=3 )
        {
            Real dtemp = x[i];
            x[i]=y[i];
            y[i]=dtemp;
            dtemp = x[i+1];
            x[i+1]=y[i+1];
            y[i+1]=dtemp;
            dtemp = x[i+2];
            x[i+2]=y[i+2];
            y[i+2]=dtemp;
            
        }
    }
    // code for unequal increments or equal increments not equal to 1
    else
    {
        Index ix = 0;
        Index iy = 0;
        if(incx < 0) ix = (-n+1)*incx;
        if(incy < 0) iy = (-n+1)*incy;
        for(Index i = 0; i<n; i++)
        {
            Real dtemp = x[ix];
            x[ix] = y[iy];
            y[iy] = dtemp;
            ix += incx;
            iy += incy;
        }
    }
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::scaleVector(Index n, Real a, Real * x, Index incx) 
{
    if(n <= 0 || incx <= 0 || a == 1) return;
    
    // code for increment equal to 1
    if(incx == 1)
    {
        Index m = n % 5;
        if(m != 0)
        {
            for(Index i = 0; i<m ; i++)
            {
                x[i] = a * x[i];
            }
            if(n < 5) return;
        }
        Index mp1 = m + 1;
        for(Index i = mp1 - 1 ; i<n ; i += 5 )
        {
            x[i] *= a;
            x[i+1] *= a;
            x[i+2] *= a;
            x[i+3] *= a;
            x[i+4] *= a;
        }
    }
    // code for increment not equal to 1
    else
    {
        Index nincx = n*incx;
        for(Index i = 0; i<nincx; i += incx)
        {
            x[i] *= a; 
        }
    }
    return ;
}



template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::computeLUBandMatrix(
    Index nbRows, Index nbCols, Index nbSubDiag, Index nbSupDiag, Matrix& LUFactorizedMatrix, Index dimLUMatrix, Vector& indexPivot)
{

    Index ju;

    Index nbSupDiagU = nbSupDiag + nbSubDiag ;
    
    sofa::type::vector<Real> TestMatrix;
    TestMatrix.resize(nbCols * (nbSupDiagU+1+nbSubDiag));
    
    Index colTestMatrix = nbSupDiagU + 1 + nbSubDiag;
    
    for(Index i=0 ; i< nbRows; i++)
    {

        for(Index j = nbSupDiagU+1; j>-1; j--)
        {
            TestMatrix[i*colTestMatrix+j+nbSubDiag-1] = LUFactorizedMatrix.element(j,i);
        }
    }

    ju = 0 ;
    
    
    Real * ptr = TestMatrix.data() ;
    for(Index j = 0; j < std::min(nbRows,nbCols); j++ )
    {
        Index nbSubDiagElements = std::min(nbSubDiag , (nbRows-1) - j);
        Index j_pivot = indexMax(nbSubDiagElements, ptr + colTestMatrix*j+nbSupDiagU ,1);
        indexPivot[j] = j_pivot + j;

        if ( TestMatrix[j*colTestMatrix + nbSupDiagU+j_pivot] != 0)
        {
            ju = std::max(ju, std::min(j+nbSupDiag+j_pivot,nbRows));
            if(j_pivot != 0)
            {
                swapVector(ju-j+1, ptr + colTestMatrix*j + (nbSupDiagU+j_pivot-1), dimLUMatrix-1, ptr + colTestMatrix*j +nbSupDiagU,dimLUMatrix-1);
            }
            if(nbSubDiagElements > 0)
            {
                Real valdiv = *(ptr + colTestMatrix*j + nbSupDiagU);
                scaleVector(nbSubDiagElements, 1./valdiv , ptr + colTestMatrix*j+(nbSupDiagU+1) ,1 );

                if(ju>j)
                {
                    
                    rank1MatrixUpdate(nbSubDiagElements, ju-j, -1, ptr + colTestMatrix*(j) + (nbSupDiagU+1), 1, ptr + colTestMatrix*(j+1) + (nbSupDiagU-1), dimLUMatrix-1, ptr + colTestMatrix*(j+1) + (nbSupDiagU));
                }
            }
        }
        else
        {
            for(Index i = 0; i < nbCols; i++)
            {
                for(Index j = 0; j < colTestMatrix; j++)
                {
                    LUFactorizedMatrix.set(j,i,TestMatrix[i*colTestMatrix + j]);
                }
            }
            return;
        }

        
    }
    
    for(Index i = 0; i < nbCols; i++)
    {
        for(Index j = 0; j < colTestMatrix; j++)
        {
            LUFactorizedMatrix.set(j,i,TestMatrix[i*colTestMatrix + j]);
        }
    }

    
    

}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::solveUxB(Index matrixOrder, Index nbSupDiag, Real * Matrix, Real * x, Index incrementX, Index nbCol ) 
{
    /*
    Set up the start point in X if the increment is not unity. This will be (N - 1)*incx too small for descending loops
    */
    Index kx;
    if(incrementX <= 0)
    {
        kx = 1 - (matrixOrder-1)*incrementX - 1;
    }
    else if(incrementX != 1)
    {
        kx = 0;
    }


    /*
    Start the operations. In this version the elements of A are accessed by sequentially with one pass through A
    */

    Index kplus1 = nbSupDiag + 1;

    if (incrementX == 1)
    {
        for (Index j=matrixOrder-1; j>=0 ; j--)
        {
            if( x[j] != 0 )
            {
                Index l = kplus1 - j -1;
                if(Matrix[nbCol*(kplus1-1)+j] != 0) x[j] = x[j] /Matrix[nbCol*(kplus1-1)+j];
                Real temp = x[j];
                for (Index i =j-1; i>= std::max(1,j-nbSupDiag)-1; i--)
                {
                    x[i] = x[i] - temp*Matrix[nbCol*(l+i) + j];
                }
            }
        } 
    }
    else
    {
        kx = kx + (matrixOrder-1)*incrementX - 1;
        Index jx = kx;
        for(Index j = matrixOrder-1; j>=0; j--)
        {
            kx = kx - incrementX;
            if(x[jx] != 0)
            {
                Index ix = kx;
                Index l = kplus1 - j;
                if (Matrix[nbCol*(kplus1-1)+j] != 0) x[jx] = x[jx] /Matrix[nbCol*(kplus1-1)+j];
                Real temp = x[j];
                for(Index i=j-1; i >= std::max(matrixOrder,j+nbSupDiag)-1; i--)
                {
                    x[ix]= x[ix] - temp*Matrix[(l+i)*nbCol+j];
                    ix = ix - incrementX; 
                }
            }
            jx = jx - incrementX;
        }
    }
}

template <class TMatrix, class TVector>
void BandMatrixSolver<TMatrix, TVector>::solveAxB(
    Index matrixOrder, Index nbSubDiag, Index nbSupDiag, Index nbColumnB, Matrix& LUFactorizedMatrix, Vector indexPivot, Matrix& B, Index dimB)
{
    Index nbRow = LUFactorizedMatrix.rowSize();
    Index nbCols = LUFactorizedMatrix.colSize();
    Index nbSupDiagU = nbSubDiag + nbSupDiag;
    
    /*
    Solve L*x=B, overwriting B with X

    L is represented as a product of permutations and unit lower triangular matrices L = P(1) * L(1)* ... *P(n-1) * L(n-1)
    where each transformation L(i) is a rank-one modification of the identity matrix
    */
    

    sofa::type::vector<Real> TestMatrix;
    TestMatrix.resize(nbCols * (nbSupDiagU+1+nbSubDiag));
    Index colTestMatrix = nbSupDiagU + 1 + nbSubDiag;    

    for(Index i=0 ; i< matrixOrder; i++)
    {
        for(Index j = nbSupDiagU+1; j>-1; j--)
        {
            TestMatrix[i*colTestMatrix+j+nbSubDiag-1] = LUFactorizedMatrix.element(j,i);
        }
    }

    Real * ptr = TestMatrix.data();

    sofa::type::vector<Real> B_v;
    B_v.resize(B.rowSize());
    
    Index SizeB_v = B_v.size();
    for(Index i=0; i< SizeB_v; i++)
    {
        B_v[i] = B.element(i,0);
    }
    Real * ptr2 = B_v.data();

    if(nbSubDiag > 0)
    {
        
        for(Index j=0; j < matrixOrder-1; j++)
        {
            
            Index lm = std::min(nbSubDiag, matrixOrder-j);
            
            Index l = indexPivot[j];
            if(l!=j)
            {
                swapVector(1, ptr2+l, dimB, ptr2 + j, dimB);

            }
            rank1MatrixUpdate(lm, nbColumnB, -1., ptr + j*colTestMatrix + nbSupDiagU+1, 1, ptr2 + j, dimB, ptr2 + j + 1);
        }
    }
    

    sofa::type::vector<Real> LUFact;
    LUFact.resize(nbCols * nbRow);
    

    for(Index i=0 ; i< nbRow; i++)
    {
        for(Index j = 0; j<nbCols; j++)
        {
            LUFact[i*nbRow+j] = LUFactorizedMatrix.element(i,j);
        }
    }
    Real * ptr3 = LUFact.data() ;

    solveUxB(matrixOrder, nbSubDiag + nbSupDiag, ptr3, ptr2 ,1,nbRow);

    for(Index i=0; i< SizeB_v; i++)
    {
        B.set(i,0,B_v[i]);
    }
}


} //namespace sofa::component::linearsolver::direct
