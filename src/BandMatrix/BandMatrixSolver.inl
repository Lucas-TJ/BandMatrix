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
    LUFactorizedMatrix = M;
    std::cout << "Matrix LU Avant = " << LUFactorizedMatrix << std::endl;
    Index nbRows = LUFactorizedMatrix.rowSize();
    Index nbCols = LUFactorizedMatrix.colSize();
    Index nbSubDiag = nbCols / 2;
    Index nbSupDiag = nbCols / 2;
    Index dimLUMatrix = nbCols ;
    
    indexPivot.resize(std::min(nbRows,nbCols));
    this->computeLUBandMatrix(nbRows-1, nbCols-1, nbSubDiag, nbSupDiag, LUFactorizedMatrix, dimLUMatrix-1, indexPivot);
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "_____________________________________________________________________"<< std::endl;
    std::cout << "Matrix LU Après = " << LUFactorizedMatrix << std::endl;
}




/// This function will call every methods from LAPACK to solve the system
template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::solve (Matrix& /*M*/, Vector& x, Vector& b)
{


    Index nbSubDiag = LUFactorizedMatrix.colSize() / 2;
    Index nbSupDiag = LUFactorizedMatrix.colSize() / 2;
    
    Index matrixOrder = LUFactorizedMatrix.rowSize();
    Index nbColumnB = 1;
    Index dimB = b.size();

    B.resize(dimB,nbColumnB);
    for(Index i = 0; i<dimB;i++)
    {
        B.add(i,nbColumnB,b[i]); 
    }
    


    solveAxB(matrixOrder, nbSubDiag,nbSupDiag, nbColumnB, LUFactorizedMatrix, indexPivot, B, dimB);


    for(Index i = 0; i<dimB;i++)
    {
        x[i] = B.element(i,0); 
    }
    
    ///std::cout << "M = " << M.element(0,0) << std::endl;
    std::cout << "LUFactor = " << LUFactorizedMatrix.element(0,0) << std::endl;


}

template<class TMatrix, class TVector>
bool BandMatrixSolver<TMatrix,TVector>::addJMInvJt(linearalgebra::BaseMatrix* result, linearalgebra::BaseMatrix* J, SReal fact)
{
    return true;
}



template<class TMatrix, class TVector>
typename BandMatrixSolver<TMatrix,TVector>::Index BandMatrixSolver<TMatrix,TVector>::indexMax(Index nbElem, const Vector& dx, Index increment)
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
        Index ix = 1;
        double dmax = std::abs(dx[0]);
        ix += increment;
        for (Index i = 1; i < nbElem; i++)
        {
            if(std::abs(dx[ix]) > dmax)
            {
                idmax = i;
                dmax= std::abs(dx[ix]);
            }
            
        } 
    }

    return idmax;

}
template<class TMatrix, class TVector>
typename BandMatrixSolver<TMatrix,TVector>::Index BandMatrixSolver<TMatrix,TVector>::indexMax2(Index nbElem, const Matrix& M, Index increment, Index currentColumn, Index indexDiagP )
{
    Index idmax = 0;
    Index nbRows = M.rowSize();
    if ((nbElem < 1) || (increment <= 0))
    {
        return idmax;
    }
    idmax = 0;
    if (nbElem == 1) return idmax;
    if(increment == 1)
    {
        ///Real dmax = std::abs(M[indexDiagP][currentColumn]);
        Real dmax = std::abs(M.element(indexDiagP,currentColumn));
        for (Index i = indexDiagP; i < indexDiagP + nbElem; i++)
        {
            if(std::abs(M.element(i,currentColumn) > dmax))
            {
                idmax = i - indexDiagP;
                dmax = std::abs(M.element(i,currentColumn));
            }
        }
    }
    else
    {
        Index ix = indexDiagP;
        Real dmax = std::abs(M.element(indexDiagP,currentColumn));
        ix += increment;
        for (Index i = indexDiagP; i < indexDiagP + nbElem; i++)
        {
            if(std::abs(M.element(ix-1,currentColumn)) > dmax)
            {
                idmax = i;
                dmax = std::abs(M.element(ix-1,currentColumn));
            }
            ix += increment;
        }
    }
    return idmax;
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::rank1MatrixUpdate(Index nbRows, Index nbCols, Real alpha, Vector x, Index incrementX, Vector y, Index incrementY, Matrix& updatedMatrix)
{
    Index ix, jy, kx;
    Real temp;
    
    if(incrementY > 0 ) jy=0;
    else jy = 1 - (nbCols - 1)*incrementY;

    if(incrementX == 1)
    {
        for(Index j=0; j<nbCols; j++)
        {
            if(y[jy] != 0 )
            {
                temp = alpha*y[jy];
                for(Index i = 0; i<nbRows; i++)
                {
                    ///a[i][j] += x[i]*temp;
                    updatedMatrix.add(i,j,x[i]*temp);
                } 
            }
            jy += incrementY;
        }
    }
    else
    {
        if(incrementX > 0) kx = 1;
        else kx = 1 - (nbRows - 1)*incrementX;

        for(Index j = 0; j<nbCols; j++)
        {
            if(y[jy] != 0 )
            {
                temp = alpha*y[jy];
                ix = kx;
                for(Index i=0; i<nbRows; i++)
                {
                    ///a[i][j] += x[ix]*temp;
                    updatedMatrix.add(i,j,x[ix]*temp);
                    ix += incrementX;
                }
            }
            jy += incrementY;
        }
    }


}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::rank1MatrixUpdate2(Index nbRows, Index nbCols, Real alpha, Matrix& M, Index rowOfM1, Index colOfM1 ,Index incrementX, Index rowOfM2, Index colOfM2, Index incrementY, Matrix& updatedMatrix)
{   // rank1MatrixUpdate(nbSubDiagElements, ju-j, -1, LUFactorizedMatrix.getVector(nbSupDiagU+1,j), 1, LUFactorizedMatrix.getVector(nbSupDiagU-1,j+1), dimLUMatrix-1, LUFactorizedMatrix )
    // rank1MatrixUpdate(Index nbRows, Index nbCols, Real alpha, Vector x, Index incrementX, Vector y, Index incrementY, Matrix& updatedMatrix)

    Index ix, jy, kx;
    Real temp;
    
    Vector x ;
    Vector y ;

    x.resize(std::abs(nbRows-rowOfM1));
    y.resize(std::abs(nbRows-rowOfM2));
    for(Index i = 0 ; i < nbRows-rowOfM1 ; i++)
    {
        x.set(i,M.element(i+rowOfM1,colOfM1));
        y.set(i,M.element(i+rowOfM2,colOfM2));
    }
    

    if(incrementY > 0 ) jy=0;
    else jy = 1 - (nbCols - 1)*incrementY ;

    if(incrementX == 1)
    {
        
        
        for(Index j=0; j<nbCols; j++)
        {
            
            //std::cout <<"l 268" <<std::endl;
            //if(y[jy] != 0 )
            if( M.element(jy+rowOfM2,colOfM2) != 0)
            {
                //std::cout <<"l 275" <<std::endl;
                
                //temp = alpha*y[jy];
                // std::cout <<"jy = " << jy << std::endl;
                // std::cout <<"y size" << nbRows-rowOfM2<<std::endl;
                temp = alpha * y.element(jy);
                // std::cout << "temp = " << temp << std::endl;
                for(Index i = 0; i<nbRows; i++)
                {
                    
                    // std::cout <<"l 282" <<std::endl;
                    // ///a[i][j] += x[i]*temp;
                    // //updatedMatrix.add(i,j,x[i]*temp);
                    // std::cout << "i =  " <<i<< std::endl;
                    // std::cout << "j =  " <<j<< std::endl;
                    // std::cout << "a(i,j) = " << M.element(i,j) << std::endl;
                    // std::cout << "x[i] = " << M.element(i,colOfM1)<<std::endl;
                    
                    //M.add(i,j, M.element(i+rowOfM1,colOfM1) * temp );
                    M.add(i,j, x.element(i) * temp );
                    // std::cout << "matrix is " << M << std::endl;
                    // std::cout <<"l 290" <<std::endl;
                } 
            }
            jy += incrementY;
        }
    }
    else
    {
        
        
        if(incrementX > 0) kx = 1;
        else kx = 1 - (nbRows - 1)*incrementX;

        for(Index j = 0; j<nbCols; j++)
        {   
            
            
            //if(y[jy] != 0 )
            if(M.element(jy+rowOfM2,colOfM2) != 0)
            {
                
                
                //temp = alpha*y[jy];
                temp = alpha * M.element(jy+rowOfM2,colOfM2);
                ix = kx;
                for(Index i=0; i<nbRows; i++)
                {
                    
                    
                    ///a[i][j] += x[ix]*temp;
                    //updatedMatrix.add(i,j,x[ix]*temp);
                    updatedMatrix.add(i,j, M.element(ix+rowOfM1,colOfM1) *temp );
                    ix += incrementX;
                }
            }
            jy += incrementY;
        }
    }

}


template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::swapVector(Index n, Vector x, Index incx, Vector y, Index incy) 
{
    for (Index i = 0; i < n; ++i) {
        std::swap(x[i * incx], y[i * incy]);
    }
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::swapVector2(Index nbElemToSwap, Matrix& M, Index currentColumn, Index indexOfM1, Index incrementM1, Index indexOfM2, Index incrementM2 ) 
{
    Index ix = 0;
    Index iy = 0; 
    if(incrementM1 < 0)
    {
        ix = (-nbElemToSwap+1)*incrementM1+1;
        iy = (-nbElemToSwap+1)*incrementM2+1;
    }
    Index ns = nbElemToSwap/2;
    for(Index i=0; i< ns;i++)
    {
        Index dtemp = M.element(indexOfM1+ix,currentColumn);
        M.set(indexOfM1+ix, currentColumn, M.element(indexOfM2+iy,currentColumn));
        M.set(indexOfM2+iy, currentColumn, dtemp);
        
        ix += incrementM1;
        iy += incrementM2;
    }

}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::scaleVector(Index n, SReal a, Vector x, Index incx) 
{
    if (n <= 0 || incx <= 0) return;
    for (Index i = 0; i < n; ++i) {
        x[i * incx] *= a;
    }
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::scaleVector2(Index nbElem, SReal a, Matrix& M, Index incx, Index currentColumn, Index indexDiagP) 
{

    if (nbElem <= 0 || incx <= 0) return;
    for (Index i = indexDiagP; i < nbElem+indexDiagP; ++i) {
        ///x[i * incx] *= a;
        M.set(i*incx, currentColumn, M.element(i*incx, currentColumn)*a);
    }
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::scaleVector3(Index nbElem, SReal a, Matrix& M, Index incx, Index currentColumn, Index indexDiagP) 
{


    // Index m;
    // if(incx == 1)
    // {
    //     std::cout<<"IF"<<std::endl;
    //     m = nbElem%5;
    //     if(m =!0)
    //     {
    //         std::cout <<"IF - IF" << std::endl;
    //         for(Index i = indexDiagP; i - indexDiagP < m; i++)
    //         {
    //             std::cout <<"IF - IF - For" << std::endl;
    //             M.set(i,currentColumn, a*M.element(i,currentColumn));
    //         }
    //         //if(nbElem < 5) return;
    //     }
    //     Index mp1 = m+5;
    //     Index nincx = nbElem*incx;
    //     for(Index i= indexDiagP; i-indexDiagP< nincx ; i += incx)
    //     {
    //         std::cout <<"IF - FOR" << std::endl;
    //         std::cout << "i = " << i<< std::endl;
    //         std::cout << "M[] = " << M.element(i,currentColumn)<<std::endl;
    //         M.set(i,currentColumn, a*M.element(i,currentColumn));
            
    //     }
    // }
    // else
    
        // std::cout<<"ELSE"<<std::endl;
        Index nincx = nbElem*incx;
        for(Index i=indexDiagP ; i - indexDiagP < nincx; i += incx)
        {
            // std::cout << "i = " << i<< std::endl;
            // std::cout << "M[] = " << M.element(i,currentColumn)<<std::endl;
            M.set(i,currentColumn, a*M.element(i,currentColumn));
        }
    
    
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::computeLUBandMatrix(Index nbRows, Index nbCols, Index nbSubDiag, Index nbSupDiag, Matrix& LUFactorizedMatrix, Index dimLUMatrix, Vector& indexPivot)
{
    Index ju;

    Index nbSupDiagU = nbSupDiag + nbSubDiag ;

    std::cout << LUFactorizedMatrix << std::endl;

    // obliger de resize
    Index oldbandwidth = nbSupDiag + nbSubDiag + 1;
    Index newbandwidth = nbSupDiagU + 1 + nbSubDiag;
    Index nbAdd = (newbandwidth - oldbandwidth) /2 ;
    LUFactorizedMatrix.resize(newbandwidth,nbCols);
    for(Index i = oldbandwidth-1; i>-1;i--)
    {
        for(Index j = 0; j<nbCols;j++)
        {
            LUFactorizedMatrix.set(i+nbAdd*2,j,LUFactorizedMatrix.element(i,j));
            LUFactorizedMatrix.set(i,j,0);
        }
    }
    //fin du resize

    std::cout << LUFactorizedMatrix << std::endl;
    // for (Index j= nbSupDiag + 1 ; j < std::min(nbSupDiagU,nbCols); j++)    /// Set fill-in elements in columns KU+2 to KV to zero
    // {
    //     for (Index i= nbSupDiagU - j + 1; i< nbSubDiag;i++)
    //     {
    //         std::cout << "i = " << i<<std::endl;
    //         std::cout<<"j = "<<j<<std::endl;
    //         std::cout<<"TEST "<<std::endl;
    //         LUFactorizedMatrix.set(i,j,55);
    //     }
    // }
    // std::cout << LUFactorizedMatrix << std::endl;
    
    ju = 1; /// JU is the index of the last column affected by the current stage of the factorization

    for (Index j=1; j < std::min(nbRows+1,nbCols+1)-1; j++)
    {
        std::cout << "PASSAGE DE J : " << j << std::endl;
        std::cout << "std::min(nbRow,nbCol) : " << std::min(nbRows,nbCols) << std::endl;
        
        /// Set fill-in elements in column J+KV to zero
        // if((j+nbSupDiagU) < nbCols)
        // {
        //     for (Index i = 0; i< nbSubDiag; i++)
        //     {
        //         std::cout << "i = " << i<<std::endl;
        //         std::cout<<"j = "<<j<<std::endl;
        //         std::cout<<" j+nbSupDiagU = "<<j+nbSupDiagU<<std::endl;
        //         ///ab[i][j+kv] = 0;
        //         LUFactorizedMatrix.set(i,j+nbSupDiagU,0);
        //     }
        // }
        
        /// Find pivot and test for singularity. KM is the number of subdiagonal elements in the current column km = std::min(kl, m-j)
        Index nbSubDiagElements = std::min(nbSubDiag, nbRows-j);
        std::cout << "TEST 1 " << std::endl;
        ///jp = idamax(km+1, ab.element(kv+1,j), 1);
        ///jp = idamax(km+1, ab.getVector(kv+1,j), 1); /// jp ? Maybe index max for the pivot

        /// Trouver un moyen de naviguer entre la diagonal principal et le bas de la matrice dans LUFactorizedMatrix 
        ///Index j_pivot = indexMax(nbSubDiagElements + 1, *(ptr + nbSupDiagU + (j-1) * dimLUMatrix) ,1); 
        Index j_pivot = indexMax2(nbSubDiagElements+1, LUFactorizedMatrix, 1, j, nbSupDiagU);
        std::cout << "TEST 2 " << std::endl;
        indexPivot[j-1] = j_pivot + j - 1 ;
        if( LUFactorizedMatrix.element(nbSupDiagU+j_pivot,j) != 0)
        {
            std::cout << "TEST 3 " << std::endl;
            ju = std::max(ju, std::min(j+nbSupDiag+j_pivot-1,nbCols));
            /// Apply interchange to columns J to JU
            if(j_pivot != 1)
            {
                std::cout << "TEST 4 " << std::endl;
                ///Uniquement dans la partie L, on échange partiellement les vecteurs :
                /// Idée : accéder à la matrice, se placer dans les deux colonnes et intervertir ?
                ///dswap(ju-j+1, ab.element(kv+jp,j) , ldab-1, ab( kv+1, j ), ldab-1 );
                ///swapVector(ju-j+1, LUFactorizedMatrix.getVector(nbSupDiagU+jp,j) , dimLUMatrix-1, LUFactorizedMatrix.getVector( kv+1, j ), dimLUMatrix-1 ); /// Forcément utiliser pointeur

                swapVector2(ju-j+1,LUFactorizedMatrix, j,nbSupDiagU+j_pivot,dimLUMatrix-1,nbSupDiagU,dimLUMatrix-1);
            }
            if(nbSubDiagElements > 0)
            {
                std::cout << "TEST 5 " << std::endl;
                ///dscal( km, 1/ ab.element( kv+1, j ), ab.element( kv+2, j ), 1);

                ///scaleVector( nbSubDiagElements, 1/ LUFactorizedMatrix.element( nbSupDiagU+1, j ), LUFactorizedMatrix.getVector( nbSupDiagU+2, j ), 1);
                scaleVector3( nbSubDiagElements, 1/ LUFactorizedMatrix.element( nbSupDiagU, j-1 ), LUFactorizedMatrix, 1, j, nbSupDiagU+1);
                if(ju>j)
                {
                    ///dger(km, ju-j, -1, ab.element(kv+2,j), 1, ab.element(kv,j+1), ldab-1, ab.element(kv+1,j+1),ldab-1 );
                    //std::cout << LUFactorizedMatrix << std::endl;
                    //std::cout << "nbSupDiagU = " << nbSupDiagU << std::endl;
                    //std::cout << "j = " << j << std::endl;
                    //std::cout << "vector x = " << LUFactorizedMatrix.getVector(nbSupDiagU+1,j) << std::endl;
                    //std::cout << "vector y = " << LUFactorizedMatrix.getVector(nbSupDiagU-1,j+1) << std::endl;
                    std::cout << "TEST 6 " << std::endl;
                    std::cout <<"MATRIX JUSTE AVANT RANK1 = "<< LUFactorizedMatrix << std::endl;
                    ///rank1MatrixUpdate(nbSubDiagElements, ju-j, -1, LUFactorizedMatrix.getVector(nbSupDiagU+1,j), 1, LUFactorizedMatrix.getVector(nbSupDiagU-1,j+1), dimLUMatrix-1, LUFactorizedMatrix );                    
                    rank1MatrixUpdate2(nbSubDiagElements, ju-j, -1, LUFactorizedMatrix, nbSupDiagU+2,j , 1, nbSupDiagU, j+1, dimLUMatrix-1, LUFactorizedMatrix );
                }
            }
        }
        std::cout << "TEST 7 je crois que c'est ici" << std::endl;
    }
} 

template<class TMatrix, class TVector>
///solveUxB(n,k,a,lda,x,incx)
void BandMatrixSolver<TMatrix,TVector>::solveUxB(Index matrixOrder, Index nbSupDiag, const Matrix& Matrix, Vector& x, Index incrementX ) 
{
    /*
    Set up the start point in X if the increment is not unity. This will be (N - 1)*incx too small for descending loops
    */
    Index kx;
    if(incrementX <= 0)
    {
        kx = 1 - (matrixOrder-1)*incrementX;
    }
    else
    {
        kx = 0;
    }

    /*
    Start the operations. In this version the elements of A are accessed by sequentially with one pass through A
    */

    Index kplus1 = nbSupDiag + 1;
    if (incrementX == 1)
    {
        for (Index j=matrixOrder; j>1 ; j--)
        {
            if( x[j] != 0 )
            {
                Index l = kplus1 - j;
                x[j] = x[j]/Matrix.element(kplus1,j);
                Index temp = x[j];
                for (Index i =j-1; i> std::max(1,j-nbSupDiag); i--)
                {
                    x[i]= x[i]-temp*Matrix.element(l+i,j);
                }
            }
        } 
    }
    else
    {
        kx = kx + (matrixOrder-1)*incrementX;
        Index jx = kx;
        for(Index j = matrixOrder; j>1; j--)
        {
            kx = kx - incrementX;
            if(x[jx] != 0)
            {
                Index ix = kx;
                Index l = kplus1 - j;
                x[jx] = x[jx]/ Matrix(kplus1,j);
                Index temp = x[j];
                for(Index i=j-1; std::max(matrixOrder,j+nbSupDiag); i--)
                {
                    x[ix]= x[ix] - temp*Matrix.element(l+i,j);
                    ix = ix - incrementX; 
                }
            }
            jx = jx - incrementX;
        }
    }
}

template <class TMatrix, class TVector>
void BandMatrixSolver<TMatrix, TVector>::solveAxB(Index matrixOrder, Index nbSubDiag, Index nbSupDiag, Index nbColumnB, Matrix& LUFactorizedMatrix, Vector indexPivot, Matrix& B, Index dimB)
{
    Index nbRow = LUFactorizedMatrix.rowSize();
    
    /*
    Solve L*x=B, overwriting B with X

    L is represented as a product of permutations and unit lower triangular matrices L = P(1) * L(1)* ... *P(n-1) * L(n-1)
    where each transformation L(i) is a rank-one modification of the identity matrix
    */
    if(nbSubDiag > 0)
    {
        for(Index j=0; j <= matrixOrder-1; j++)
        {
            Index lm = std::min(nbSubDiag, matrixOrder-j);
            Index l = indexPivot[j];
            if(l!=j)
            {
                swapVector2(nbColumnB, B, 0, l, dimB, j, dimB);
            }
            rank1MatrixUpdate(lm, nbColumnB, -1, LUFactorizedMatrix.getVector(nbRow+1,j), 1, B.getVector(j,0), dimB, B);
        }
    }
    for(Index i=1; i < nbColumnB; i++)
    {
        Vector BVec = B.getVector(0,i);
        solveUxB(matrixOrder, nbSubDiag + nbSupDiag, LUFactorizedMatrix, BVec,1);
    }
}



} //namespace sofa::component::linearsolver::direct
