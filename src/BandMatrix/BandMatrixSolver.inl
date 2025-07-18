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

}




/// This function will call every methods from LAPACK to solve the system
template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::solve (Matrix& /*M*/, Vector& x, Vector& b)
{
    
    

}

template<class TMatrix, class TVector>
bool BandMatrixSolver<TMatrix,TVector>::addJMInvJt(linearalgebra::BaseMatrix* result, linearalgebra::BaseMatrix* J, SReal fact)
{
    
    return true;
}


/* 
    IDAMAX finds the index of the first element having maximum absolute value
    Parameters :
        N is INTEGER: number of elements in input vector(s)
        DX is DOUBLE PRECISION array, dimension (1 + (N - 1)* abs(INCX))
        INCX is INTEGER : storage spacing between elements of DX
*/
template<class TMatrix, class TVector>
typename BandMatrixSolver<TMatrix,TVector>::Index BandMatrixSolver<TMatrix,TVector>::idamax(Index n, Vector dx, Index incx)
{
    Index idmax = 0;
    if ((n < 1) || (incx <= 0))
    {
        idmax = 0;
        return idmax;
    }
    idmax = 1;
    if(n == 1) return idmax;
    if(incx == 1) /// code for increment equal to 1
    {
        double dmax = std::abs(dx[1]);
        for (Index i = 2; i < n; i++)
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
        double dmax = std::abs(dx[1]);
        ix += incx;
        for (Index i = 2; i < n; i++)
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

/*
    DGER performs the rank 1 operation

        A:= alpha * x * y**T + A

    where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix

    Parameters:
        M is INTEGER: On entry, M specifies the number of rows of the matrix A. M must be at least zero.
        N is INTEGER: On entry, N specifies the number of columns of the matrix A. N must be at least zero.
        ALPHA is DOUBLE PRECISION : On entry, ALPHA specifies the scalar alpha.
        X is DOUBLE PRECISION array, dimension at least
            ( 1 + ( m - 1 )*abs( INCX ) )
            Before entry, the incremented array X must contain the m element vector x
        INCX is INTEGER: On entry, INCX specifies the increment for the elements of X. INCX must not be zero.
        Y is DOUBLE PRECISION array, dimension at least
            ( 1 + ( n - 1 )*abs( INCY ) )
            Before entry, the incremented array Y must contain the n element vector y.
        INCY is INTEGER: On entry, INCY specifies the increment for the elements of Y. INCY must not be zero.
        A is DOUBLE PRECISION array, dimension ( LDA, N ): Before entry, the leading m by n part of the array A must
           contain the matrix of coefficients. On exit, A is overwritten by the updated matrix.
        LDA is INTEGER: On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
            LDA must be at least max( 1, m ).
 
 
*/
template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::dger(Index m, Index n, double alpha, Vector x, Index incx, Vector y, Index incy, Matrix a, Index lda)
{
    Index ix, jy, kx;
    double temp;
    
    if(incy > 0 ) jy=1;
    else jy = 1 - (n - 1)*incy;

    if(incx == 1)
    {
        for(Index j=1; j<n; j++)
        {
            if(y[jy] != 0 )
            {
                temp = alpha*y[jy];
                for(Index i = 1; i<m; i++)
                {
                    ///a[i][j] += x[i]*temp;
                    a.add(i,j,x[i]*temp);
                } 
            }
            jy += incy;
        }
    }
    else
    {
        if(incx > 0) kx = 1;
        else kx = 1 - (m - 1)*incx;

        for(Index j = 1; j<n; j++)
        {
            if(y[jy] != 0 )
            {
                temp = alpha*y[jy];
                ix = kx;
                for(Index i=1; i<m; i++)
                {
                    ///a[i][j] += x[ix]*temp;
                    a.add(i,j,x[ix]*temp);
                    ix += incx;
                }
            }
            jy += incy;
        }
    }


}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::dswap(Index n, Vector x, Index incx, Vector y, Index incy) 
{
    for (Index i = 0; i < n; ++i) {
        std::swap(x[i * incx], y[i * incy]);
    }
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::dscal(Index n, SReal a, Vector x, Index incx) {
    if (n <= 0 || incx <= 0) return;
    for (Index i = 0; i < n; ++i) {
        x[i * incx] *= a;
    }
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::dgbtf2(Index m, Index n, Index kl, Index ku, Matrix& ab, Index ldab, Vector& ipiv)
{
    Index jp, ju, km, kv;

    kv = ku + kl;

    for (Index j= ku + 2; j < std::min(kv,n); j++)
    {
        for (Index i= kv - j + 2; i<kl;i++)
        {
            ab.set(i,j,0);
        }
    }
    
    ju = 1; /// JU is the index of the last column affected by the current stage of the factorization

    for (Index j=1; j < std::min(m,n); j++)
    {
        /// Set fill-in elements in column J+KV to zero
        if((j+kv) <= n)
        {
            for (Index i = 1; i<kl; i++)
            {
                ///ab[i][j+kv] = 0;
                ab.set(i,j+kv,0);
            }
        }

        /// Find pivot and test for singularity. KM is the number of subdiagonal elements in the current column
        km = std::min(kl, m-j);
        ///jp = idamax(km+1, ab.element(kv+1,j), 1);
        jp = idamax(km+1, ab.getVector(kv+1,j), 1);
        ipiv[j] = jp + j - 1;
        if( ab.element(kv+jp,j) != 0)
        {
            ju = std::max(ju, std::min(j+ku+jp-1,n));
            /// Apply interchange to columns J to JU
            if(jp != 1)
            {
                ///dswap(ju-j+1, ab.element(kv+jp,j) , ldab-1, ab( kv+1, j ), ldab-1 );
                dswap(ju-j+1, ab.getVector(kv+jp,j) , ldab-1, ab.getVector( kv+1, j ), ldab-1 );

            }
            if(km > 0)
            {
                ///dscal( km, 1/ ab.element( kv+1, j ), ab.element( kv+2, j ), 1);
                dscal( km, 1/ ab.element( kv+1, j ), ab.getVector( kv+2, j ), 1);
            }
            if(ju>j)
            {
                Matrix AB = ab.getMatrix();
                ///dger(km, ju-j, -1, ab.element(kv+2,j), 1, ab.element(kv,j+1), ldab-1, ab.element(kv+1,j+1),ldab-1 );
                dger(km, ju-j, -1, ab.getVector(kv+2,j), 1, ab.getVector(kv,j+1), ldab-1, AB[kv+1][j+1],ldab-1 );
            }
        }
        
    }
} 



} //namespace sofa::component::linearsolver::direct
