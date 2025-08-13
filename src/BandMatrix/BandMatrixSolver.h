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
#include <BandMatrix/config.h>
///#include <sofa/component/linearsolver/direct/config.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/component/linearsolver/iterative/MatrixLinearSolver.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/linearalgebra/SparseMatrix.h>
#include <BandMatrix/BandMatrix.h>

#include <sofa/linearalgebra/FullMatrix.h>
#include <cmath>
#include <sofa/type/Mat.h>
#include <sofa/type/vector.h>

namespace sofa::component::linearsolver::direct
{
/// Linear system solver using Gaussian elimination for Banded matrices

template<class TMatrix, class TVector>
class BandMatrixSolver : public sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector>
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE2(BandMatrixSolver,TMatrix,TVector), 
        SOFA_TEMPLATE2(sofa::component::linearsolver::MatrixLinearSolver, TMatrix, TVector));

    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef typename Vector::Real Real;
    typedef sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector> Inherit;

    typedef typename linearalgebra::BaseMatrix::Index Index;

    BandMatrixSolver();

    void reset() override {}

    ///typename Matrix::InvMatrixType Minv;  
    ////////////////////////////

protected:
    
    Matrix LUFactorizedMatrix;
    Matrix B;
    Vector indexPivot;
    
public:

    void invert(Matrix& M) override;

    /// Solve Mx=b
    void solve (Matrix& /*M*/, Vector& x, Vector& b) override;

    /// Multiply the inverse of the system matrix by the transpose of the given matrix, and multiply the result with the given matrix J
    ///
    /// @param result the variable where the result will be added
    /// @param J the matrix J to use
    /// @return false if the solver does not support this operation, of it the system matrix is not invertible
    bool addJMInvJt(linearalgebra::BaseMatrix* result, linearalgebra::BaseMatrix* J, SReal fact) override;

    /// @brief Finds the index of the first element having maximum absolute value
    /// @param nbElem the number of elements in vector
    /// @param dx the vector
    /// @param increment storage spacing between elements of dx
    /// @return The index max of the vector
    static Index indexMax(Index nbElem,  Real * dx, Index increment);

    /// @brief This function performs the rank 1 operation.
    ///        A := alpha * x * y**T + A
    ///        where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix
    /// @param nbRows the number of Rows
    /// @param nbCols the number of Columns
    /// @param alpha a scalar
    /// @param x vector x (as pointer)
    /// @param incrementX increment of x 
    /// @param y vector y (as pointer)
    /// @param incrementY increment of y 
    /// @param updatedMatrix Matrix which we want to update (as pointer)
    static void rank1MatrixUpdate(Index nbRows, Index nbCols, Real alpha, Real * x, Index incrementX, Real * y, Index incrementY, Real * updatedMatrix);

    /// @brief This function performs the rank 1 operation.
    ///        A := alpha * x * y**T + A
    ///        where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix
    /// @param n number of elements
    /// @param x vector x (as a pointer)
    /// @param incx increment of x
    /// @param y vector x (as a pointer)
    /// @param incy increment of y 
    static void swapVector(Index n, Real * x, Index incx, Real * y, Index incy); 

    /// @brief Scales a vector by a constant
    /// @param n number of elements in input vector
    /// @param a scalar a
    /// @param x vector which we want to scale
    /// @param incx storage spacing between elements of x
    static void scaleVector(Index n, Real a, Real * x, Index incx); 

    /// @brief This function computes an LU factorization of a real band matrix using partial pivoting with row interchanges
    /// @param nbRows the number of rows
    /// @param nbCols  the number of columns
    /// @param nbSubDiag the number of subdiagonals
    /// @param nbSupDiag the number of superdiagonals
    /// @param LUFactorizedMatrix Matrix contains LU factorization with :  U is stored as an upper triangular band matrix,
    /// and the multipliers used during the factorization are stored below
    /// @param dimLUMatrix the dimension of LUFactorizedMatrix
    /// @param indexPivot the pivot indices    
    static void computeLUBandMatrix(Index nbRows, Index nbCols, Index nbSubDiag, Index nbSupDiag, Matrix& LUFactorizedMatrix, Index dimLUMatrix, Vector& indexPivot);

    /// @brief This function solves the system of equations A*x=B. Where b and x are n elements vectors and A is as n by n band matrix with k+1 diagonals
    /// @param matrixOrder the order of the matrix A
    /// @param nbSupDiag the number of super-diagonals of the matrix A
    /// @param Matrix the matrix A
    /// @param dimMatrix the first dimension of matrix A
    /// @param x the vector x
    /// @param incrementX the increment for the elements of x
    static void solveUxB(Index matrixOrder, Index nbSupDiag, const Matrix& Matrix, Real * x, Index incrementX );

    static void solveAxB(Index matrixOrder, Index nbSubDiag, Index nbSupDiag, Index nbColumnB, Matrix& LUFactorizedMatrix, Vector indexPivot, Vector& B, Index dimB);
};

#if !defined(SOFA_COMPONENT_LINEARSOLVER_BANDMATRIXSOLVER_CPP)
extern template class BANDMATRIX_API BandMatrixSolver< linearalgebra::BandMatrix<SReal>, linearalgebra::FullVector<SReal> >;
#endif
} //namespace sofa::component::linearsolver::direct
