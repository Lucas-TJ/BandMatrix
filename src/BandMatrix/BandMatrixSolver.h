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
    SOFA_CLASS(SOFA_TEMPLATE2(BandMatrixSolver,TMatrix,TVector), SOFA_TEMPLATE2(sofa::component::linearsolver::MatrixLinearSolver, TMatrix, TVector));

    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef typename Vector::Real Real;
    typedef sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector> Inherit;

    typedef typename linearalgebra::BaseMatrix::Index Index;

    BandMatrixSolver();

    ///typename Matrix::InvMatrixType Minv;  

    ////////////////////////////


    Vector Y;
protected:
    

    Index incx;
    Vector dx;
    Matrix ab;
    Matrix a;
    Vector x;
    Vector y;

    
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

    Index idamax(Index n, Vector dx, Index incx);

    void dger(Index m, Index n, double alpha, Vector x, Index incx, Vector y, Index incy, Matrix a, Index lda);

    void dswap(Index n, Vector x, Index incx, Vector y, Index incy);

    void dscal(Index n, SReal a, Vector x, Index incx);

    void dgbtf2(Index m, Index n, Index kl, Index ku, Matrix& ab, Index ldab, Vector& ipiv);
    

};

#if !defined(SOFA_COMPONENT_LINEARSOLVER_BANDMATRIXSOLVER_CPP)
extern template class BANDMATRIX_API BandMatrixSolver< linearalgebra::BandMatrix<SReal>, linearalgebra::FullVector<SReal> >;
#endif
} //namespace sofa::component::linearsolver::direct
