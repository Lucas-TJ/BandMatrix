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
#include <sofa/component/linearsolver/direct/config.h>

#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/component/linearsolver/iterative/MatrixLinearSolver.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/linearalgebra/SparseMatrix.h>
#include <sofa/linearalgebra/FullMatrix.h>
#include <BandMatrix.h>
#include <sofa/helper/map.h>

#include <cmath>

namespace sofa::component::linearsolver::direct
{

/// Direct linear solver for banded matrices
template<class TMatrix, class TVector>
class BandMatrixSolver : public sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(BandMatrixSolver,TMatrix,TVector),SOFA_TEMPLATE2(sofa::component::linearsolver::MatrixLinearSolver,TMatrix,TVector));

    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef typename Vector::Real Real;
    typedef sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector> Inherit;

    BandMatrixSolver();

    /// Computes an LU factorization (DGBTF2)
    void computeLU(Matrix& M) override;

    /// Compute x such as Mx=b. M is not used, it must have been factored before using method invert(Matrix& M)
    /// Correspond à DGBTRS pour l'instant
    void solve (Matrix& M, Vector& x, Vector& b) override;

    /// Factors the matrix. Must be done before solving
    void invert(Matrix& M) override;

private :
    linearalgebra::BandMatrix<typename Vector::Real> L;
};

#if !defined(SOFA_COMPONENT_LINEARSOLVER_BANDMATRIXSOLVER_CPP)
extern template class SOFA_COMPONENT_LINEARSOLVER_DIRECT_API BandMatrixSolver< linearalgebra::SparseMatrix<SReal>, linearalgebra::FullVector<SReal> >;
extern template class SOFA_COMPONENT_LINEARSOLVER_DIRECT_API BandMatrixSolver< linearalgebra::FullMatrix<SReal>, linearalgebra::FullVector<SReal> >;
extern template class SOFA_COMPONENT_LINEARSOLVER_DIRECT_API BandMatrixSolver< linearalgebra::BandMatrix<SReal>, linearalgebra::BandMatrix<SReal> >;
#endif

} //namespace sofa::component::linearsolver::direct
