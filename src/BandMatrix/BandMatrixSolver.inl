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
#include <BandMatrixSolver.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/linearalgebra/SparseMatrix.h>
#include <BandMatrix.h>
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



template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::computeLU(Matrix& M)
{
    Index mu;
    Index ml;
    getBandLU(mu,ml);
    Index mv ;
    mv = mu + ml ; /// mu et ml sont determinés dans BandMatrix.inl avec la fonction getBandLU
    for (Index j = mu + 2; std::min(mv,N))
    {
        for (Index i = mv- j + 2; ml)
        {

        }
    }

    for (Index )
    
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::solve (Matrix& /*M*/, Vector& z, Vector& r)
{
    this->mat
}

template<class TMatrix, class TVector>
void BandMatrixSolver<TMatrix,TVector>::invert(Matrix& M)
{
    
}

} //namespace sofa::component::linearsolver::direct