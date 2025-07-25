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
///#include <Sofa.LinearAlgebra.Testing/BaseMatrix_test.h>
#include <BandMatrix/tests/BandMatrixSolver_test.h>
#include <BandMatrix/BandMatrixSolver.h>



namespace sofa
{
using namespace testing;



template<sofa::Index TNbRows, sofa::Index TNbCols, class TReal >
using TestBandMatrixTypes = ::testing::Types<
    TestBaseMatrixTraits< linearalgebra::BandMatrix<TReal>, TNbRows, TNbCols, TReal>
>;

template<class TReal>
using TestBandMatrix9x9Types = TestBandMatrixTypes<9,9,TReal>;

INSTANTIATE_TYPED_TEST_SUITE_P(
    SofaBandLinearSolver_test_indexmax,
    TestBaseMatrix,
    TestBandMatrix9x9Types<double>
);

}




