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
#include <sofa/testing/BaseTest.h>

#include <BandMatrix/BandMatrixSolver.h>

#include <sofa/simulation/Node.h>
#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/simpleapi/SimpleApi.h>

#include <sofa/testing/NumericTest.h>

namespace testing
{

/*
TEST(BandMatrixSolver, indexMax)
{
    // required to be able to use EXPECT_MSG_NOEMIT and EXPECT_MSG_EMIT
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;

    using MatrixType = sofa::linearalgebra::FullMatrix<SReal>;

    MatrixType matrix;
    matrix.resize(4, 4);
    matrix.compress();

    for(auto i = 0; i < 4; i++)
    {
        for(auto j = 0 ; j < 4 ; j++)
        {
            matrix.set(i,j,i+j);
        }
    }

    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
    const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();


    EXPECT_EQ(solver->indexMax2(4,matrix,1,0,0),3);
    
}
*/

/// Type-parameterized tests for classes derived from BaseMatrix
template<class T>
class TestBandMatrixSolver : public sofa::testing::NumericTest<typename T::Real>
{
public:
    using Matrix = typename T::Matrix;
    using Real = typename T::Real;
    using Inherit = sofa::testing::NumericTest<typename T::Real>;

    void checkindexMax2()
    {
        using MatrixType = sofa::linearalgebra::FullMatrix<SReal>;

        MatrixType matrix;
        matrix.resize(4, 4);
        matrix.compress();

        for(auto i = 0; i < 4; i++)
        {
            for(auto j = 0 ; j < 4 ; j++)
            {
                matrix.set(i,j,i+j);
            }
        }

        using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
        const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();


        EXPECT_EQ(solver->indexMax2(4,matrix,1,0,0),3);

    }
};

TYPED_TEST_SUITE_P(TestBandMatrixSolver);

TYPED_TEST_P(TestBandMatrixSolver, indexMax2)
{
    this->checkindexMax2();
}



REGISTER_TYPED_TEST_SUITE_P(TestBandMatrixSolver,
                            indexMax2
);


}