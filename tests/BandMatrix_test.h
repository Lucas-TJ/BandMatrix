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
#include <gtest/gtest.h>
#include <sofa/type/Mat.h>
#include <sofa/linearalgebra/config.h>
#include <sofa/testing/NumericTest.h>
#include <sofa/linearalgebra/BaseMatrix.h>
#include <BandMatrix/BandMatrix.h>
#include <sofa/type/vector.h>
#include <map>




namespace sofa::linearalgebra::testing
{

template<class TMatrix, sofa::Index TNbRows, sofa::Index TNbCols, class TReal >
struct TestBandMatrixTraits
{
    using Matrix = TMatrix;
    static constexpr sofa::Index NbRows = TNbRows;
    static constexpr sofa::Index NbCols = TNbCols;
    using Real = TReal;
};

template<class T>
class TestBandMatrix : public sofa::testing::NumericTest<typename T::Real>
{
    
public:
    using Matrix = typename T::Matrix;  
    using Real = typename T::Real;
    using Inherit = sofa::testing::NumericTest<typename T::Real>;

    static constexpr sofa::Index NbRows = T::NbRows;
    static constexpr sofa::Index NbCols = T::NbCols;


    void doSetUp() override
    {
        m_testedMatrix = std::make_unique<Matrix>();
        m_testedMatrix->resize(NbRows, NbCols);

        m_modelMatrix.clear();
    }

    void doTearDown() override
    {
        m_testedMatrix.reset();
    }

  
    void checkAdd()
    {
        m_testedMatrix->clear();

        Real value = (Real)0;
        for (sofa::Index i = 0 ; i < NbRows; ++i)
        {
            for (sofa::Index j = 0 ; j < NbCols; ++j)
            {
                m_modelMatrix(i, j) = ++value;
                m_testedMatrix->add(i, j, value);
            }
        }

        m_testedMatrix->compress(); //this operation is not required for all types of BaseMatrix

        EXPECT_LT(Inherit::matrixMaxDiff(m_modelMatrix, *m_testedMatrix), 100 * Inherit::epsilon());
        
    }


protected:

    sofa::type::Mat<NbRows, NbCols, Real> m_modelMatrix;
    std::unique_ptr<sofa::linearalgebra::BaseMatrix> m_testedMatrix {nullptr};



};

TYPED_TEST_SUITE_P(TestBandMatrix);


TYPED_TEST_P(TestBandMatrix, add)
{
    this->checkAdd();
}

REGISTER_TYPED_TEST_SUITE_P(TestBandMatrix,
                            add
);


} /// namespace sofa::linearalgebra::testing
