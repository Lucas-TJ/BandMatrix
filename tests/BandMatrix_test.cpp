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
#include <Sofa.LinearAlgebra.Testing/BaseMatrix_test.h>
#include <BandMatrix/BandMatrix.h>



namespace sofa
{
using namespace sofa::linearalgebra::testing;



template<sofa::Index TNbRows, sofa::Index TNbCols, class TReal >
using TestBandMatrixTypes = ::testing::Types<
    TestBaseMatrixTraits< linearalgebra::BandMatrix<TReal>, TNbRows, TNbCols, TReal>
>;

template<class TReal>
using TestBandMatrix9x9Types = TestBandMatrixTypes<9,9,TReal>;

INSTANTIATE_TYPED_TEST_SUITE_P(
    SofaBandLinearSolver_test_matrixformat,
    TestBaseMatrix,
    TestBandMatrix9x9Types<double>
);































/* 
/// faire un generateMatrix basé sur bandmatrix ? car ici c'est basé sur crsMatrix
template<typename TBlock>
void generateMatrix(sofa::linearalgebra::CompressedRowSparseMatrix<TBlock>& matrix,
    sofa::SignedIndex nbRows, sofa::SignedIndex nbCols,
    typename sofa::linearalgebra::CompressedRowSparseMatrix<TBlock>::Real sparsity,
    long seed)
{
    using Real = typename sofa::linearalgebra::CompressedRowSparseMatrix<TBlock>::Real;
    const auto nbNonZero = static_cast<sofa::SignedIndex>(sparsity * static_cast<Real>(nbRows*nbCols));

    sofa::helper::RandomGenerator randomGenerator;
    randomGenerator.initSeed(seed);

    matrix.resize(nbRows, nbCols);

    for (sofa::SignedIndex i = 0; i < nbNonZero; ++i)
    {
        const auto value = static_cast<Real>(sofa::helper::drand(1));
        const auto row = randomGenerator.random<sofa::Index>(0, nbRows);
        const auto col = randomGenerator.random<sofa::Index>(0, nbCols);
        matrix.add(row, col, value);
    }
    matrix.compress();
}
    
TEST(BandMatrix, ShapeMatrix)
{
    ///Qu'est-ce que qu'il me faudrait pour tester add et resize ?... 
    resize doit me donner la forme souhaitée. Par exemple pour bandwidth = 5 et un nbCol=6
    *     *   a02  a13  a24  a35
    *    a01  a12  a23  a34  a45
    a00  a11  a22  a33  a44  a55
    a10  a21  a32  a43  a54   *
    a20  a31  a42  a53   *    *

    add doit donc remplir cette matrice correctement.


    

    

}
*/

}
