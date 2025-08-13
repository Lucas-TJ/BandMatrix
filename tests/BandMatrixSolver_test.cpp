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
#include <sofa/testing/BaseTest.h>
#include <BandMatrix/config.h>
///#include <BandMatrix/tests/BandMatrixSolver_test.h>
#include <BandMatrix/BandMatrixSolver.h>
#include <gtest/gtest.h>
#include <sofa/testing/NumericTest.h>
#include <sofa/linearalgebra/BaseMatrix.h>

#include <sofa/linearalgebra/FullMatrix.h>



// TEST(BandMatrixSolver, rank1MatrixUpdate)
// {
//     ///rank1MatrixUpdate( nbRows,  nbCols,  alpha,  x,  incrementX,  y,  incrementY,  updatedMatrix)
//     sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
//     using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

//     MatrixType matrix;
//     matrix.resize(4, 4);
//     matrix.compress();

//     for(auto i = 0; i < 4; i++)
//     {
//         for(auto j = 0 ; j < 4 ; j++)
//         {
//             matrix.set(i,j,i+j);
//         }
//     }

//     using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
//     //const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();

//     using VectorType = sofa::linearalgebra::FullVector<SReal>;
//     VectorType x;
//     VectorType y;
//     x.resize(4);
//     y.resize(4);
//     for (auto i = 0; i<4;i++)
//     {
//         x.set(i,i+1);
//         y.set(i,i+5);
//     }

//     Solver::rank1MatrixUpdate(4,4,1,x,1,y,1,matrix);

//     EXPECT_EQ(matrix.element(0,0),5);
//     EXPECT_EQ(matrix.element(1,1),14);
//     EXPECT_EQ(matrix.element(2,2),25);
//     EXPECT_EQ(matrix.element(3,3),38);

// }

TEST(BandMatrixSolver, computeLUBandMatrix)
{
    ///computeLUBandMatrix(nbRows, nbCols, nbSubDiag, nbSupDiag, LUFactorizedMatrix, dimLUMatrix, indexPivot)
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    MatrixType matrix;

    int nbRows = 4, nbCols = 4, nbSubDiag = 1, nbSupDiag = 1;
    int dimLUMatrix = 2*nbSubDiag + nbSupDiag + 1 ;
    matrix.resize(3, 4);
    matrix.compress();

    // First example 
    // matrix.add(0,1,2);
    // matrix.add(0,2,4);
    // matrix.add(0,3,5);
    // matrix.add(1,0,10);
    // matrix.add(1,1,10);
    // matrix.add(1,2,7);
    // matrix.add(1,3,8);
    // matrix.add(2,0,0);
    // matrix.add(2,1,3);
    // matrix.add(2,2,1);
    // matrix.add(2,3,2);
    
    // Second example
    matrix.add(0,1,3);
    matrix.add(0,2,48);
    matrix.add(0,3,1);
    matrix.add(1,0,58);
    matrix.add(1,1,22);
    matrix.add(1,2,95);
    matrix.add(1,3,14);
    matrix.add(2,0,0);
    matrix.add(2,1,3);
    matrix.add(2,2,15);
    matrix.add(2,3,33);
    
    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType indexPivot;
    indexPivot.resize(std::min(nbCols,nbRows));

    Solver::computeLUBandMatrix(nbRows,nbCols,nbSubDiag,nbSupDiag,matrix,dimLUMatrix,indexPivot);
}

TEST(BandMatrixSolver, SolveAxB)
{
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    MatrixType matrix;

    int nbRows = 4, nbCols = 4, nbSubDiag = 1, nbSupDiag = 1;
    int dimLUMatrix = 2*nbSubDiag + nbSupDiag + 1 ;
    matrix.resize(3, 4);
    matrix.compress();

    // First example
    // matrix.add(0,1,2);
    // matrix.add(0,2,4);
    // matrix.add(0,3,5);
    // matrix.add(1,0,10);
    // matrix.add(1,1,10);
    // matrix.add(1,2,7);
    // matrix.add(1,3,8);
    // matrix.add(2,0,0);
    // matrix.add(2,1,3);
    // matrix.add(2,2,1);
    // matrix.add(2,3,2);
    
    // Second example
    matrix.add(0,1,3);
    matrix.add(0,2,48);
    matrix.add(0,3,1);
    matrix.add(1,0,58);
    matrix.add(1,1,22);
    matrix.add(1,2,95);
    matrix.add(1,3,14);
    matrix.add(2,0,0);
    matrix.add(2,1,3);
    matrix.add(2,2,15);
    matrix.add(2,3,33);
    
    
    
    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType indexPivot;
    indexPivot.resize(std::min(nbCols,nbRows));

    std::cout << matrix << std::endl;
    Solver::computeLUBandMatrix(nbRows,nbCols,nbSubDiag,nbSupDiag,matrix,dimLUMatrix,indexPivot);
    
}



