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

TEST(BandMatrixSolver,indexMax)
{
    /*
    This unit test is about the function indexMax2. indexMax2 returns the index from a column in a matrix which has the maximum value.
    We construct a matrix like this :
        0 | 1 | 2 | 3
        1 | 2 | 3 | 4
        2 | 3 | 4 | 5
        3 | 4 | 5 | 6

    Then, we test indexMax with different parameter (nbElem, Matrix, increment, currentColumn, DiagP)
    */
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

    sofa::type::vector<double> matrix;
    matrix.resize(4*4);


    for(auto i = 0; i < 4; i++)
    {
        for(auto j = 0 ; j < 4 ; j++)
        {
            matrix[4*i+j]= i+j;
        }
    }
    double * ptr = matrix.data();


    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

    /// Test when increment == 1 : 
    EXPECT_EQ(Solver::indexMax(4,ptr,1),3);
    EXPECT_EQ(Solver::indexMax(4,ptr+4,1),3);
    EXPECT_EQ(Solver::indexMax(4,ptr+4*2,1),3);
    /// Test when increment != 1 :
    EXPECT_EQ(Solver::indexMax(2,ptr,2),1);
};

TEST(BandMatrixSolver, swapVector)
{
    /*
    This unit test is about the function swapVector. SwapVector interchanges two vectors from a place to an other using parameters "indexOfM1" abd "indexOfM2"
    We construct a matrix like this :
        0 | 1 | 2 | 3
        1 | 2 | 3 | 4
        2 | 3 | 4 | 5
        3 | 4 | 5 | 6

    swapVector is defined like this : (cf .h for documentation) 
    swapVector2(Index nbElemToSwap, Matrix& M, Index currentColumn, Index indexOfM1, Index incrementM1, Index indexOfM2, Index incrementM2 )
    static void swapVector(Index n, Real * x, Index incx, Real * y, Index incy);
    */
    

    ///
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

    sofa::type::vector<double> matrix;
    matrix.resize(4*4);

    for(auto i = 0; i < 4; i++)
    {
        for(auto j = 0 ; j < 4 ; j++)
        {
            matrix[4*i + j]= i+j;
        }
    }


    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
    //const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();
    double * ptr = matrix.data();
    
    /*Then, we swap the two first rows with the last two in the first column 
    Before :
        0 | 1 | 2 | 3
        1 | 2 | 3 | 4
        2 | 3 | 4 | 5
        3 | 4 | 5 | 6
    After : 
        2 | 1 | 2 | 3
        3 | 2 | 3 | 4
        0 | 3 | 4 | 5
        1 | 4 | 5 | 6

    */
    Solver::swapVector(2,ptr,4,ptr+4*2,4);
    EXPECT_EQ(matrix[4*2],0);
    EXPECT_EQ(matrix[4*3],1);

}


TEST(BandMatrixSolver, scaleVector)
{
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    sofa::type::vector<double> matrix;
    matrix.resize(4*4);

    for(auto i = 0; i < 4; i++)
    {
        for(auto j = 0 ; j < 4 ; j++)
        {
            matrix[4*i + j]= i+j;
        }
    }

    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

    double * ptr = matrix.data();
    Solver::scaleVector(4,2,ptr,1);
    Solver::scaleVector(4,3,ptr + 4,2);
    Solver::scaleVector(4,4,ptr + 4*2 + 1,3);

    EXPECT_EQ(matrix[0],0);
    EXPECT_EQ(matrix[1],2);
    EXPECT_EQ(matrix[2],4);
    EXPECT_EQ(matrix[3],6);
    
    EXPECT_EQ(matrix[4],3);
    EXPECT_EQ(matrix[4+2],9);

    EXPECT_EQ(matrix[4*2+1],12);
    EXPECT_EQ(matrix[4*2+1+3],12);
    
}

TEST(BandMatrixSolver, rank1MatrixUpdate)
{
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

    sofa::type::vector<double>  matrix;
    matrix.resize(4* 4);
    

    for(auto i = 0; i < 4; i++)
    {
        for(auto j = 0 ; j < 4 ; j++)
        {
            matrix[i*4+j]=i+j;
        }
    }
    SReal * ptr = matrix.data();
    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType x;
    VectorType y;
    x.resize(4);
    y.resize(4);
    for (auto i = 0; i<4;i++)
    {
        x.set(i,i+1);
        y.set(i,i+5);
    }

    Solver::rank1MatrixUpdate(4,4,1,x.ptr(),1,y.ptr(),1,ptr);

    EXPECT_EQ(matrix[0],5);
    EXPECT_EQ(matrix[1*4+1],14);
    EXPECT_EQ(matrix[2*4+2],25);
    EXPECT_EQ(matrix[3*4+3],38);

}

TEST(BandMatrixSolver, computeLUBandMatrix)
{
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    MatrixType matrix;

    int nbRows = 4, nbCols = 4, nbSubDiag = 1, nbSupDiag = 1;
    int dimLUMatrix = 2*nbSubDiag + nbSupDiag + 1 ;
    matrix.resize(4, 4);
    matrix.compress();

    // Example 
    matrix.add(1,1,3);
    matrix.add(1,2,48);
    matrix.add(1,3,1);
    matrix.add(2,0,58);
    matrix.add(2,1,22);
    matrix.add(2,2,95);
    matrix.add(2,3,14);
    matrix.add(3,0,0);
    matrix.add(3,1,3);
    matrix.add(3,2,15);
    matrix.add(3,3,33);
    
    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType indexPivot;
    indexPivot.resize(std::min(nbCols,nbRows));

    Solver::computeLUBandMatrix(nbRows,nbCols,nbSubDiag,nbSupDiag,matrix,dimLUMatrix,indexPivot);
    
    EXPECT_EQ(matrix.element(2,2), 88.454545454545453);
    EXPECT_EQ(matrix.element(2,3), 13.830421377183967);


    EXPECT_EQ(matrix.element(3,1), 0.13636363636363635);
    EXPECT_EQ(matrix.element(3,2), 0.16957862281603289);



}

TEST(BandMatrixSolver, SolveUxB)
{
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    sofa::type::vector<double>  matrix;
    matrix.resize(4* 4);
    
    // Remplissage de la matrice
    matrix[0+ 3] = 12;
    matrix[1*4+ 2]=4;
    matrix[1*4+ 3]=5;
    matrix[2*4+1]=1;
    matrix[2*4+2]=7;
    matrix[2*4+3]=5;
    matrix[3*4+0]=3;
    matrix[3*4+1]=89;
    matrix[3*4+2]=2;
    matrix[3*4+3]=8;

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType vectorx;
    vectorx.resize(4);
    vectorx.add(0,12);
    vectorx.add(1,3);
    vectorx.add(2,5);
    vectorx.add(3,47);

    double * ptr = matrix.data();
    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
    Solver::solveUxB(4,3,ptr,vectorx.ptr(),1,4);

    EXPECT_EQ(vectorx[0], -3.4707397003745317);
    EXPECT_EQ(vectorx[1], 0.6622191011235955);
    EXPECT_EQ(vectorx[2], -12.1875);
    EXPECT_EQ(vectorx[3], 5.875);

    
}

TEST(BandMatrixSolver, SolveAxB)
{
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    MatrixType matrix;

    int nbRows = 4, nbCols = 4, nbSubDiag = 1, nbSupDiag = 1;
    int dimLUMatrix = 2*nbSubDiag + nbSupDiag + 1 ;
    matrix.resize(4, 4);
    matrix.compress();
    
    // example 
    matrix.add(1,1,2);
    matrix.add(1,2,3);
    matrix.add(1,3,2);
    
    matrix.add(2,0,4);
    matrix.add(2,1,4);
    matrix.add(2,2,4);
    matrix.add(2,3,4);

    matrix.add(3,0,1);
    matrix.add(3,1,1);
    matrix.add(3,2,1);
    

    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType indexPivot;
    indexPivot.resize(std::min(nbCols,nbRows));

    Solver::computeLUBandMatrix(nbRows,nbCols,nbSubDiag,nbSupDiag,matrix,dimLUMatrix,indexPivot);

    EXPECT_EQ(matrix.element(2,1), 3.5);
    EXPECT_EQ(matrix.element(2,2), 3.1428571428571428);
    EXPECT_EQ(matrix.element(2,3), 3.3636363636363638);

    EXPECT_EQ(matrix.element(3,0), 0.25);
    EXPECT_EQ(matrix.element(3,1), 0.2857142857142857);
    EXPECT_EQ(matrix.element(3,2), 0.31818181818181818);
    EXPECT_EQ(matrix.element(3,3), 0);

    int matrixOrder = matrix.rowSize();
    int nbColumnB = 1;
    int dimB = 4;
    MatrixType B;
    B.resize(4,nbColumnB);
    B.set(0,0,5);
    B.set(1,0,10);
    B.set(2,0,15);
    B.set(3,0,10);

    Solver::solveAxB(matrixOrder, nbSubDiag, nbSupDiag, nbColumnB, matrix, indexPivot, B, dimB);

    EXPECT_EQ(B.element(0,0), 1.2162162162162162);
    EXPECT_EQ(B.element(1,0), 0.067567567567567419);
    EXPECT_EQ(B.element(2,0), 2.8378378378378382);
    EXPECT_EQ(B.element(3,0), 1.7905405405405406);
    
}



