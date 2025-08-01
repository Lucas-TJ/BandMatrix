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

TEST(BandMatrixSolver,indexMax2)
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
    ///const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();

    /// Test when increment == 1 : 
    EXPECT_EQ(Solver::indexMax2(4,matrix,1,0,0),3);
    EXPECT_EQ(Solver::indexMax2(2,matrix,1,0,0),1);
    EXPECT_EQ(Solver::indexMax2(2,matrix,1,2,1),1);
    EXPECT_EQ(Solver::indexMax2(3,matrix,1,1,0),2);    
    EXPECT_EQ(Solver::indexMax2(2,matrix,1,3,2),1);
    /// Test when increment != 1 :
    EXPECT_EQ(Solver::indexMax2(2,matrix,2,0,1),1);
    EXPECT_EQ(Solver::indexMax2(2,matrix,2,0,0),1);
    EXPECT_EQ(Solver::indexMax2(2,matrix,2,2,1),1);
};

TEST(BandMatrixSolver, swapVector2)
{
    /*
    This unit test is about the function swapVector2. SwapVector2 interchanges two vectors from a place to an other using parameters "indexOfM1" abd "indexOfM2"
    We construct a matrix like this :
        0 | 1 | 2 | 3
        1 | 2 | 3 | 4
        2 | 3 | 4 | 5
        3 | 4 | 5 | 6

    swapVector2 is defined like this : (cf .h for documentation) 
    swapVector2(Index nbElemToSwap, Matrix& M, Index currentColumn, Index indexOfM1, Index incrementM1, Index indexOfM2, Index incrementM2 )
    */
    

    ///
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

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
    //const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();
    
    
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
    Solver::swapVector2(4,matrix,0,0,1,2,1);
    EXPECT_EQ(matrix.element(2,0),0);
    EXPECT_EQ(matrix.element(3,0),1);

    /*And, we just swap two elements (the first && the last one) in the second column 
    Before :
        0 | 1 | 2 | 3
        1 | 2 | 3 | 4
        2 | 3 | 4 | 5
        3 | 4 | 5 | 6
    After : 
        0 | 4 | 2 | 3
        1 | 2 | 3 | 4
        2 | 3 | 4 | 5
        3 | 1 | 5 | 6

    */
    Solver::swapVector2(2,matrix,1,0,1,3,1);
    EXPECT_EQ(matrix.element(0,1),4);
    EXPECT_EQ(matrix.element(3,1),1);

}

TEST(BandMatrixSolver, rank1MatrixUpdate)
{
    ///rank1MatrixUpdate( nbRows,  nbCols,  alpha,  x,  incrementX,  y,  incrementY,  updatedMatrix)
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

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
    //const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();

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

    Solver::rank1MatrixUpdate(4,4,1,x,1,y,1,matrix);

    EXPECT_EQ(matrix.element(0,0),5);
    EXPECT_EQ(matrix.element(1,1),14);
    EXPECT_EQ(matrix.element(2,2),25);
    EXPECT_EQ(matrix.element(3,3),38);

}

// TEST(BandMatrixSolver,scaleVector2)
// {
//     ///scaleVector2(Index nbElem, SReal a, Matrix& M, Index incx, Index currentColumn) 
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
//     ///const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();

//     // Solver::scaleVector2(4,2,matrix,1,0,0);

//     // EXPECT_EQ(matrix.element(0,0),0);
//     // EXPECT_EQ(matrix.element(1,0),2);
//     // EXPECT_EQ(matrix.element(2,0),4);
//     // EXPECT_EQ(matrix.element(3,0),6);

//     // Solver::scaleVector2(2,2,matrix,1,3,0);
//     Solver::scaleVector2(2,2,matrix,1,3,2);
//     EXPECT_EQ(matrix.element(0,3),6);
//     EXPECT_EQ(matrix.element(1,3),8);
//     EXPECT_EQ(matrix.element(2,3),10);
//     EXPECT_EQ(matrix.element(3,3),12);
//     std::cout<<matrix<<std::endl;

// }

TEST(BandMatrixSolver,scaleVector3)
{
    //Index nbElem, SReal a, Matrix& M, Index incx, Index currentColumn, Index indexDiagP
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

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
    ///const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();


    Solver::scaleVector3(4,2,matrix,1,0,0);
    EXPECT_EQ(matrix.element(0,0),0);
    EXPECT_EQ(matrix.element(1,0),2);
    EXPECT_EQ(matrix.element(2,0),4);
    EXPECT_EQ(matrix.element(3,0),6);

    
    Solver::scaleVector3(2,2,matrix,2,2,1);
    EXPECT_EQ(matrix.element(1,2),6);
    EXPECT_EQ(matrix.element(3,2),10);
    
    Solver::scaleVector3(2,2,matrix,1,3,0);
    EXPECT_EQ(matrix.element(0,3),6);
    EXPECT_EQ(matrix.element(1,3),8);


}

TEST(BandMatrixSolver, rank1MatrixUpdate2)
{
    ///rank1MatrixUpdate2( nbRows,  nbCols,  alpha,  matrix,  rowOfM1,  colOfM1 , incrementX,  rowOfM2,  colOfM2,  incrementY, updatedMatrix);
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;

    MatrixType matrix;
    matrix.resize(6, 6);
    matrix.compress();

    for(auto i = 0; i < 6; i++)
    {
        for(auto j = 0 ; j < 6 ; j++)
        {
            matrix.set(i,j,i+j);
        }
    }

    ///std::cout << matrix << std::endl;

    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
    //const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();

    // using VectorType = sofa::linearalgebra::FullVector<SReal>;
    // VectorType x ;
    // VectorType y ;

    // x.resize(6);
    // y.resize(6);
    // for(auto i = 0; i<6;i++)
    // {
    //     for(auto j = 0; j<6;j++)
    //     {
    //         auto tempx = matrix.element(i,0);
    //         auto tempy = matrix.element(i,1);
    //         x.set(i,tempx);
    //         y.set(i,tempy);
    //     }
    // }

    // std::cout << "vecteur x = " << x << std::endl;
    // std::cout << "vecteur y = " << y << std::endl;

    Solver::rank1MatrixUpdate2(6,6,2,matrix,0,0,1,0,1,1,matrix);

    //std::cout << matrix << std::endl;
    EXPECT_EQ(matrix.element(0,0),0);
    EXPECT_EQ(matrix.element(1,1),6);
    EXPECT_EQ(matrix.element(2,2),16);
    EXPECT_EQ(matrix.element(3,3),30);
    EXPECT_EQ(matrix.element(4,4),48);
    EXPECT_EQ(matrix.element(5,5),70);

}


TEST(BandMatrixSolver, computeLUBandMatrix) /// MAINTENANT
{
    ///computeLUBandMatrix(nbRows, nbCols, nbSubDiag, nbSupDiag, LUFactorizedMatrix, dimLUMatrix, indexPivot)
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    MatrixType matrix;

    int nbRows = 5, nbCols = 5, nbSubDiag = 2, nbSupDiag = 2;
    ///std::cout<<"nbSupDiag = " << nbSupDiag << std::endl;
    int dimLUMatrix = 2*nbSubDiag + nbSupDiag + 1 ;
    ///std::cout<<"dimLUMatrix = " << dimLUMatrix << std::endl;
    matrix.resize(nbRows, nbCols);
    matrix.compress();

    for(auto i = 0; i < dimLUMatrix; i++)
    {
        for(auto j = 0 ; j < nbRows ; j++)
        {
            if ((i==2) || ((i==3)&&(j<4)) || ((i==1)&&(j>0)) || ((i==0)&&(j>1)) || ((i==4)&&(j<3)))
            matrix.add(i,j,i+j);
        }
    }
    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
    //const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();
    //solver->init();

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType indexPivot;
    indexPivot.resize(std::min(matrix.colSize(),matrix.rowSize()));

    

    //std::cout << matrix << std::endl;
    Solver::computeLUBandMatrix(nbRows,nbCols,nbSubDiag,nbSupDiag,matrix,dimLUMatrix,indexPivot);
    std::cout << matrix << std::endl;
    std::cout << indexPivot << std::endl;
}






