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
    ///const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();

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
    for(auto i = 0; i < 4; i++)
    {
        for(auto j = 0 ; j < 4 ; j++)
        {
            std::cout << matrix[4*i+j] << " |";
        }
        std::cout<<"___" <<std::endl;
    }
    EXPECT_EQ(matrix[4*2],0);
    EXPECT_EQ(matrix[4*3],1);

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
    // Solver::swapVector(2,matrix,1,0,1,3,1);
    // EXPECT_EQ(matrix.element(0,1),4);
    // EXPECT_EQ(matrix.element(3,1),1);
}


TEST(BandMatrixSolver, scaleVector)
{
    //scaleVector(Index n, Real a, Real * x, Index incx) 
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
    ///const Solver::SPtr solver = sofa::core::objectmodel::New<Solver>();

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

// TEST(BandMatrixSolver, computeLUBandMatrix)
// {
//     ///computeLUBandMatrix(nbRows, nbCols, nbSubDiag, nbSupDiag, LUFactorizedMatrix, dimLUMatrix, indexPivot)
//     sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
//     using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
//     MatrixType matrix;

//     int nbRows = 4, nbCols = 4, nbSubDiag = 1, nbSupDiag = 1;
//     int dimLUMatrix = 2*nbSubDiag + nbSupDiag + 1 ;
//     matrix.resize(3, 4);
//     matrix.compress();

//     // First example 
//     // matrix.add(0,1,2);
//     // matrix.add(0,2,4);
//     // matrix.add(0,3,5);
//     // matrix.add(1,0,10);
//     // matrix.add(1,1,10);
//     // matrix.add(1,2,7);
//     // matrix.add(1,3,8);
//     // matrix.add(2,0,0);
//     // matrix.add(2,1,3);
//     // matrix.add(2,2,1);
//     // matrix.add(2,3,2);
    
//     // Second example
//     matrix.add(0,1,3);
//     matrix.add(0,2,48);
//     matrix.add(0,3,1);
//     matrix.add(1,0,58);
//     matrix.add(1,1,22);
//     matrix.add(1,2,95);
//     matrix.add(1,3,14);
//     matrix.add(2,0,0);
//     matrix.add(2,1,3);
//     matrix.add(2,2,15);
//     matrix.add(2,3,33);
    
//     using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

//     using VectorType = sofa::linearalgebra::FullVector<SReal>;
//     VectorType indexPivot;
//     indexPivot.resize(std::min(nbCols,nbRows));

//     std::cout<< "Voici la matrice avant : " << std::endl;
//     std::cout << matrix << std::endl;
//     std::cout<< "Voici la ipiv avant : " << std::endl;
//     std::cout << indexPivot << std::endl;
//     Solver::computeLUBandMatrix(nbRows,nbCols,nbSubDiag,nbSupDiag,matrix,dimLUMatrix,indexPivot);

//     std::cout<< "Voici la matrice après : " << std::endl;
//     std::cout << matrix << std::endl;
//     std::cout<< "Voici la ipiv après : " << std::endl;
//     std::cout << indexPivot << std::endl;

// }

TEST(BandMatrixSolver, SolveUxB)
{
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    MatrixType matrix;
    matrix.resize(4, 4);
    matrix.compress();

    // Remplissage de la matrice
    matrix.add(0,3,12);
    matrix.add(1,2,4);
    matrix.add(1,3,5);
    matrix.add(2,1,1);
    matrix.add(2,2,7);
    matrix.add(2,3,5);
    matrix.add(3,0,3);
    matrix.add(3,1,89);
    matrix.add(3,2,2);
    matrix.add(3,3,8);
    std::cout<< "Voici la matrice : " << std::endl;
    std::cout << matrix << std::endl;

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType vectorx;
    vectorx.resize(4);
    vectorx.add(0,12);
    vectorx.add(1,3);
    vectorx.add(2,5);
    vectorx.add(3,47);
    std::cout<< "Voici le second membre : " << std::endl;
    std::cout<< vectorx << std::endl;
    
    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;
    Solver::solveUxB(4,3,matrix,vectorx.ptr(),1);

    std::cout<< "Voici le resultat : " << std::endl;
    std::cout<< vectorx << std::endl;
    
}

TEST(BandMatrixSolver, SolveAxB)
{
        ///computeLUBandMatrix(nbRows, nbCols, nbSubDiag, nbSupDiag, LUFactorizedMatrix, dimLUMatrix, indexPivot)
    sofa::helper::logging::MessageDispatcher::addHandler(sofa::testing::MainGtestMessageHandler::getInstance() ) ;
    using MatrixType = sofa::linearalgebra::BandMatrix<SReal>;
    MatrixType matrix;

    int nbRows = 4, nbCols = 4, nbSubDiag = 1, nbSupDiag = 1;
    int dimLUMatrix = 2*nbSubDiag + nbSupDiag + 1 ;
    matrix.resize(4, 4);
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
    // matrix.add(0,1,3);
    // matrix.add(0,2,48);
    // matrix.add(0,3,1);
    // matrix.add(1,0,58);
    // matrix.add(1,1,22);
    // matrix.add(1,2,95);
    // matrix.add(1,3,14);
    // matrix.add(2,0,0);
    // matrix.add(2,1,3);
    // matrix.add(2,2,15);
    // matrix.add(2,3,33);

    // Third example (thx gpt)
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
    

    std::cout <<"Voici matrix: " << std::endl;
    std::cout << matrix << std::endl;

    using Solver = sofa::component::linearsolver::direct::BandMatrixSolver<MatrixType, sofa::linearalgebra::FullVector<SReal> >;

    using VectorType = sofa::linearalgebra::FullVector<SReal>;
    VectorType indexPivot;
    indexPivot.resize(std::min(nbCols,nbRows));

    Solver::computeLUBandMatrix(nbRows,nbCols,nbSubDiag,nbSupDiag,matrix,dimLUMatrix,indexPivot);

    std::cout <<"Voici matrix LU : " << std::endl;
    std::cout << matrix << std::endl;

    int matrixOrder = matrix.rowSize();
    int nbColumnB = 1;
    int dimB = 4;
    MatrixType B;
    B.resize(4,nbColumnB);
    B.set(0,0,5);
    B.set(1,0,10);
    B.set(2,0,15);
    B.set(3,0,10);
    
    std::cout << "Voici B" << std::endl;
    std::cout << B << std::endl;

    // using VectorType = sofa::linearalgebra::FullVector<SReal>;
    // VectorType B;
    // // sofa::type::vector<double> B;
    // B.resize(4);
    // // B.add(0,5);
    // // B.add(1,10);
    // // B.add(2,15);
    // // B.add(3,10);

    // B[0] = 5;
    // B[1] = 10;
    // B[2] = 15;
    // B[3] = 10 ;

    Solver::solveAxB(matrixOrder, nbSubDiag, nbSupDiag, nbColumnB, matrix, indexPivot, B, dimB);

    std::cout << "VOICI LE RESULTAT FINAL : " << std::endl;
    std::cout << B << std::endl;

    std::cout << "ipiv = " << std::endl;
    std::cout << indexPivot << std::endl;
    
}



