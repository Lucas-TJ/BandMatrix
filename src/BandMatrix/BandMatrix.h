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
#include <sofa/linearalgebra/config.h>

#include <sofa/linearalgebra/BaseMatrix.h>
#include <sofa/linearalgebra/BlockFullMatrix.h>
#include <sofa/linearalgebra/FullVector.h>
#include <sofa/linearalgebra/matrix_bloc_traits.h>

namespace sofa::linearalgebra
{
    
    /// Simple banded matrix container (with block ?)
    template<typename T>
    class BandMatrix : public linearalgebra::BaseMatrix
    {

        public:
            enum { category = MATRIX_BAND };
            typedef T Real;
            typedef typename linearalgebra::BaseMatrix::Index Index;



        protected:
            SReal* data;
            Index nbRow, nbCol;
            
        
        public:
            BandMatrix();

            ~BandMatrix() override;

            void resize(Index nbRow, Index nbCol) override;

            Index rowSize(void) const override;

            Index colSize(void) const override;

            SReal element(Index i, Index j) const override;

            void set(Index i, Index j, double v) override
            {
            }

            void add(Index i, Index j, double v) override
            {
            }
            void clear(Index i, Index j);
            void clearRow(Index i);
            void clearCol(Index j);
            void clearRowCol(Index i);
            void clear() ;
            Index getBandWidth() const {return -1;}



            class StoredMatrix : public type::Mat
            {
                public:
                   const T& element(Index i, Index j)
                    {

                    }
                    void set(Index i, Index j, const T& v)
                    {
                        (*this)[i][j] = v;
                    }
                    void add(Index i, Index j, const T& v)
                    {

                    }
                    
            }
            

    };


}