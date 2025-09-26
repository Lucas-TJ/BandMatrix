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
#include <BandMatrix/config.h>
#include <sofa/type/Mat.h>
#include <sofa/linearalgebra/BaseMatrix.h>
#include <sofa/linearalgebra/BlockFullMatrix.h>
#include <sofa/linearalgebra/FullVector.h>
#include <sofa/linearalgebra/matrix_bloc_traits.h>
#include <sofa/type/vector.h>

#include <map>

namespace sofa::linearalgebra
{
    using sofa::type::vector;
    using namespace sofa::linearalgebra;
    /// Simple banded matrix container
    template<typename T>
    class BandMatrix : public linearalgebra::BaseMatrix
    {

        public:
            enum { category = MATRIX_BAND };
            typedef T Real;
            typedef FullVector<Real> Vector;
            using linearalgebra::BaseMatrix::Index;
            


        protected:
            sofa::type::vector<sofa::type::vector<T>> data;
            Index nbRow, nbCol;
            Index bandwidth;

            /// Get the index corresponding to the place in the band
            Index getBandIndex(Index i, Index j, int bandwidth) const;

            /// Check if the data is inside the band. 
            /// If the data is not inside the band, the bandwidth will adapt his width in order to fit with the data
            void insideBand(Index i, Index j, int band);
            
        
        public:
            BandMatrix();

            ~BandMatrix() override;

            /// Resize the matrix and reset all values to 0
            void resize(Index nbRow, Index nbCol) override;

            /// Number of rows
            Index rowSize(void) const override;

            /// Number of columns
            Index colSize(void) const override;

            /// Read the value of the element at row i, column j (using 0-based indices)
            SReal element(Index i, Index j) const override;

            /// Write the value of the element at row i, column j (using 0-based indices)
            void set(Index i, Index j, double v) override;

            /// Add v to the existing value of the element at row i, column j (using 0-based indices)
            void add(Index i, Index j, double v) override;

            /// Reset all values to 0
            void clear() override;
            
            Vector getVector(Index i, Index j);

            static const char* Name()
            {
                static std::string name { "BandMatrix" };
                return name.c_str();
            }

    };

#if !defined(BANDMATRIX_CPP)
extern template class BANDMATRIX_API linearalgebra::BandMatrix<SReal>;
#endif

}