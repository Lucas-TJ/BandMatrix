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
    /// Simple banded matrix container
    template<typename T>
    class BandMatrix : public linearalgebra::BaseMatrix
    {

        public:
            enum { category = MATRIX_BAND };
            typedef T Real;
            typedef typename linearalgebra::BaseMatrix::Index Index;
            ///typedef typename sofa::type::vector vector;
            


        protected:
            sofa::type::vector<sofa::type::vector<T>> data;
            SReal* StoredData;
            Index nbRow, nbCol;
            Index bandwidth;
            Index getBandIndex(Index i, Index j, int bandwidth) const;
            void insideBand(Index i, Index j, int band);
            
        
        public:
            BandMatrix();

            ~BandMatrix() override;

            void resize(Index nbRow, Index nbCol) override;

            Index rowSize(void) const override;

            Index colSize(void) const override;

            SReal element(Index i, Index j) const override;

            void set(Index i, Index j, double v) override;

            void add(Index i, Index j, double v) override;

            void clear() override;
            

           
            
            


        /// Petit test/variante pour apprendre : ici le but est de créer une classe qui calculera automatiquement la taille de la bande
        
            
            

    };


}