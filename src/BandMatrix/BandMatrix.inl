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
#include <BandMatrix/BandMatrix.h>

namespace sofa::linearalgebra
{

    template<typename T>
    BandMatrix<T>::BandMatrix()
        : data(), nbRow(0), nbCol(0), bandwidth(0)
    {
    }

    template<typename T>
    BandMatrix<T>::~BandMatrix()
    {
        
    }

    
    template<typename T>
    void BandMatrix<T>::resize(Index nr, Index nc)
    {
        nbRow = nr;
        nbCol = nc;
        data.resize(bandwidth);
        for (Index i = 0; i < bandwidth; ++i)
        {
            data[i].resize(nbCol, 0);
            /// data[i].assign(nbCol, 0);
        }
        
    }

    template<typename T>
    void BandMatrix<T>::insideBand(Index i, Index j, int banded_i)
    {     
        if (( banded_i < 0 ) || (banded_i >= bandwidth) )
        {
            Index ki; 
            const int oldbandwidth = bandwidth;
            bandwidth = std::max(2*i + 1, 2*j + 1); 

            int nbAdd; 
            nbAdd = (bandwidth - oldbandwidth)/2;
            
            resize(nbRow,nbCol);
            
            for (ki = oldbandwidth-1; ki>-1; ki--)
            {
                data[ki+nbAdd] = data[ki];
            }
            
            for (ki = 0; ki < nbAdd; ki++)
            {
                data[ki].clear();
                data[ki].resize(nbCol);
                
            }

            
        }
    }

    template<typename T>
    void BandMatrix<T>::add(Index i, Index j, double v)
    {
        int banded_i = getBandIndex(i,j,bandwidth);
        ///std::cout << "Before insideBand: i=" << i << " j=" << j << " bandwidth=" << bandwidth << " banded_i=" << banded_i << std::endl;
        insideBand(i,j,banded_i);
        ///std::cout << "After insideBand: i=" << i << " j=" << j << " bandwidth=" << bandwidth << std::endl;
        banded_i = getBandIndex(i,j,bandwidth);
        ///std::cout << "After recompute: banded_i=" << banded_i << std::endl;
        data[banded_i][j] += v;
        ///std::cout << "===============> data[b_i][j] += v, banded_i =  " << banded_i << ", j= " << j << ", v = " << v << std::endl;
        
    } 

    template<typename T>
    SReal BandMatrix<T>::element(Index i, Index j) const
    {
        Index banded_i = getBandIndex(i,j,bandwidth);
        if (( banded_i >= 0 ) && (banded_i < bandwidth)) return (SReal)data[banded_i][j];
        else return 0;
            
    }

    template<typename T>
    typename BandMatrix<T>::Index BandMatrix<T>::getBandIndex(Index i, Index j, int bandwidth) const
    {
        return (bandwidth-1)/2 + i - j;
    }

    template<typename T>
    void BandMatrix<T>::set(Index i, Index j, double v)
    {
        Index banded_i = getBandIndex(i,j,bandwidth);
        insideBand(i,j,banded_i);
        data[banded_i][j] = v;
        
    } 

    template<typename T>
    typename BandMatrix<T>::Vector BandMatrix<T>::getVector(Index i, Index j)
    {
        Vector dataVec;
        dataVec.resize(nbRow);
        for (Index ki = i; ki<nbRow;i++)
        {

            dataVec[ki] = this->element(ki,j);
        }

        return dataVec;
    }

    template<typename T>
    typename BandMatrix<T>::Matrix BandMatrix<T>::getMatrix()
    {
        int rows = rowSize();
        int cols = colSize();
        for(Index j = 0; j < cols; j++)
        {
            for(Index i = 0; i < rows; i++)
            {
                matrix[i + j * rows]=data[i][j];
            }
        }
        return matrix;
    }

    template<typename T>
    typename BandMatrix<T>::Index BandMatrix<T>::rowSize() const
    {
        return nbRow;
    }
    template<typename T>
    typename BandMatrix<T>::Index BandMatrix<T>::colSize() const
    {
        return nbCol;
    }
  
    template<typename T>
    void BandMatrix<T>::clear()
    {
        // nbRow = 0;
        // nbCol = 0;
        // bandwidth = 0;
        // data.clear();

        for (auto& d : data)
            for (auto& e : d)
                e = 0;
            
    }
    

    
}