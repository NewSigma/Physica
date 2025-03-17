/*
 * Copyright 2024-2025 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/SparseVector.h"
#include "Physica/Core/Utils/Container/ArrayND.h"

namespace Physica {
    /**
     * Force constant(FC) swap vector is the constraint vector $w$ as introduced in [1].
     * 
     * References:
     * [1] Mounet, N. Structural, vibrational and thermodynamic properties of carbon allotropes from ﬁrst-principles: diamond, graphite, and nanotubes. Master’s thesis, Massachusetts Institute of Technology (2005); http://hdl.handle.net/1721.1/33400
     */
    template<Scalar T>
    class FCSwapVector : public SparseVector<T> {
        using Base = SparseVector<T>;

        Index3D superSize;
        Index3D cellIndex;
        size_t numDOF;
        size_t dof1;
        size_t dof2;
    public:
        FCSwapVector(Index3D superSize_, Index3D cellIndex_, size_t numDOF_, size_t dof1_, size_t dof2_);
        FCSwapVector(const FCSwapVector&) = default;
        FCSwapVector(FCSwapVector&&) noexcept = default;
        ~FCSwapVector() = default;
        /* Operators */
        FCSwapVector& operator=(FCSwapVector obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(FCSwapVector& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] static Index5D index1DTo5D(size_t numDOF, const Index3D& superSize, size_t index1D);
        [[nodiscard]] static size_t index5DTo1D(size_t numDOF, const Index3D& superSize, const Index5D& index5D);
    private:
        [[nodiscard]] inline static size_t calcLength(size_t numDOF, const Index3D& superSize);
        inline static void checkParam(size_t numDOF, const Index3D& superSize, const Index5D& index);
    };

    template<Scalar T>
    FCSwapVector<T>::FCSwapVector(Index3D superSize_, Index3D cellIndex_, size_t numDOF_, size_t dof1_, size_t dof2_)
            : superSize(std::move(superSize_))
            , cellIndex(std::move(cellIndex_))
            , numDOF(numDOF_)
            , dof1(dof1_)
            , dof2(dof2_) {
        assert(!(cellIndex[0] == 0 && cellIndex[1] == 0 && cellIndex[2] == 0 && dof1 == dof2)
              && "[Error]: Vector will be empty using the given params");
        checkParam(numDOF, superSize, Index5D{cellIndex[0], cellIndex[1], cellIndex[2], dof1, dof2});
        Base::resize(calcLength(numDOF, superSize));
        Base::reserve(2);

        Index5D temp{};
        for (int i = 0; i < 3; ++i)
            temp[i] = cellIndex[i];
        temp[3] = dof1;
        temp[4] = dof2;
        Base::operator[](index5DTo1D(numDOF, superSize, temp)) = T(M_SQRT1_2);

        for (int i = 0; i < 3; ++i)
            temp[i] = (superSize[i] - cellIndex[i]) % superSize[i];
        temp[3] = dof2;
        temp[4] = dof1;
        Base::operator[](index5DTo1D(numDOF, superSize, temp)) = T(-M_SQRT1_2);
    }

    template<Scalar T>
    void FCSwapVector<T>::swap(FCSwapVector& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(dof1, obj.dof1);
        std::swap(dof2, obj.dof2);
        std::swap(numDOF, obj.numDOF);
        cellIndex.swap(obj.cellIndex);
        superSize.swap(obj.superSize);
    }

    template<Scalar T>
    Index5D FCSwapVector<T>::index1DTo5D(size_t numDOF, const Index3D& superSize, size_t index1D) {
        assert(index1D < calcLength(numDOF, superSize) && "[Error]: Index out of range");
        Index5D result{};
        result[4] = index1D % numDOF;
        index1D /= numDOF;
        result[3] = index1D % numDOF;
        index1D /= numDOF;
        result[2] = index1D % superSize[2];
        index1D /= superSize[2];
        result[1] = index1D % superSize[1];
        result[0] = index1D / superSize[1];
        return result;
    }

    template<Scalar T>
    size_t FCSwapVector<T>::index5DTo1D(size_t numDOF, const Index3D& superSize, const Index5D& index5D) {
        checkParam(numDOF, superSize, index5D);
        return (((index5D[0] * superSize[1] + index5D[1]) * superSize[2] + index5D[2]) * numDOF + index5D[3]) * numDOF + index5D[4];
    }

    template<Scalar T>
    inline size_t FCSwapVector<T>::calcLength(size_t numDOF, const Index3D& superSize) {
        return numDOF * numDOF * superSize[0] * superSize[1] * superSize[2];
    }

    template<Scalar T>
    void FCSwapVector<T>::checkParam(
            [[maybe_unused]] size_t numDOF,
            [[maybe_unused]] const Index3D& superSize,
            [[maybe_unused]] const Index5D& index) {
        assert(index[0] < superSize[0] && "[Error]: Cell index out of range");
        assert(index[1] < superSize[1] && "[Error]: Cell index out of range");
        assert(index[2] < superSize[2] && "[Error]: Cell index out of range");
        assert(index[3] < numDOF && index[4] < numDOF && "[Error]: DOF out of range");
    }
}

namespace Physica {
    template<class T>
    class Traits<FCSwapVector<T>> {
    public:
        using ScalarType = T;
        constexpr static size_t SizeAtCompile = Dynamic;
    };
}
