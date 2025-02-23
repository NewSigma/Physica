/*
 * Copyright 2021-2025 Weibo He.
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

#include "TensorImpl/TensorStorage.h"
#include "TensorImpl/LValueTensor.h"

namespace Physica {
    template<Scalar T>
    class DenseTensor : public LValueTensor<DenseTensor<T>>, private TensorStorage<T> {
        using This = DenseTensor<T>;
        using Base = LValueTensor<This>;
        using Storage = TensorStorage<T>;
    public:
        DenseTensor() = default;
        template<class... Args>
        DenseTensor(Index3D index, Args&&... args);
        DenseTensor(const This&) = default;
        DenseTensor(This&&) noexcept = default;
        ~DenseTensor() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(This grid) noexcept;
        using Storage::operator();
        template<class U> friend std::ostream& operator<<(std::ostream& os, const DenseTensor<U>& grid);
        template<class U> friend std::istream& operator>>(std::istream& is, DenseTensor<U>& grid);
        /* Operations */
        using Base::random_normal;
        template<class... Args>
        inline void resize(Index3D index, Args&&... args);
        inline void swap(This& __restrict grid) noexcept;
        /* Getters */
        using Storage::data_ptr;
        using Storage::getDim;
        using Storage::getDimX;
        using Storage::getDimY;
        using Storage::getDimZ;
        using Storage::getSize;
        /* Static members */
        using Base::forIndexInTensor;
        using Base::forPointIndexInTensor;
        using Base::forPointInTensor;
        template<RNG R>
        static DenseTensor random_uniform(Index3D size);
        template<RNG R>
        static DenseTensor random_normal(Index3D size);
    };

    template<Scalar T>
    template<class... Args>
    DenseTensor<T>::DenseTensor(Index3D index, Args&&... args) : Storage(index, std::forward<Args>(args)...) {}

    template<Scalar T>
    DenseTensor<T>& DenseTensor<T>::operator=(DenseTensor grid) noexcept {
        swap(grid);
        return *this;
    }

    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const DenseTensor<T>& grid) {
        const Index3D dim = grid.getDim();
        os.write(reinterpret_cast<const char*>(&dim), sizeof(Index3D));
        os.write(reinterpret_cast<const char*>(grid.asArray().data()), grid.getSize() * sizeof(T));
        return os;
    }

    template<Scalar T>
    std::istream& operator>>(std::istream& is, DenseTensor<T>& grid) {
        Index3D dim;
        is.read(reinterpret_cast<char*>(&dim), sizeof(Index3D));
        grid.resize(dim);
        is.read(reinterpret_cast<char*>(grid.asArray().data()), grid.getSize() * sizeof(T));
        return is;
    }

    template<Scalar T>
    template<class... Args>
    inline void DenseTensor<T>::resize(Index3D index, Args&&... args) {
        Storage::resize(index, std::forward<Args>(args)...);
    }

    template<Scalar T>
    inline void DenseTensor<T>::swap(This& __restrict grid) noexcept {
        assert(this != &grid && "[Error]: Self swap is likely a bug");
        Storage::swap(grid);
    }

    template<Scalar T>
    template<RNG R>
    inline DenseTensor<T> DenseTensor<T>::random_uniform(Index3D size) {
        auto result = DenseTensor<T>(size);
        result.flatten().template random_uniform<R>();
        return result;
    }

    template<Scalar T>
    template<RNG R>
    inline DenseTensor<T> DenseTensor<T>::random_normal(Index3D size) {
        auto result = DenseTensor<T>(size);
        result.flatten().template random_normal<R>();
        return result;
    }

    template<Scalar T>
    inline void swap(DenseTensor<T>& __restrict grid1, DenseTensor<T>& __restrict grid2) noexcept {
        grid1.swap(grid2);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<DenseTensor<T>> {
    public:
        using ScalarType = T;
        constexpr static int Dim = 3;
    };
}
