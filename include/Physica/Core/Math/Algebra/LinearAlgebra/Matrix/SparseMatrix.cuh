/*
 * Copyright 2025-2026 Weibo He.
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

#include "SparseMatrix.h"
#include "MatrixImpl/RValueMatrix.cuh"
#include "Physica/Core/Utils/Container/Array.cuh"

namespace Physica {
    template<Scalar T, int Major>
    class device_obj<SparseMatrix<T, Major>> : public device_obj<RValueMatrix<SparseMatrix<T, Major>>> {
        using host_obj = SparseMatrix<T, Major>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    private:
        device_obj<Array<T>> elements;
        device_obj<Array<size_t>> minorIndexes;
        device_obj<Array<size_t>> majorStarts;
        size_t maxMinor = 0;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t col);
        device_obj(const host_obj& obj);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t, size_t) const { noImpl(); }

        using Base::resize;
        void resize(size_t row, size_t col);
        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;
        
        void zeros();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] size_t getMaxMajor() const noexcept;
        [[nodiscard]] size_t getMaxMinor() const noexcept;
        [[nodiscard]] size_t getNumNonZero() const noexcept { return elements.getLength(); }
        [[nodiscard]] auto&& getElements(this auto&&) noexcept;
        [[nodiscard]] const auto& getMinorIndexes() const noexcept { return minorIndexes; }
        [[nodiscard]] const auto& getMajorStarts() const noexcept { return majorStarts; }
        /* Friends */
        friend class SparseMatrix<T, Major>;
    };

    template<Scalar T, int Major>
    device_obj<SparseMatrix<T, Major>>::device_obj(size_t row, size_t col)
            : majorStarts(MatrixMajor::selectMajor<This>(row, col) + 1, 0)
            , maxMinor(MatrixMajor::selectMinor<This>(row, col)) {
        size_t size = std::max(row, col);
        elements.reserve(size);
        minorIndexes.reserve(size);
    }

    template<Scalar T, int Major>
    device_obj<SparseMatrix<T, Major>>::device_obj(const host_obj& obj) : device_obj(obj.getRow(), obj.getCol()) {
        obj.toDevice(*this);
    }

    template<Scalar T, int Major>
    void device_obj<SparseMatrix<T, Major>>::resize(size_t row, size_t col) {
        zeros();
        majorStarts.resize(MatrixMajor::selectMajor<This>(row, col) + 1, 0);
        maxMinor = MatrixMajor::selectMinor<This>(row, col);
    }

    template<Scalar T, int Major>
    auto device_obj<SparseMatrix<T, Major>>::toHost() const -> host_obj {
        auto result = toHostAsync();
        CUDAExecutor::wait();
        return result;
    }

    template<Scalar T, int Major>
    auto device_obj<SparseMatrix<T, Major>>::toHostAsync() const -> host_obj {
        host_obj result{};
        toHostAsync(result);
        return result;
    }

    template<Scalar T, int Major>
    void device_obj<SparseMatrix<T, Major>>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAExecutor::wait();
    }

    template<Scalar T, int Major>
    void device_obj<SparseMatrix<T, Major>>::toHostAsync(host_obj& obj) const {
        obj.resize(getRow(), getCol());
        elements.toHostAsync(obj.elements);
        majorStarts.toHostAsync(obj.majorStarts);
        minorIndexes.toHostAsync(obj.minorIndexes);
    }

    template<Scalar T, int Major>
    void device_obj<SparseMatrix<T, Major>>::zeros() {
        elements.resize(0);
        minorIndexes.resize(0);
        majorStarts.zeros();
    }

    template<Scalar T, int Major>
    void device_obj<SparseMatrix<T, Major>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        elements.swap(obj.elements);
        minorIndexes.swap(obj.minorIndexes);
        majorStarts.swap(obj.majorStarts);
        std::swap(maxMinor, obj.maxMinor);
    }

    template<Scalar T, int Major>
    size_t device_obj<SparseMatrix<T, Major>>::getRow() const noexcept {
        return MatrixMajor::isColMatrix<This>() ? getMaxMinor() : getMaxMajor();
    }

    template<Scalar T, int Major>
    size_t device_obj<SparseMatrix<T, Major>>::getCol() const noexcept {
        return MatrixMajor::isColMatrix<This>() ? getMaxMajor() : getMaxMinor();
    }

    template<Scalar T, int Major>
    size_t device_obj<SparseMatrix<T, Major>>::getMaxMajor() const noexcept {
        return majorStarts.getLength() - 1;
    }

    template<Scalar T, int Major>
    size_t device_obj<SparseMatrix<T, Major>>::getMaxMinor() const noexcept {
        return maxMinor;
    }

    template<Scalar T, int Major>
    auto&& device_obj<SparseMatrix<T, Major>>::getElements(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).elements;
    }

    template<Scalar T, int Major>
    auto SparseMatrix<T, Major>::toDevice() const -> device_obj<This> {
        auto result = toDeviceAsync();
        CUDAExecutor::wait();
        return result;
    }

    template<Scalar T, int Major>
    auto SparseMatrix<T, Major>::toDeviceAsync() const -> device_obj<This> {
        device_obj<This> result{};
        toDeviceAsync(result);
        return result;
    }

    template<Scalar T, int Major>
    void SparseMatrix<T, Major>::toDevice(device_obj<This>& obj) const {
        toDeviceAsync(obj);
        CUDAExecutor::wait();
    }

    template<Scalar T, int Major>
    void SparseMatrix<T, Major>::toDeviceAsync(device_obj<This>& obj) const {
        obj.resize(getRow(), getCol());
        elements.toDeviceAsync(obj.elements);
        majorStarts.toDeviceAsync(obj.majorStarts);
        minorIndexes.toDeviceAsync(obj.minorIndexes);
    }
}

namespace Physica {
    template<Scalar T, int Op>
    class Traits<device_obj<SparseMatrix<T, Op>>> : public Traits<SparseMatrix<T, Op>> {};
}

