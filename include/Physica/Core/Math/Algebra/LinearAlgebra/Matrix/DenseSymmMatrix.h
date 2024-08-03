/*
 * Copyright 2020-2024 Weibo He.
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

#include "DenseMatrix.h"
#include "DenseMatrixImpl/HalfDenseMatrixStorage.h"

namespace Physica::Core {
    template<class T, size_t Order = Dynamic, size_t MaxOrder = Order>
    class DenseSymmMatrix : public LValueMatrix<DenseSymmMatrix<T, Order, MaxOrder>>
                          , private HalfDenseMatrixStorage<T, Order, MaxOrder> {
        using This = DenseSymmMatrix<T, Order, MaxOrder>;
        using Base = LValueMatrix<This>;
        using Storage = HalfDenseMatrixStorage<T, Order, MaxOrder>;
        using VectorStorage = typename Storage::ArrayType;
    public:
        using typename Base::ScalarType;
        using ColMatrix = This;
        using RowMatrix = This;
        using RealMatrix = DenseSymmMatrix<typename ScalarType::RealType, Order, MaxOrder>;
    public:
        template<class OtherMatrix>
        DenseSymmMatrix(const RValueMatrix<OtherMatrix>& mat);
        using Storage::Storage;
        DenseSymmMatrix(const This&) = default;
        DenseSymmMatrix(This&&) noexcept = default;
        ~DenseSymmMatrix() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator();
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        inline This& operator=(const ScalarType& value);
        template<class VectorType>
        [[nodiscard]] inline MatrixVectorProduct<This, VectorType> operator*(const RValueVector<VectorType>& vec) const noexcept;
        /* Operations */
        using Base::assignTo;
        using Base::calc;
        using Base::format;
        using Base::hermite;
        using Storage::resize;
        [[nodiscard]] ScalarType max() const { return asVector().max(); }
        [[nodiscard]] ScalarType min() const { return asVector().min(); }
        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        void swap(This& __restrict m) noexcept;

        template<class RandomGenerator>
        void random_uniform(RandomGenerator& gen) { asVector().random_uniform(gen); }
        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen) { asVector().random_normal(gen); }
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen) { asVector().random_any(dist, gen); }

        inline const H5DataSet<1> read(const H5Location& loc, const char* name, const H5::DSetMemXferPropList& xfer_plist = H5::DSetMemXferPropList::DEFAULT);
        inline H5DataSet<1> write(H5Location& loc, const char* name, const H5::DSetMemXferPropList& xfer_plist = H5::DSetMemXferPropList::DEFAULT) const;
        /* Getters */
        using Storage::data_ptr;
        using Storage::getColumn;
        using Storage::getOrder;
        using Storage::getRow;
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorStorage& asVector() noexcept { return Storage::getArray(); }
        [[nodiscard]] const VectorStorage& asVector() const noexcept { return Storage::getArray(); }
        /* Static members */
        [[nodiscard]] static DenseSymmMatrix unitMatrix(size_t order);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t order, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t order, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t order, Distribution& dist, RandomGenerator& gen);
    };
    /**
     * Assuming mat is a symmetric matrix, if it is not the case, only half of the elements is saved correctly
     */
    template<class T, size_t Order, size_t MaxOrder>
    template<class OtherMatrix>
    DenseSymmMatrix<T, Order, MaxOrder>::DenseSymmMatrix(const RValueMatrix<OtherMatrix>& mat)
            : DenseSymmMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getColumn());
        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = i; j < mat.getRow(); ++j)
                Base::operator()(i, j) = mat.calc(i, j);
    }

    template<class T, size_t Order, size_t MaxOrder>
    inline DenseSymmMatrix<T, Order, MaxOrder>& DenseSymmMatrix<T, Order, MaxOrder>::operator=(const ScalarType& value) {
        asVector() = value;
        return *this;
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class VectorType>
    inline MatrixVectorProduct<DenseSymmMatrix<T, Order, MaxOrder>, VectorType>
    DenseSymmMatrix<T, Order, MaxOrder>::operator*(const RValueVector<VectorType>& vec) const noexcept {
        return static_cast<const Base&>(*this) * vec;
    }

    template<class T, size_t Order, size_t MaxOrder>
    void DenseSymmMatrix<T, Order, MaxOrder>::swap(This& __restrict m) noexcept {
        assert(this != &m && "[Error]: Self swap is likely a bug");
        Storage::swap(m);
    }

    template<class T, size_t Order, size_t MaxOrder>
    DenseSymmMatrix<T, Order, MaxOrder> DenseSymmMatrix<T, Order, MaxOrder>::unitMatrix(size_t order) {
        DenseSymmMatrix<T, Order, MaxOrder> result(order);
        result.toUnitMatrix();
        return result;
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class RandomGenerator>
    [[nodiscard]] inline DenseSymmMatrix<T, Order, MaxOrder>
    DenseSymmMatrix<T, Order, MaxOrder>::random_uniform(size_t order, RandomGenerator& gen) {
        This result(order);
        result.random_uniform(gen);
        return result;
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class RandomGenerator>
    [[nodiscard]] inline DenseSymmMatrix<T, Order, MaxOrder>
    DenseSymmMatrix<T, Order, MaxOrder>::random_normal(size_t order, RandomGenerator& gen) {
        This result(order);
        result.random_normal(gen);
        return result;
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class Distribution, class RandomGenerator>
    [[nodiscard]] inline DenseSymmMatrix<T, Order, MaxOrder>
    DenseSymmMatrix<T, Order, MaxOrder>::random_any(size_t order, Distribution& dist, RandomGenerator& gen) {
        This result(order);
        result.random_any(dist, gen);
        return result;
    }

    template<class T, size_t Order, size_t MaxOrder>
    inline const H5DataSet<1> DenseSymmMatrix<T, Order, MaxOrder>::read(
            const H5Location& loc, const char* name, const H5::DSetMemXferPropList& xfer_plist) {
        return asVector().read(loc, name, xfer_plist);
    }

    template<class T, size_t Order, size_t MaxOrder>
    inline H5DataSet<1> DenseSymmMatrix<T, Order, MaxOrder>::write(
            H5Location& loc, const char* name, const H5::DSetMemXferPropList& xfer_plist) const {
        return asVector().write(loc, name, xfer_plist);
    }
}

namespace Physica {
    using namespace Core;

    template<class T, size_t Order, size_t MaxOrder>
    class Traits<DenseSymmMatrix<T, Order, MaxOrder>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColumnAtCompile = Order;
        constexpr static size_t MaxRowAtCompile = MaxOrder;
        constexpr static size_t MaxColumnAtCompile = MaxOrder;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
        constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
    };
}

namespace std {
    template<class T, size_t Order, size_t MaxOrder>
    inline void swap(Physica::Core::DenseSymmMatrix<T, Order, MaxOrder>& __restrict m1,
                     Physica::Core::DenseSymmMatrix<T, Order, MaxOrder>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}
