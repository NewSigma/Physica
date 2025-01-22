/*
 * Copyright 2020-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"
#include "DenseSymmImpl/HalfDenseMatrixStorage.h"

namespace Physica::Core {
    template<Scalar T, size_t Order = Dynamic>
    class DenseSymmMatrix : public LValueMatrix<DenseSymmMatrix<T, Order>>
                          , private HalfDenseMatrixStorage<T, Order> {
        using This = DenseSymmMatrix<T, Order>;
        using Base = LValueMatrix<This>;
        using Storage = HalfDenseMatrixStorage<T, Order>;
        using VectorStorage = Storage::ArrayType;
    public:
        using typename Base::ScalarType;
        using ColMatrix = This;
        using RowMatrix = This;
    public:
        template<Matrix M>
        DenseSymmMatrix(const M& mat);
        using Storage::Storage;
        DenseSymmMatrix(const This&) = default;
        DenseSymmMatrix(This&&) noexcept = default;
        ~DenseSymmMatrix() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator();
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        template<Scalar U> inline This& operator=(const U& x);
        template<Scalar U> inline void operator+=(const U& x);
        template<Scalar U> inline void operator-=(const U& x);
        template<Scalar U> inline void operator*=(const U& x);
        template<Scalar U> inline void operator/=(const U& x);
        /* Operations */
        using Base::assign;
        using Base::calc;
        using Base::format;
        using Storage::resize;
        [[nodiscard]] ScalarType max() const { return asVector().max(); }
        [[nodiscard]] ScalarType min() const { return asVector().min(); }
        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        [[nodiscard]] inline auto hermite() const noexcept;
        void swap(This& __restrict m) noexcept;

        template<RandomGenerator R>
        void random_uniform() { asVector().template random_uniform<R>(); }
        template<RandomGenerator R>
        void random_normal() { asVector().template random_normal<R>(); }
        template<class Distribution, RandomGenerator R>
        void random_any(Distribution& dist) { asVector().template random_any<Distribution, R>(dist); }

        inline const H5DataSet<1> read(const H5Location& loc, const char* name);
        inline H5DataSet<1> write(H5Location& loc, const char* name) const;
        /* Getters */
        using Storage::data_ptr;
        using Storage::getCol;
        using Storage::getOrder;
        using Storage::getRow;
        using Storage::toIndex1D;
        [[nodiscard]] VectorStorage& asVector() noexcept { return Storage::asArray(); }
        [[nodiscard]] const VectorStorage& asVector() const noexcept { return Storage::asArray(); }
        /* Static members */
        [[nodiscard]] static DenseSymmMatrix unitMatrix(size_t order);
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_uniform(size_t order);
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_normal(size_t order);
        template<class Distribution, RandomGenerator R>
        [[nodiscard]] inline static This random_any(size_t order, Distribution& dist);
    };
    /**
     * Assuming mat is a symmetric matrix, if it is not the case, only half of the elements is saved correctly
     */
    template<Scalar T, size_t Order>
    template<Matrix M>
    DenseSymmMatrix<T, Order>::DenseSymmMatrix(const M& mat)
            : DenseSymmMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getCol());
        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = i; j < mat.getRow(); ++j)
                Base::operator()(i, j) = T(mat.calc(i, j));
    }

    template<Scalar T, size_t Order>
    template<Scalar U>
    inline DenseSymmMatrix<T, Order>& DenseSymmMatrix<T, Order>::operator=(const U& x) {
        asVector() = x;
        return *this;
    }

    template<Scalar T, size_t Order>
    template<Scalar U>
    inline void DenseSymmMatrix<T, Order>::operator+=(const U& x) {
        asVector() += x;
    }

    template<Scalar T, size_t Order>
    template<Scalar U>
    inline void DenseSymmMatrix<T, Order>::operator-=(const U& x) {
        asVector() -= x;
    }

    template<Scalar T, size_t Order>
    template<Scalar U>
    inline void DenseSymmMatrix<T, Order>::operator*=(const U& x) {
        asVector() *= x;
    }

    template<Scalar T, size_t Order>
    template<Scalar U>
    inline void DenseSymmMatrix<T, Order>::operator/=(const U& x) {
        asVector() /= x;
    }

    template<Scalar T, size_t Order>
    inline auto DenseSymmMatrix<T, Order>::hermite() const noexcept {
        if constexpr (T::isComplex)
            return Base::hermite();
        else
            return *this;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::swap(This& __restrict m) noexcept {
        assert(this != &m && "[Error]: Self swap is likely a bug");
        Storage::swap(m);
    }

    template<Scalar T, size_t Order>
    DenseSymmMatrix<T, Order> DenseSymmMatrix<T, Order>::unitMatrix(size_t order) {
        DenseSymmMatrix<T, Order> result(order);
        result.toUnitMatrix();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RandomGenerator R>
    [[nodiscard]] inline DenseSymmMatrix<T, Order>
    DenseSymmMatrix<T, Order>::random_uniform(size_t order) {
        This result(order);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RandomGenerator R>
    [[nodiscard]] inline DenseSymmMatrix<T, Order>
    DenseSymmMatrix<T, Order>::random_normal(size_t order) {
        This result(order);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<class Distribution, RandomGenerator R>
    [[nodiscard]] inline DenseSymmMatrix<T, Order>
    DenseSymmMatrix<T, Order>::random_any(size_t order, Distribution& dist) {
        This result(order);
        result.template random_any<Distribution, R>(dist);
        return result;
    }
#ifdef PHYSICA_HDF5
    template<Scalar T, size_t Order>
    inline const H5DataSet<1> DenseSymmMatrix<T, Order>::read(const H5Location& loc, const char* name) {
        return asVector().read(loc, name);
    }

    template<Scalar T, size_t Order>
    inline H5DataSet<1> DenseSymmMatrix<T, Order>::write(H5Location& loc, const char* name) const {
        return asVector().write(loc, name);
    }
#endif
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<DenseSymmMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

namespace std {
    template<Physica::Core::Scalar T, size_t Order>
    inline void swap(Physica::Core::DenseSymmMatrix<T, Order>& __restrict m1,
                     Physica::Core::DenseSymmMatrix<T, Order>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

#include "DenseSymmImpl/MatrixVectorProduct.h"
