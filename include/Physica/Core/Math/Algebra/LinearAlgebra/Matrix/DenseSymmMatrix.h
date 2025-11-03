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
#include "Physica/Core/Utils/Container/SymmArray.h"

namespace Physica {
    template<Scalar T, size_t Order, Vector V>
    class SyMV;

    template<Scalar T, size_t Order = Dynamic>
    class DenseSymmMatrix : public LValueMatrix<DenseSymmMatrix<T, Order>>
                          , private SymmArray<T, Order> {
        using This = DenseSymmMatrix<T, Order>;
        using Base = LValueMatrix<This>;
        using Storage = SymmArray<T, Order>;
        using VectorStorage = Storage::ArrayType;
    public:
        using typename Base::ScalarType;
        using ColMatrix = This;
        using RowMatrix = This;
    private:
        using HermiteRtnTy = std::conditional<ScalarType::isComplex, Hermite<This>, const This&>::type;
    public:
        DenseSymmMatrix() = default;
        DenseSymmMatrix(const Matrix auto& mat);
        using Storage::Storage;
        DenseSymmMatrix(const This&) = default;
        DenseSymmMatrix(This&&) noexcept = default;
        ~DenseSymmMatrix() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator();
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        This& operator=(const Scalar auto& x);
        void operator+=(const Scalar auto& x);
        void operator-=(const Scalar auto& x);
        void operator*=(const Scalar auto& x);
        void operator/=(const Scalar auto& x);

        template<Vector V>
        [[nodiscard]] auto operator*(const V& v) const noexcept;
        /* Operations */
        using Base::assign;
        using Base::calc;
        using Storage::resize;
        void resize(const Matrix auto& m, auto&&... args);

        [[nodiscard]] ScalarType max() const { return asVector().max(); }
        [[nodiscard]] ScalarType min() const { return asVector().min(); }
        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return *this; }
        void swap(This& __restrict m) noexcept;

        using Base::format;

        using Storage::zeros;
        template<RNG R = Random<>>
        void random_uniform() { asVector().template random_uniform<R>(); }
        template<RNG R = Random<>>
        void random_normal() { asVector().template random_normal<R>(); }
        template<RNG R = Random<>>
        void random_any(auto& distribution) { asVector().template random_any<R>(distribution); }

        using Storage::read;
        using Storage::write;
        /* Getters */
        using Storage::data_ptr;
        using Storage::getCol;
        using Storage::getRow;
        using Storage::getOrder;
        using Storage::toIndex1D;
        [[nodiscard]] VectorStorage& asVector() noexcept { return Storage::asArray(); }
        [[nodiscard]] const VectorStorage& asVector() const noexcept { return Storage::asArray(); }
        /* Static members */
        [[nodiscard]] static DenseSymmMatrix unitMatrix(size_t order);
        template<RNG R = Random<>>
        [[nodiscard]] static This random_uniform(size_t order);
        template<RNG R = Random<>>
        [[nodiscard]] static This random_normal(size_t order);
        template<RNG R = Random<>>
        [[nodiscard]] static This random_any(size_t order, auto& distribution);
    };
    /**
     * Assuming mat is a symmetric matrix, if it is not the case, only half of the elements is saved correctly
     */
    template<Scalar T, size_t Order>
    DenseSymmMatrix<T, Order>::DenseSymmMatrix(const Matrix auto& mat) : DenseSymmMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getCol());
        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = i; j < mat.getRow(); ++j)
                Base::operator()(i, j) = mat.calc(i, j);
    }

    template<Scalar T, size_t Order>
    DenseSymmMatrix<T, Order>& DenseSymmMatrix<T, Order>::operator=(const Scalar auto& x) {
        asVector() = x;
        return *this;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator+=(const Scalar auto& x) {
        asVector() += x;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator-=(const Scalar auto& x) {
        asVector() -= x;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator*=(const Scalar auto& x) {
        asVector() *= x;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator/=(const Scalar auto& x) {
        asVector() /= x;
    }

    template<Scalar T, size_t Order>
    template<Vector V>
    auto DenseSymmMatrix<T, Order>::operator*(const V& v) const noexcept {
        return SyMV<T, Order, V>(*this, v);
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::resize(const Matrix auto& m, auto&&... args) {
        Base::resize(m, std::forward<decltype(args)>(args)...);
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
    template<RNG R>
    auto DenseSymmMatrix<T, Order>::random_uniform(size_t order) -> This {
        This result(order);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseSymmMatrix<T, Order>::random_normal(size_t order) -> This {
        This result(order);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseSymmMatrix<T, Order>::random_any(size_t order, auto& distribution) -> This {
        This result(order);
        result.template random_any<R>(distribution);
        return result;
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<DenseSymmMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

namespace std {
    template<Physica::Scalar T, size_t Order>
    void swap(Physica::DenseSymmMatrix<T, Order>& __restrict m1,
                     Physica::DenseSymmMatrix<T, Order>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

#include "MatrixImpl/MatrixProduct/SyMV.h"
