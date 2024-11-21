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
        using RealMatrix = DenseSymmMatrix<typename ScalarType::RealType, Order>;
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
        inline This& operator=(const ScalarType& value);
        template<Vector V>
        [[nodiscard]] inline auto operator*(const V& vec) const noexcept;
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

        template<class RandomType>
        void random_uniform(RandomType& gen) { asVector().random_uniform(gen); }
        template<class RandomType>
        void random_normal(RandomType& gen) { asVector().random_normal(gen); }
        template<class Distribution, class RandomType>
        void random_any(Distribution& dist, RandomType& gen) { asVector().random_any(dist, gen); }

        inline const H5DataSet<1> read(const H5Location& loc, const char* name);
        inline H5DataSet<1> write(H5Location& loc, const char* name) const;
        /* Getters */
        using Storage::data_ptr;
        using Storage::getCol;
        using Storage::getOrder;
        using Storage::getRow;
        [[nodiscard]] VectorStorage& asVector() noexcept { return Storage::asArray(); }
        [[nodiscard]] const VectorStorage& asVector() const noexcept { return Storage::asArray(); }
        /* Static members */
        [[nodiscard]] static DenseSymmMatrix unitMatrix(size_t order);
        template<class RandomType>
        [[nodiscard]] inline static This random_uniform(size_t order, RandomType& gen);
        template<class RandomType>
        [[nodiscard]] inline static This random_normal(size_t order, RandomType& gen);
        template<class Distribution, class RandomType>
        [[nodiscard]] inline static This random_any(size_t order, Distribution& dist, RandomType& gen);
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
                Base::operator()(i, j) = mat.calc(i, j);
    }

    template<Scalar T, size_t Order>
    inline DenseSymmMatrix<T, Order>& DenseSymmMatrix<T, Order>::operator=(const ScalarType& value) {
        asVector() = value;
        return *this;
    }

    template<Scalar T, size_t Order>
    template<Vector V>
    inline auto DenseSymmMatrix<T, Order>::operator*(const V& vec) const noexcept {
        return MatrixVectorProduct<This, V>(*this, vec);
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
    template<class RandomType>
    [[nodiscard]] inline DenseSymmMatrix<T, Order>
    DenseSymmMatrix<T, Order>::random_uniform(size_t order, RandomType& gen) {
        This result(order);
        result.random_uniform(gen);
        return result;
    }

    template<Scalar T, size_t Order>
    template<class RandomType>
    [[nodiscard]] inline DenseSymmMatrix<T, Order>
    DenseSymmMatrix<T, Order>::random_normal(size_t order, RandomType& gen) {
        This result(order);
        result.random_normal(gen);
        return result;
    }

    template<Scalar T, size_t Order>
    template<class Distribution, class RandomType>
    [[nodiscard]] inline DenseSymmMatrix<T, Order>
    DenseSymmMatrix<T, Order>::random_any(size_t order, Distribution& dist, RandomType& gen) {
        This result(order);
        result.random_any(dist, gen);
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
