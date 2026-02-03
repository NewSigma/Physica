/*
 * Copyright 2021-2026 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<class T>
    struct remove_transpose {
        using Type = T;
    };

    template<class T>
    struct remove_transpose<Transpose<T>> {
        using Type = T;
    };

    template<class T>
    struct is_transpose {
        constexpr static bool value = !std::is_same<T, typename remove_transpose<T>::Type>::value;
    };

    template<Matrix M>
    class Transpose<M> : public RValueMatrix<Transpose<M>> {
        using This = Transpose<M>;
        using Base = RValueMatrix<This>;

        LazyDestroy<M> mat;
    public:        
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
    public:
        Transpose(M&& mat_) : mat(std::forward<M>(mat_)) {}
        Transpose(const This&) = default;
        Transpose(This&&) noexcept = default;
        ~Transpose() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& target) const;
        void assign_mkl(Matrix auto&& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const { return mat.calc(col, row); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return mat.calc_value(col, row); }

        [[nodiscard]] const M& transpose() const noexcept { return mat; }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getRow() const noexcept { return mat.getCol(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getRow(); }
    };

    template<Matrix M>
    template<ExecutePolicy P>
    void Transpose<M>::assign(Matrix auto&& target) const {
        constexpr bool LargeMatrix = Traits<This>::SizeAtCompile == Dynamic;
        if constexpr (LargeMatrix && Internal::EnableMKL<M, decltype(target)>::value && MatrixMajor::isSameMajor<M, decltype(target)>()) {
            if (Base::getSize() <= 16)
                Base::template assign<P>(target);
            else
                assign_mkl(target);
        }
        else
            Base::template assign<P>(target);
    }

    template<Matrix M>
    auto&& Transpose<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Vector V>
    class Transpose<V> : public RValueMatrix<Transpose<V>> {
        using This = Transpose<V>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<V> vec;
    public:
        explicit Transpose(V&& vec_) : vec(std::forward<V>(vec_)) {}
        Transpose(const This&) = default;
        Transpose(This&&) noexcept = default;
        ~Transpose() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col); }
        [[nodiscard]] Tv calc_value([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc_value(col); }

        [[nodiscard]] const V& transpose() const noexcept { return vec; }
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getCol() const noexcept { return vec.getLength(); }
    };

    template<Vector V>
    void Transpose<V>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(target.rowFromMajorMinor(0, i), target.colFromMajorMinor(0, i));
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Transpose<M>> {
        using M1 = std::remove_cvref_t<M>;
        constexpr static int OtherMajor = MatrixMajor::isColMatrix<M>() ? MatrixMajor::Row : MatrixMajor::Col;
        constexpr static int Major = MatrixMajor::isBothMajor<M>() ? MatrixMajor::BothMajor : OtherMajor;
    public:
        using ScalarType = M1::ScalarType;
        constexpr static int Option = Major;
        constexpr static size_t RowAtCompile = M1::ColAtCompile;
        constexpr static size_t ColAtCompile = M1::RowAtCompile;
        constexpr static size_t SizeAtCompile = M1::SizeAtCompile;
    };

    template<Vector V>
    class Traits<Transpose<V>> {
        using V1 = std::remove_cvref_t<V>;
    public:
        using ScalarType = V1::ScalarType;
        constexpr static int Option = MatrixMajor::Row;
        constexpr static size_t RowAtCompile = 1;
        constexpr static size_t ColAtCompile = V1::SizeAtCompile;
        constexpr static size_t SizeAtCompile = V1::SizeAtCompile;
    };
}

#ifdef PHYSICA_MKL
    #include "Transpose_MKL.h"
#endif
