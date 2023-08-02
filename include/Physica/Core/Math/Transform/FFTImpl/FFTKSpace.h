/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class Derived, size_t Dim> class FFTKSpace;

    namespace Internal {
        template<class Derived>
        class Traits<FFTKSpace<Derived, 1>> {
            static_assert(Traits<Derived>::Dim == 1, "[Error]: Inconsistent template param");
            using T = typename Traits<Derived>::ScalarType;
        public:
            using ScalarType = typename T::ComplexType;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = Dynamic;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class Derived>
        class Traits<FFTKSpace<Derived, 2>> {
            static_assert(Traits<Derived>::Dim == 2, "[Error]: Inconsistent template param");
            using T = typename Traits<Derived>::ScalarType;
        public:
            using ScalarType = typename T::ComplexType;
            constexpr static int Option = MatrixOption::Row | MatrixOption::Element;
            constexpr static size_t RowAtCompile = Dynamic;
            constexpr static size_t ColumnAtCompile = Dynamic;
            constexpr static size_t MaxRowAtCompile = RowAtCompile;
            constexpr static size_t MaxColumnAtCompile = ColumnAtCompile;
            constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
        };

        template<class Derived>
        class Traits<FFTKSpace<Derived, 3>> {
            static_assert(Traits<Derived>::Dim == 3, "[Error]: Inconsistent template param");
            using T = typename Traits<Derived>::ScalarType;
        public:
            using ScalarType = typename T::ComplexType;
        };
    }

    template<class Derived>
    class FFTKSpace<Derived, 1>
            : public Utils::CRTPBase<Derived, 1>
            , public ContinuousVector<FFTKSpace<Derived, 1>> {
        using This = FFTKSpace<Derived, 1>;
        using Base = Utils::CRTPBase<Derived, 1>;
        using VectorBase = ContinuousVector<This>;
        using RealType = typename Internal::Traits<This>::ScalarType::RealType;
    public:
        using typename VectorBase::ScalarType;
    public:
        /* Operators */
        using VectorBase::operator=;
        [[nodiscard]] inline ScalarType& operator[](size_t index);
        [[nodiscard]] inline const ScalarType& operator[](size_t index) const;
        /* Operations */
        template<class VectorType>
        inline void invTransform(const RValueVector<VectorType>& data);
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Derived::rSizeToKSize(Base::getDerived().getRSpaceSize()); }
        [[nodiscard]] size_t getSize() const noexcept { return getLength(); }
        [[nodiscard]] RealType getDelta() const noexcept { return reciprocal(Base::getDerived().getRSpaceDelta() * Base::getDerived().getRSpaceSize()); }
    };

    template<class Derived>
    inline typename FFTKSpace<Derived, 1>::ScalarType& FFTKSpace<Derived, 1>::operator[](size_t index) {
        return const_cast<ScalarType&>(const_cast<const This&>(*this).operator[](index));
    }

    template<class Derived>
    inline const typename FFTKSpace<Derived, 1>::ScalarType& FFTKSpace<Derived, 1>::operator[](size_t index) const {
        assert(index < getLength());
        return Base::getDerived().asComplexBuffer()[index];
    }

    template<class Derived>
    template<class VectorType>
    inline void FFTKSpace<Derived, 1>::invTransform(const RValueVector<VectorType>& data) {
        assert(data.getLength() == getLength());
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    class FFTKSpace<Derived, 2>
            : public Utils::CRTPBase<Derived, 1>
            , public ContinuousMatrix<FFTKSpace<Derived, 2>> {
        using This = FFTKSpace<Derived, 2>;
        using Base = Utils::CRTPBase<Derived, 1>;
        using MatrixBase = ContinuousMatrix<This>;
    public:
        using typename MatrixBase::ScalarType;
    public:
        /* Operators */
        using MatrixBase::operator=;
        [[nodiscard]] inline ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] inline const ScalarType& operator()(size_t row, size_t col) const;
        /* Operations */
        template<class MatrixType>
        inline void invTransform(const RValueMatrix<MatrixType>& data);
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow()); assert(col == getColumn()); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getDerived().getKSpaceSize()[0]; }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::getDerived().getKSpaceSize()[1]; }
    };

    template<class Derived>
    inline typename FFTKSpace<Derived, 2>::ScalarType& FFTKSpace<Derived, 2>::operator()(size_t row, size_t col) {
        return const_cast<ScalarType&>(const_cast<const This&>(*this).operator()(row, col));
    }

    template<class Derived>
    inline const typename FFTKSpace<Derived, 2>::ScalarType& FFTKSpace<Derived, 2>::operator()(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getColumn());
        return Base::getDerived().asComplexBuffer()[row * getColumn() + col];
    }

    template<class Derived>
    template<class MatrixType>
    inline void FFTKSpace<Derived, 2>::invTransform(const RValueMatrix<MatrixType>& data) {
        assert(data.getRow() == getRow());
        assert(data.getColumn() == getColumn());
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    class FFTKSpace<Derived, 3>
            : public Utils::CRTPBase<Derived, 1>
            , public LValueGrid<FFTKSpace<Derived, 3>> {
        using This = FFTKSpace<Derived, 3>;
        using Base = Utils::CRTPBase<Derived, 1>;
        using GridBase = LValueGrid<This>;
    public:
        using typename GridBase::ScalarType;
    public:
        /* Operators */
        [[nodiscard]] inline ScalarType& operator()(size_t x, size_t y, size_t z);
        [[nodiscard]] inline const ScalarType& operator()(size_t x, size_t y, size_t z) const;
        /* Operations */
        template<class GridType>
        inline void invTransform(const LValueGrid<GridType>& data);
        /* Getters */
        [[nodiscard]] size_t getDimX() const noexcept { return Base::getDerived().getKSpaceSize()[0]; }
        [[nodiscard]] size_t getDimY() const noexcept { return Base::getDerived().getKSpaceSize()[1]; }
        [[nodiscard]] size_t getDimZ() const noexcept { return Base::getDerived().getKSpaceSize()[2]; }
    };

    template<class Derived>
    inline typename FFTKSpace<Derived, 3>::ScalarType& FFTKSpace<Derived, 3>::operator()(size_t x, size_t y, size_t z) {
        return const_cast<ScalarType&>(const_cast<const This&>(*this).operator()(x, y, z));
    }

    template<class Derived>
    inline const typename FFTKSpace<Derived, 3>::ScalarType& FFTKSpace<Derived, 3>::operator()(size_t x, size_t y, size_t z) const {
        assert(x < getDimX());
        assert(y < getDimY());
        assert(z < getDimZ());
        return Base::getDerived().asComplexBuffer()[(x * getDimY() + y) * getDimZ() + z];
    }

    template<class Derived>
    template<class GridType>
    inline void FFTKSpace<Derived, 3>::invTransform(const LValueGrid<GridType>& data) {
        assert(data.getDimX() == getDimX());
        assert(data.getDimY() == getDimY());
        assert(data.getDimZ() == getDimZ());
        *this = data;
        Base::getDerived().invTransform();
    }
}
