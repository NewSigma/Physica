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
    template<class Derived, size_t Dim> class FFTRSpace;

    namespace Internal {
        template<class Derived>
        class Traits<FFTRSpace<Derived, 1>> {
            static_assert(Traits<Derived>::Dim == 1, "[Error]: Inconsistent template param");
        public:
            using ScalarType = typename Traits<Derived>::ScalarType;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = Dynamic;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class Derived>
        class Traits<FFTRSpace<Derived, 2>> {
            static_assert(Traits<Derived>::Dim == 2, "[Error]: Inconsistent template param");
        public:
            using ScalarType = typename Traits<Derived>::ScalarType;
            constexpr static int Option = MatrixOption::Row | MatrixOption::Vector;
            constexpr static size_t RowAtCompile = Dynamic;
            constexpr static size_t ColumnAtCompile = Dynamic;
            constexpr static size_t MaxRowAtCompile = RowAtCompile;
            constexpr static size_t MaxColumnAtCompile = ColumnAtCompile;
            constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
        };

        template<class Derived>
        class Traits<FFTRSpace<Derived, 3>> {
            static_assert(Traits<Derived>::Dim == 3, "[Error]: Inconsistent template param");
            using T = typename Traits<Derived>::ScalarType;
        public:
            using ScalarType = typename Traits<Derived>::ScalarType;
        };
    }

    template<class Derived>
    class FFTRSpace<Derived, 1>
            : public Utils::CRTPBase<Derived, 0>
            , public ContinuousVector<FFTRSpace<Derived, 1>> {
        using This = FFTRSpace<Derived, 1>;
        using Base = Utils::CRTPBase<Derived, 0>;
        using VectorBase = ContinuousVector<This>;
    public:
        using typename VectorBase::ScalarType;
    private:
        static constexpr bool isComplex = ScalarType::isComplex;
    public:
        /* Operators */
        using VectorBase::operator=;
        [[nodiscard]] inline ScalarType& operator[](size_t index);
        [[nodiscard]] inline const ScalarType& operator[](size_t index) const;
        /* Operations */
        template<class VectorType>
        inline void transform(const RValueVector<VectorType>& data);
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().rSpaceSize; }
        [[nodiscard]] size_t getSize() const noexcept { return getLength(); }
    };

    template<class Derived>
    inline typename FFTRSpace<Derived, 1>::ScalarType& FFTRSpace<Derived, 1>::operator[](size_t index) {
        return const_cast<ScalarType&>(const_cast<const This&>(*this).operator[](index));
    }

    template<class Derived>
    inline const typename FFTRSpace<Derived, 1>::ScalarType& FFTRSpace<Derived, 1>::operator[](size_t index) const {
        if constexpr (isComplex)
            return Base::getDerived().complex_buffer[index];
        else
            return Base::getDerived().real_buffer[index];
    }

    template<class Derived>
    template<class VectorType>
    inline void FFTRSpace<Derived, 1>::transform(const RValueVector<VectorType>& data) {
        assert(data.getLength() == getLength());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    class FFTRSpace<Derived, 2>
            : public Utils::CRTPBase<Derived, 0>
            , public ContinuousMatrix<FFTRSpace<Derived, 2>> {
        using This = FFTRSpace<Derived, 2>;
        using Base = Utils::CRTPBase<Derived, 0>;
        using MatrixBase = ContinuousMatrix<This>;
    public:
        using typename MatrixBase::ScalarType;
    private:
        static constexpr bool isComplex = ScalarType::isComplex;
    public:
        /* Operators */
        using MatrixBase::operator=;
        [[nodiscard]] inline ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] inline const ScalarType& operator()(size_t row, size_t col) const;
        /* Operations */
        template<class MatrixType>
        inline void transform(const RValueMatrix<MatrixType>& data);
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow()); assert(col == getColumn()); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getDerived().getRSpaceSize()[0]; }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::getDerived().getRSpaceSize()[1]; }
    };

    template<class Derived>
    inline typename FFTRSpace<Derived, 2>::ScalarType& FFTRSpace<Derived, 2>::operator()(size_t row, size_t col) {
        return const_cast<ScalarType&>(const_cast<const This&>(*this).operator()(row, col));
    }

    template<class Derived>
    inline const typename FFTRSpace<Derived, 2>::ScalarType& FFTRSpace<Derived, 2>::operator()(size_t row, size_t col) const {
        const size_t numColumn = getColumn();
        if constexpr (isComplex)
            return Base::getDerived().complex_buffer[row * numColumn + col];
        else {
            const size_t shift = (numColumn / 2 + 1) * 2;
            return Base::getDerived().real_buffer[row * shift + col];
        }
    }

    template<class Derived>
    template<class MatrixType>
    inline void FFTRSpace<Derived, 2>::transform(const RValueMatrix<MatrixType>& data) {
        assert(data.getRow() == getRow());
        assert(data.getColumn() == getColumn());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    class FFTRSpace<Derived, 3>
            : public Utils::CRTPBase<Derived, 0>
            , public LValueGrid<FFTRSpace<Derived, 3>> {
        using This = FFTRSpace<Derived, 3>;
        using Base = Utils::CRTPBase<Derived, 0>;
        using GridBase = LValueGrid<This>;
    public:
        using typename GridBase::ScalarType;
    private:
        static constexpr bool isComplex = ScalarType::isComplex;
    public:
        /* Operators */
        [[nodiscard]] inline ScalarType& operator()(size_t x, size_t y, size_t z);
        [[nodiscard]] inline const ScalarType& operator()(size_t x, size_t y, size_t z) const;
        /* Operations */
        template<class GridType>
        inline void transform(const LValueGrid<GridType>& data);
        /* Getters */
        [[nodiscard]] size_t getDimX() const noexcept { return Base::getDerived().getRSpaceSize()[0]; }
        [[nodiscard]] size_t getDimY() const noexcept { return Base::getDerived().getRSpaceSize()[1]; }
        [[nodiscard]] size_t getDimZ() const noexcept { return Base::getDerived().getRSpaceSize()[2]; }
    };

    template<class Derived>
    inline typename FFTRSpace<Derived, 3>::ScalarType& FFTRSpace<Derived, 3>::operator()(size_t x, size_t y, size_t z) {
        return const_cast<ScalarType&>(const_cast<const This&>(*this).operator()(x, y, z));
    }

    template<class Derived>
    inline const typename FFTRSpace<Derived, 3>::ScalarType& FFTRSpace<Derived, 3>::operator()(size_t x, size_t y, size_t z) const {
        if constexpr (isComplex)
            return Base::getDerived().complex_buffer[(x * getDimY() + y) * getDimZ() + z];
        else {
            const size_t shift = (getDimZ() / 2 + 1) * 2;
            return Base::getDerived().real_buffer[(x * getDimY() + y) * shift + z];
        }
    }

    template<class Derived>
    template<class GridType>
    inline void FFTRSpace<Derived, 3>::transform(const LValueGrid<GridType>& data) {
        assert(data.getDimX() == getDimX());
        assert(data.getDimY() == getDimY());
        assert(data.getDimZ() == getDimZ());
        *this = data;
        Base::getDerived().transform();
    }
}
