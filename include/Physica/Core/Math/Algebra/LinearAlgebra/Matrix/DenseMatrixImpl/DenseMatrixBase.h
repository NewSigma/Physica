/*
 * Copyright 2021 WeiBo He.
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

#include "DenseMatrixStorage/DenseMatrixStorage.h"

namespace Physica::Core {
    /**
     * The \class DenseDenseMatrixBase provide algorithms that a matrix should support.
     * 
     * TODO: merge DenseMatrix into DenseMatrixBase
     * 
     * \tparam Derived
     * A class that contains data structure for a matrix.
     */
    template<class Derived>
    class DenseMatrixBase : public DenseMatrixStorage<typename Internal::Traits<Derived>::ScalarType
                                                      , Internal::Traits<Derived>::MatrixType
                                                      , Internal::Traits<Derived>::RowAtCompile
                                                      , Internal::Traits<Derived>::ColumnAtCompile
                                                      , Internal::Traits<Derived>::MaxRowAtCompile
                                                      , Internal::Traits<Derived>::MaxColumnAtCompile>
                            , public Utils::CRTPBase<Derived> {
        using Base = DenseMatrixStorage<typename Internal::Traits<Derived>::ScalarType
                                        , Internal::Traits<Derived>::MatrixType
                                        , Internal::Traits<Derived>::RowAtCompile
                                        , Internal::Traits<Derived>::ColumnAtCompile
                                        , Internal::Traits<Derived>::MaxRowAtCompile
                                        , Internal::Traits<Derived>::MaxColumnAtCompile>;
        using Utils::CRTPBase<Derived>::getDerived;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
    public:
        using Base::Base;
        /* Operators */
        template<class OtherDerived>
        inline void operator+=(const DenseMatrixBase<OtherDerived>& m) { getDerived() = *this + m; }
        template<class AnyScalar>
        inline void operator+=(const ScalarBase<AnyScalar>& s) { getDerived() = *this + s.getDerived(); }
        template<class OtherDerived>
        inline void operator-=(const DenseMatrixBase<OtherDerived>& m) { getDerived() = *this - m; }
        template<class AnyScalar>
        inline void operator-=(const ScalarBase<AnyScalar>& s) { getDerived() = *this - s.getDerived(); }
        template<class OtherDerived>
        inline void operator*=(const DenseMatrixBase<OtherDerived>& m) { getDerived() = *this * m; }
        template<class AnyScalar>
        inline void operator*=(const ScalarBase<AnyScalar>& s) { getDerived() = *this * s.getDerived(); }
        /* Operations */
        template<class OtherDerived>
        void assignTo(DenseMatrixBase<OtherDerived>& mat) const;
        ScalarType determinate() const;
        void rowReduce(size_t r1, size_t r2, size_t elementIndex);
        void rowReduce(size_t r1, size_t r2, const ScalarType& factor);
        void columnReduce(size_t c1, size_t c2, size_t elementIndex);
        void columnReduce(size_t c1, size_t c2, const ScalarType& factor);
        inline void majorReduce(size_t v1, size_t v2, size_t elementIndex);
        inline void majorReduce(size_t v1, size_t v2, const ScalarType& factor);
        void rowMulScalar(size_t r, const ScalarType& factor);
        void columnMulScalar(size_t c, const ScalarType& factor);
        inline void majorMulScalar(size_t v, const ScalarType& factor);
        inline void majorSwap(size_t v1, size_t v2);
        /* Getters */
        [[nodiscard]] ScalarType& getElementFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] const ScalarType& getElementFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] inline size_t getOrder() const noexcept;
        [[nodiscard]] inline size_t getMaxMajor() const noexcept;
        [[nodiscard]] inline size_t getMaxMinor() const noexcept;
        /* Setters */
        void toUnitMatrix();
        /* Static members */
        [[nodiscard]] inline size_t rowFromMajorMinor(size_t major, size_t minor) const noexcept;
        [[nodiscard]] inline size_t columnFromMajorMinor(size_t major, size_t minor) const noexcept;
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class Derived>
    Derived reciprocal(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived sqrt(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived factorial(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived ln(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived log(const DenseMatrixBase<Derived>& m, const MultiScalar& a);
    
    template<class Derived>
    Derived exp(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived pow(const DenseMatrixBase<Derived>& m, const MultiScalar& a);
    
    template<class Derived>
    Derived cos(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived sin(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived tan(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived sec(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived csc(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived cot(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arccos(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arcsin(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arctan(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arcsec(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arccsc(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arccot(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived cosh(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived sinh(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived tanh(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived sech(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived csch(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived coth(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arccosh(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arcsinh(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arctanh(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arcsech(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arccsch(const DenseMatrixBase<Derived>& m);
    
    template<class Derived>
    Derived arccoth(const DenseMatrixBase<Derived>& m);
}

#include "DenseMatrixBaseImpl.h"