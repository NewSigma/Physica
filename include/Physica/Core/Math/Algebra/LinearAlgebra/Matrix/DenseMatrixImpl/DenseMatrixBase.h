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
     * The \class DenseDenseMatrixBase provide algorithms that a matrix should support, but the data is provided
     * by \tparam Derived.
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
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        ScalarType determinate() const;
        void rowReduce(size_t r1, size_t r2, size_t elementIndex);
        /* Getters */
        [[nodiscard]] inline size_t getOrder() const noexcept;
    };
    /* Operators */
    template<class Derived>
    Derived operator+(const DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2);

    template<class Derived>
    Derived operator-(const DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2);

    template<class Derived>
    Derived operator*(const DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2);

    template<class Derived>
    Derived operator*(const DenseMatrixBase<Derived>& m, const MultiScalar& n);

    template<class Derived>
    Derived operator-(const DenseMatrixBase<Derived>& m);
    /* Inline Implementations */
    template<class Derived>
    inline void operator+=(DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2) { m1 = m1 + m2; }

    template<class Derived>
    inline void operator-=(DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2) { m1 = m1 - m2; }

    template<class Derived>
    inline void operator*=(DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2) { m1 = m1 * m2; }

    template<class Derived>
    inline void operator*=(DenseMatrixBase<Derived>& m, const MultiScalar& n) { m = m * n; }
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