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

#include <cassert>

namespace Physica::Core {
    template<class Derived>
    template<class OtherDerived>
    void RValueMatrix<Derived>::assignTo(LValueMatrix<OtherDerived>& target) const {
        assert(getRow() == target.getRow() && getColumn() == target.getColumn());
        Base::getDerived().assignTo(target);
    }

    template<class Derived>
    void RValueMatrix<Derived>::toUnitMatrix() {
        assert(getRow() == getColumn());
        auto matIterator = Base::getDerived().begin();
        auto eleIterator = Base::getDerived().ebegin(matIterator);
        const size_t order = getRow();
        for (size_t i = 0; i < order; ++i) {
            for (size_t j = 0; j < order; ++j) {
                *eleIterator = i == j;
                ++eleIterator;
            }
            Derived::updateIterator(matIterator, eleIterator);
        }
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class Derived>
    Derived reciprocal(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(reciprocal(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived sqrt(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sqrt(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived factorial(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.allocate(factorial(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived ln(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(ln(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived log(const RValueMatrix<Derived>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(log(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived exp(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(exp(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived pow(const RValueMatrix<Derived>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(pow(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cos(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cos(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sin(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sin(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived tan(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(tan(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sec(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sec(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived csc(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(csc(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cot(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cot(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccos(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccos(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsin(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsin(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arctan(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arctan(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsec(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsec(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccsc(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccsc(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccot(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccot(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cosh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cosh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sinh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sinh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived tanh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(tanh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sech(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sech(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived csch(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(csch(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived coth(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(coth(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccosh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccosh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsinh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsinh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arctanh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arctanh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsech(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsech(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccsch(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccsch(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccoth(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccoth(m[i]), i);
        result.setLength(length);
        return result;
    }
}