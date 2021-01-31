/*
 * Copyright 2020 WeiBo He.
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

#include <list>
#include "Physica/Core/MultiPrecision/MultiPrecisionType.h"
/**
 * This file contains implementation of Scalar add and sub expressions,
 * which is created to provide support for associative law of addition of float numbers.
 */
namespace Physica::Core {
    template<ScalarType type, bool errorTrack> class Scalar;
}


namespace Physica::Core::Internal {
    class AbstractScalarAddSubExpression {
    public:
        /**
         * Format: (Operation between current node and the node before it) + (errorTrack).
         * e.g.
         * AddWithError stands for the current node carries its accuracy and it should add with the value of the node before it.
         */
        enum NodeType {
            PositiveWithoutError,
            PositiveWithError,
            NegativeWithoutError,
            NegativeWithError
        };
    protected:
        inline static void changeNodeSign(NodeType& type) {
            switch (type) {
                case PositiveWithoutError:
                    type = NegativeWithoutError;
                    break;
                case PositiveWithError:
                    type = NegativeWithError;
                    break;
                case NegativeWithoutError:
                    type = PositiveWithoutError;
                    break;
                case NegativeWithError:
                    type = PositiveWithError;
                    break;
            }
        }
    };

    template<ScalarType type>
    class ScalarAddSubExpressionHelper : public AbstractScalarAddSubExpression {
    private:
        struct Node {
            const Scalar<type, false>* value;
            NodeType nodeType;
        };
    private:
        std::list<Node> terms;
    protected:
        template<bool errorTrack1, bool errorTrack2>
        ScalarAddSubExpressionHelper(const Scalar<type, errorTrack1>& s1, const Scalar<type, errorTrack2>& s2, bool minus);
        template<bool errorTrack>
        ScalarAddSubExpressionHelper(const Scalar<type, errorTrack>& s, ScalarAddSubExpressionHelper<type>&& exp, bool minus);
        template<bool errorTrack>
        ScalarAddSubExpressionHelper(ScalarAddSubExpressionHelper<type>&& exp, const Scalar<type, errorTrack>& s, bool minus);
        ScalarAddSubExpressionHelper(const ScalarAddSubExpressionHelper<type>& exp1, const ScalarAddSubExpressionHelper<type>& exp2, bool minus);
        /* Operations */
        [[nodiscard]] Scalar<type, true> calculate() const;
    };
    /**
     * Optimize: errorTracks of scalars are known at compile, possible make use it and speed up.
     */
    template<ScalarType type, bool errorTrack>
    class ScalarAddSubExpression : public ScalarAddSubExpressionHelper<type> {
        using Base = ScalarAddSubExpressionHelper<type>;
    public:
        template<bool errorTrack1, bool errorTrack2>
        ScalarAddSubExpression(const Scalar<type, errorTrack1>& s1, const Scalar<type, errorTrack2>& s2, bool minus)
                : Base(s1, s2, minus) {}
        template<bool errorTrack1, bool errorTrack2>
        ScalarAddSubExpression(const Scalar<type, errorTrack1>& s, ScalarAddSubExpression<type, errorTrack2>&& exp, bool minus)
                : Base(s, std::move(exp), minus) {}
        template<bool errorTrack1, bool errorTrack2>
        ScalarAddSubExpression(ScalarAddSubExpression<type, errorTrack1>&& exp, const Scalar<type, errorTrack2>& s, bool minus)
                : Base(std::move(exp), s, minus) {}
        template<bool errorTrack1, bool errorTrack2>
        ScalarAddSubExpression(const ScalarAddSubExpression<type, errorTrack1>& exp1, const ScalarAddSubExpression<type, errorTrack2>& exp2, bool minus)
                : Base(exp1, exp2, minus) {}
        /* Operators */
        Scalar<type, errorTrack> operator<<(int bits) { return operator Scalar<type, errorTrack>() << bits; }
        Scalar<type, errorTrack> operator>>(int bits) { return operator Scalar<type, errorTrack>() >> bits; }
        Scalar<type, errorTrack> operator-() { return -operator Scalar<type, errorTrack>(); }
        /* Operators */
        operator Scalar<type, errorTrack>() const { return Base::calculate(); } //NOLINT Safe cast
        explicit operator double() const { return double(Base::calculate()); }
    };

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator+(ScalarAddSubExpression<type, errorTrack1>&& exp, const Scalar<type, errorTrack2>& s) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(std::move(exp), s, false);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator-(ScalarAddSubExpression<type, errorTrack1>&& exp, const Scalar<type, errorTrack2>& s) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(std::move(exp), s, true);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    Scalar<type, errorTrack1 || errorTrack2>
    inline operator*(const ScalarAddSubExpression<type, errorTrack1>& exp, const Scalar<type, errorTrack2>& s) {
        return Scalar<type, errorTrack1>(exp) * s;
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    Scalar<type, errorTrack1 || errorTrack2>
    inline operator/(const ScalarAddSubExpression<type, errorTrack1>& exp, const Scalar<type, errorTrack2>& s) {
        return Scalar<type, errorTrack1>(exp) / s;
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator+(ScalarAddSubExpression<type, errorTrack1>&& exp1, ScalarAddSubExpression<type, errorTrack2>&& exp2) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(std::move(exp1), std::move(exp2), false);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator-(ScalarAddSubExpression<type, errorTrack1>&& exp1, ScalarAddSubExpression<type, errorTrack2>&& exp2) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(std::move(exp1), std::move(exp2), true);
    }
}

#include "ScalarAddSubExpressionImpl.h"
