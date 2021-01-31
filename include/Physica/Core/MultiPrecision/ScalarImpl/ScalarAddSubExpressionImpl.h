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

//Forward declaration
namespace Physica::Core {
    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    Internal::ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator+(const Scalar<type, errorTrack1>& s1, const Scalar<type, errorTrack2>& s2);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    Internal::ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator-(const Scalar<type, errorTrack1>& s1, const Scalar<type, errorTrack2>& s2);
}

namespace Physica::Core::Internal {
    template<ScalarType type>
    template<bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpressionHelper<type>::ScalarAddSubExpressionHelper(
            const Scalar<type, errorTrack1>& s1, const Scalar<type, errorTrack2>& s2, bool minus) {
        constexpr NodeType nodeType1 = errorTrack1 ? PositiveWithError : PositiveWithoutError;
        constexpr NodeType nodeType2helper1 = errorTrack2 ? NegativeWithError : NegativeWithoutError;
        constexpr NodeType nodeType2helper2 = errorTrack2 ? PositiveWithError : PositiveWithoutError;
        NodeType nodeType2 = minus ? nodeType2helper1 : nodeType2helper2;
        terms.push_back(Node{&s1, nodeType1});
        terms.push_back(Node{&s2, nodeType2});
    }

    template<ScalarType type>
    template<bool errorTrack>
    ScalarAddSubExpressionHelper<type>::ScalarAddSubExpressionHelper(
            const Scalar<type, errorTrack>& s, ScalarAddSubExpressionHelper<type>&& exp, bool minus)
            : terms(std::move(exp.terms)) {
        constexpr NodeType typeToInsert = errorTrack ? PositiveWithError : PositiveWithoutError;
        auto ite = terms.begin();
        if(minus) {
            for(; ite != terms.end(); ++ite) {
                Node& node = *ite;
                changeNodeSign(node.nodeType);
                if(absCompare(*node.value, s)) {
                    ite = terms.insert(ite, Node{&s, typeToInsert});
                    ++ite;
                    goto inserted;
                }
            }
            terms.insert(ite, Node{&s, typeToInsert});
            return; //We have reached the end of the list.
        inserted:
            //Change the sign of the other nodes.
            for(; ite != terms.end(); ++ite) {
                Node& node = *ite;
                changeNodeSign(node.nodeType);
            }
        }
        else {
            for(; ite != terms.end(); ++ite) {
                Node& node = *ite;
                if(absCompare(*node.value, s)) {
                    terms.insert(ite, Node{&s, typeToInsert});
                    return;
                }
            }
            terms.insert(ite, Node{&s, typeToInsert});
        }
    }

    template<ScalarType type>
    template<bool errorTrack>
    ScalarAddSubExpressionHelper<type>::ScalarAddSubExpressionHelper(
            ScalarAddSubExpressionHelper<type>&& exp, const Scalar<type, errorTrack>& s, bool minus)
            : terms(std::move(exp.terms)) {
        constexpr NodeType typeToInsertHelper1 = errorTrack ? PositiveWithError : PositiveWithoutError;
        constexpr NodeType typeToInsertHelper2 = errorTrack ? NegativeWithError : NegativeWithoutError;
        NodeType typeToInsert = minus ? typeToInsertHelper2 : typeToInsertHelper1;

        auto ite = terms.begin();
        for(; ite != terms.end(); ++ite) {
            Node& node = *ite;
            if(absCompare(*node.value, s)) {
                terms.insert(ite, Node{&s, typeToInsert});
                return;
            }
        }
        terms.insert(ite, Node{&s, typeToInsert});
    }

    template<ScalarType type>
    ScalarAddSubExpressionHelper<type>::ScalarAddSubExpressionHelper(
            const ScalarAddSubExpressionHelper<type>& exp1, const ScalarAddSubExpressionHelper<type>& exp2, bool minus) {
        auto ite1 = exp1.terms.begin();
        auto ite2 = exp2.terms.begin();
        for (; ite1 != exp1.terms.end() && ite2 != exp2.terms.end(); ++ite1, ++ite2) {
            Node& node1 = *ite1;
            Node& node2 = *ite2;
            if(absCompare(*node1.value, *node2.value))
                terms.push_back(node2);
            else
                terms.push_back(node1);
        }
        if(ite1 == exp1.terms.end()){
            for(; ite2 != exp2.terms.end(); ++ite2)
                terms.push_back(*ite2);
        }
        else if(ite2 == exp1.terms.end()) {
            for(; ite1 != exp1.terms.end(); ++ite1)
                terms.push_back(*ite1);
        }
    }

    template<>
    Scalar<MultiPrecision, true> ScalarAddSubExpressionHelper<MultiPrecision>::calculate() const;

    template<ScalarType type>
    Scalar<type, true> ScalarAddSubExpressionHelper<type>::calculate() const {
        Scalar<type, true> result(0);
        bool withError = false; //Whether we have encountered a with-error scalar or not.
        for(auto node : terms) {
            const auto& value = *node.value;
            switch (node.nodeType) {
                case PositiveWithoutError:
                    if (withError)
                        result = Scalar<type, true>(result.getTrivial() + value.getTrivial(), result.getA());
                    else
                        result = Scalar<type, true>(result.getTrivial() + value.getTrivial());
                    break;
                case PositiveWithError: {
                    const auto& temp = static_cast<const Scalar<type, true>&>(value);
                    if (withError)
                        result = Scalar<type, true>(result.getTrivial() + value.getTrivial()
                                , result.getA() + temp.getA());
                    else
                        result = Scalar<type, true>(result.getTrivial() + value.getTrivial()
                                , temp.getA());
                    withError = true;
                }
                    break;
                case NegativeWithoutError:
                    if (withError)
                        result = Scalar<type, true>(result.getTrivial() - value.getTrivial(), result.getA());
                    else
                        result = Scalar<type, true>(result.getTrivial() - value.getTrivial());
                    break;
                case NegativeWithError:
                    const auto& temp = static_cast<const Scalar<type, true>&>(value);
                    if (withError)
                        result = Scalar<type, true>(result.getTrivial() - value.getTrivial()
                                , result.getA() + temp.getA());
                    else
                        result = Scalar<type, true>(result.getTrivial() - value.getTrivial()
                                , temp.getA());
                    withError = true;
                    break;
            }
        }
        return result;
    }
}
