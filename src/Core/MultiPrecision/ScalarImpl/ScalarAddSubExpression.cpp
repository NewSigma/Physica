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
#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core::Intenal {
    template<>
    Scalar<MultiPrecision, true> ScalarAddSubExpressionHelper<MultiPrecision>::calculate() const {
        Scalar<MultiPrecision, true> result(0);
        bool withError = false; //Whether we have encountered a with-error scalar or not.
        for(auto node : terms) {
            const auto& value = *node.value;
            switch (node.nodeType) {
                case PositiveWithoutError:
                    if (withError) {
                        result = Scalar<MultiPrecision, false>::addWithError(result, value);
                        Scalar<MultiPrecision, false>::cutLength(result);
                        if(result.getA() != 0)
                            result.applyError(result.getAccuracy());
                    }
                    else {
                        result = Scalar<MultiPrecision, false>::addNoError(result, value);
                        Scalar<MultiPrecision, false>::cutLength(result);
                    }
                    break;
                case PositiveWithError: {
                    result = Scalar<MultiPrecision, false>::addWithError(result, value);
                    Scalar<MultiPrecision, false>::cutLength(result);
                    const auto& temp = static_cast<const Scalar<MultiPrecision, true>&>(value); //NOLINT Safe cast
                    if (withError) {
                        if(result.getA() != 0 || temp.getA() != 0) {
                            if(result.getA() == 0)
                                result.applyError(temp.getAccuracy());
                            else if(temp.getA() == 0)
                                result.applyError(result.getAccuracy());
                            else
                                result.applyError(result.getAccuracy() + temp.getAccuracy());
                        }
                    }
                    else {
                        if(temp.getA() != 0)
                            result.applyError(temp.getAccuracy());
                    }
                    withError = true;
                }
                    break;
                case NegativeWithoutError:
                    if (withError) {
                        result = Scalar<MultiPrecision, false>::subWithError(result, value);
                        Scalar<MultiPrecision, false>::cutLength(result);
                        if(result.getA() != 0)
                            result.applyError(result.getAccuracy());
                    }
                    else {
                        result = Scalar<MultiPrecision, false>::subNoError(result, value);
                        Scalar<MultiPrecision, false>::cutLength(result);
                    }
                    break;
                case NegativeWithError:
                    result = Scalar<MultiPrecision, false>::subWithError(result, value);
                    Scalar<MultiPrecision, false>::cutLength(result);
                    const auto& temp = static_cast<const Scalar<MultiPrecision, true>&>(value); //NOLINT Safe cast
                    if (withError) {
                        if(result.getA() != 0 || temp.getA() != 0) {
                            if(result.getA() == 0)
                                result.applyError(temp.getAccuracy());
                            else if(temp.getA() == 0)
                                result.applyError(result.getAccuracy());
                            else
                                result.applyError(result.getAccuracy() + temp.getAccuracy());
                        }
                    }
                    else {
                        if(temp.getA() != 0)
                            result.applyError(temp.getAccuracy());
                    }
                    withError = true;
                    break;
            }
        }
        return result;
    }
}