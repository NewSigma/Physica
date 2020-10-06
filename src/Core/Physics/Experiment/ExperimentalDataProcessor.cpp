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
#include "Physica/Core/Physics/Experiment/ExperimentalDataProcessor.h"

namespace Physica::Core::Physics {
    ExperimentalDataProcessor& ExperimentalDataProcessor::operator=(const ExperimentalDataProcessor& processor) {
        if(this != &processor)
            data = processor.data;
        return *this;
    }

    ExperimentalDataProcessor& ExperimentalDataProcessor::operator=(ExperimentalDataProcessor&& processor) noexcept {
        data = std::move(processor.data);
        return *this;
    }

    void ExperimentalDataProcessor::compensate(const Scalar<MultiPrecision, false>& s) {
        for(auto& item : data)
            item[0] -= s;
    }

    void ExperimentalDataProcessor::updateInfo() {
        /* Calc total */ {
            Scalar<MultiPrecision, false> total(0);
            for(auto& item : data)
                total += item[0];
            info.total = total;
        }
        const auto column = data.getColumn();
        /* Calc average */
        info.average = info.total / Scalar<MultiPrecision, false>(static_cast<SignedScalarUnit>(column));
        /* Calc standardDeviation */ {
            Scalar<MultiPrecision, false> temp(0);
            for(auto& item : data)
                temp += square(item[0] - info.average);
            info.standardDeviation = sqrt(temp
                                          / Scalar<MultiPrecision, false>(static_cast<SignedScalarUnit>(column - 1)));
        }
        /* Calc averageStandardDeviation */ {
            info.averageStandardDeviation =
                    info.standardDeviation / Scalar<MultiPrecision, false>(static_cast<SignedScalarUnit>(column));
        }
    }
}