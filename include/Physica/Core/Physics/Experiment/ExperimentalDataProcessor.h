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
#ifndef PHYSICA_EXPERIMENTALDATAPROCESSOR_H
#define PHYSICA_EXPERIMENTALDATAPROCESSOR_H

#include <utility>
#include <Physica/Core/MultiPrecision/Scalar.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>

namespace Physica::Core::Physics {
    /*!
     * This class handles experimental data, including identify bad data and calculate some useful values.
     *
     * Number of data cannot be more than maximum of @typedef MPUnit.
     */
    class ExperimentalDataProcessor {
        using ScalarType = Scalar<Double, false>;
        using DataMatrix = DenseMatrix<ScalarType>;
    public:
        struct ExperimentalDataInfo {
            ScalarType total;
            ScalarType average;
            ScalarType standardDeviation;
            ScalarType averageStandardDeviation;
        };
    private:
        /*
         * Use column matrix for the convenience of add and remove data.
         *
         * The first row of the matrix contains the measured values.
         * The second row of the matrix contains uncertainty of the measured values.
         */
        DataMatrix data;
        ExperimentalDataInfo info;
    public:
        inline explicit ExperimentalDataProcessor(DataMatrix m);
        ~ExperimentalDataProcessor() = default;
        ExperimentalDataProcessor(const ExperimentalDataProcessor& processor) = default;
        ExperimentalDataProcessor(ExperimentalDataProcessor&& processor)  noexcept : data(std::move(processor.data)) {}
        /* Operators */
        ExperimentalDataProcessor& operator=(const ExperimentalDataProcessor& processor);
        ExperimentalDataProcessor& operator=(ExperimentalDataProcessor&& processor) noexcept;
        /* Operations */
        void compensate(const ScalarType& s);
        void updateInfo();
        /* Getters */
        const ExperimentalDataInfo& getInfo() noexcept { return info; }
    private:
    };

    inline ExperimentalDataProcessor::ExperimentalDataProcessor(DataMatrix m) : data(std::move(m)) {
        assert(data.getColumn() < MPUnitMax); //Too much data is not supported.
        updateInfo();
    }
}

#endif
