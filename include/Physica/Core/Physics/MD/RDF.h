/*
 * Copyright 2022 WeiBo He.
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

#include "MDCell.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Daan Frenkel and Berend Smit. Understanding Molecular Simulation[M]. Academic Press. 2002:85-86
     */
    template<class ScalarType>
    class RDF {
        constexpr static unsigned int Dim = 3;
    private:
        Utils::Array<bool> isFromParticle;
        Utils::Array<bool> isToParticle;
        Utils::Array<uint64_t> particleBucket;
        ScalarType stepSize;
        unsigned int numStep;
        unsigned int numSample;
    public:
        RDF() = default;
        RDF(Utils::Array<bool> isFromParticle_, Utils::Array<bool> isToParticle_, ScalarType stepSize_, unsigned int numStep_);
        RDF(const RDF&) = default;
        RDF(RDF&&) noexcept = default;
        ~RDF() = default;
        /* Operators */
        RDF& operator=(RDF rdf) noexcept;
        /* Operations */
        template<class T>
        void sample(const MDCell<T, ScalarType>& cell);
        Vector<ScalarType> makeDists() const;
        Vector<ScalarType> makeRDF(ScalarType volume) const;
        void swap(RDF& rdf) noexcept;
        [[nodiscard]] bool checkRadius(const PeriodicCell<ScalarType, 3>& cell) const;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return isFromParticle.getLength(); }
        [[nodiscard]] ScalarType getMaxRadius() const noexcept { return stepSize * (numStep + 1); }
    private:
        uint64_t sumFromParticle() const;
        uint64_t sumToParticle() const;
    };

    template<class ScalarType>
    RDF<ScalarType>::RDF(Utils::Array<bool> isFromParticle_,
                         Utils::Array<bool> isToParticle_,
                         ScalarType stepSize_,
                         unsigned int numStep_)
            : isFromParticle(std::move(isFromParticle_))
            , isToParticle(std::move(isToParticle_))
            , particleBucket(numStep_, 0)
            , stepSize(std::move(stepSize_))
            , numStep(numStep_)
            , numSample(0) {
        assert(isFromParticle.getLength() == isToParticle.getLength());
        assert(stepSize.isPositive());
    }

    template<class ScalarType>
    RDF<ScalarType>& RDF<ScalarType>::operator=(RDF rdf) noexcept {
        swap(rdf);
        return *this;
    }
    /**
     * Optimize: If isFromParticle == isToParticle, we may make use of symmetry
     */
    template<class ScalarType>
    template<class T>
    void RDF<ScalarType>::sample(const MDCell<T, ScalarType>& cell) {
        using CellType = MDCell<T, ScalarType>;
        using VectorType = Vector<ScalarType, Dim>;

        assert(getNumParticle() == cell.getNumParticle());
        const auto range = CellType::estimateRange(cell.getLattice(), getMaxRadius());
        numSample++;
        CellType::forCellInRange(range, cell.getLattice(), [this, &cell](VectorType delta) {
            for (size_t from = 0; from < getNumParticle(); ++from) {
                if (!isFromParticle[from])
                    continue;
                const VectorType from_pos = cell.getPos().row(from) + delta;
                for (size_t to = 0; to < getNumParticle(); ++to) {
                    if (!isToParticle[to])
                        continue;
                    const ScalarType r = (from_pos - cell.getPos().row(to)).norm();
                    const uint64_t index = double(r / stepSize);
                    if (0 < index && index <= numStep)
                        ++particleBucket[index - 1];
                }
            }
        });
    }

    template<class ScalarType>
    Vector<ScalarType> RDF<ScalarType>::makeDists() const {
        Vector<ScalarType> dists(numStep);
        for (size_t i = 0; i < dists.getLength(); ++i)
            dists[i] = stepSize * (i + 1.5);
        return dists;
    }

    template<class ScalarType>
    Vector<ScalarType> RDF<ScalarType>::makeRDF(ScalarType volume) const {
        Vector<ScalarType> rdf(numStep);
        const uint64_t numFromParticle = sumFromParticle();
        const uint64_t numToParticle = sumToParticle();
        for (size_t i = 0; i < numStep; ++i) {
            const ScalarType temp = (stepSize * stepSize * stepSize) * ((i + 2) * (i + 2) * (i + 2) - (i + 1) * (i + 1) * (i + 1));
            const ScalarType numParticleIdealGas = ScalarType(4.0 / 3 * M_PI) * temp / volume * numToParticle;
            const ScalarType factor = reciprocal(numParticleIdealGas * (numFromParticle * numSample));
            rdf[i] = factor * particleBucket[i];
        }
        return rdf;
    }

    template<class ScalarType>
    void RDF<ScalarType>::swap(RDF& rdf) noexcept {
        assert(this != &rdf && "[Error]: Self swap is likely a bug");
        isFromParticle.swap(rdf.isFromParticle);
        isToParticle.swap(rdf.isToParticle);
        particleBucket.swap(rdf.particleBucket);
        stepSize.swap(rdf.stepSize);
        std::swap(numStep, rdf.numStep);
        std::swap(numSample, rdf.numSample);
    }

    template<class ScalarType>
    bool RDF<ScalarType>::checkRadius(const PeriodicCell<ScalarType, 3>& cell) const {
        const ScalarType radius = getMaxRadius();
        const bool isBadRadius = cell.minDistVector(0, 0).norm() < radius;
        return isBadRadius;
    }

    template<class ScalarType>
    uint64_t RDF<ScalarType>::sumFromParticle() const {
        uint64_t sum = 0;
        for (bool elem : isFromParticle)
            sum += elem;
        return sum;
    }

    template<class ScalarType>
    uint64_t RDF<ScalarType>::sumToParticle() const {
        uint64_t sum = 0;
        for (bool elem : isToParticle)
            sum += elem;
        return sum;
    }
}
