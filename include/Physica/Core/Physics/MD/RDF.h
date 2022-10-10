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
    template<class ScalarType, class PosScalarType>
    class RDF {
        using CellType = PeriodicCell<PosScalarType, 3>;

        Utils::Array<bool> isFromParticle;
        Utils::Array<bool> isToParticle;
        Utils::Array<uint64_t> particleBucket;
        PosScalarType stepSize;
        unsigned int numStep;
        unsigned int numSample;
    public:
        RDF();
        RDF(const RDF&) = default;
        RDF(RDF&&) noexcept = default;
        ~RDF() = default;
        /* Operators */
        RDF& operator=(RDF rdf) noexcept;
        /* Operations */
        void sample(const CellType& cell);
        Vector<PosScalarType> makeDists() const;
        Vector<PosScalarType> makeRDF(PosScalarType volume) const;
        void swap(RDF& rdf) noexcept;
        [[nodiscard]] bool checkRadius(const CellType& cell) const;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return isFromParticle.getLength(); }
        [[nodiscard]] PosScalarType getMaxRadius() const noexcept { return stepSize * (numStep + 1); }
    private:
        uint64_t sumParticleBucket() const;
    };

    template<class ScalarType, class PosScalarType>
    RDF<ScalarType, PosScalarType>::RDF(Utils::Array<bool> isFromParticle_,
                                        Utils::Array<bool> isToParticle_,
                                        PosScalarType stepSize_,
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

    template<class ScalarType, class PosScalarType>
    RDF<ScalarType, PosScalarType>& RDF<ScalarType, PosScalarType>::operator=(RDF rdf) noexcept {
        swap(rdf);
        return *this;
    }
    /**
     * Optimize: If isFromParticle == isToParticle, we may make use of symmetry
     */
    template<class ScalarType, class PosScalarType>
    void RDF<ScalarType, PosScalarType>::sample(const CellType& cell) {
        using VectorType = Vector<PosScalarType, CellType::Dim>;

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
                    const PosScalarType r2 = (from_pos - cell.getPos().row(to)).squaredNorm();
                    const uint64_t index = double(r2 / square(stepSize));
                    if (0 < index && index <= numStep)
                        ++particleBucket[index - 1];
                }
            }
        });
    }

    template<class ScalarType, class PosScalarType>
    Vector<PosScalarType> RDF<ScalarType, PosScalarType>::makeDists() const {
        Vector<PosScalarType> dists(numStep);
        for (size_t i = 0; i < dists.getLength(); ++i)
            dists[i] = stepSize * (i + 1.5);
        return dists;
    }

    template<class ScalarType, class PosScalarType>
    Vector<ScalarType> RDF<ScalarType, PosScalarType>::makeRDF(PosScalarType volume) const {
        Vector<ScalarType> rdf;
        for (size_t i = 0; i < numStep; ++i) {
            const PosScalarType temp = (numStep * numStep * numStep) * ((i + 2) * (i + 2) * (i + 2) - (i + 1) * (i + 1) * (i + 1));
            const PosScalarType numParticleIdealGas = (4.0 / 3 * M_PI) * temp / volume * getNumParticle();
            const PosScalarType factor = reciprocal(numParticleIdealGas * (numStep * numSample));
            rdf[i] = factor * particleBucket[i];
        }
        return rdf;
    }

    template<class ScalarType, class PosScalarType>
    void RDF<ScalarType, PosScalarType>::swap(RDF& rdf) noexcept {
        isFromParticle.swap(rdf.isFromParticle);
        isToParticle.swap(rdf.isToParticle);
        particleBucket.swap(rdf.particleBucket);
        stepSize.swap(rdf.stepSize);
        std::swap(numStep, rdf.numStep);
        std::swap(numSample, rdf.numSample);
    }

    template<class ScalarType, class PosScalarType>
    bool RDF<ScalarType, PosScalarType>::checkRadius(const CellType& cell) const {
        const PosScalarType radius = getMaxRadius();
        const bool isBadRadius = cell.minDistVector(0, 0).norm() < radius;
        return isBadRadius;
    }

    template<class ScalarType, class PosScalarType>
    uint64_t RDF<ScalarType, PosScalarType>::sumParticleBucket() const {
        uint64_t sum = 0;
        for (uint64_t elem : particleBucket)
            sum += elem;
        return sum;
    }
}
