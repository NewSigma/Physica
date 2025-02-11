/*
 * Copyright 2022-2024 Weibo He.
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

#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica {
    /**
     * Reference:
     * [1] Daan Frenkel and Berend Smit. Understanding Molecular Simulation[M]. Academic Press. 2002:85-86
     */
    template<Scalar T>
    class RDF {
        constexpr static unsigned int Dim = 3;
    private:
        Array<bool> isFromParticle;
        Array<bool> isToParticle;
        Array<uint64_t> particleBucket;
        T cellVolume;
        T stepSize;
        unsigned int numStep;
        unsigned int numSample;
    public:
        RDF() = default;
        RDF(Array<bool> isFromParticle_,
            Array<bool> isToParticle_,
            T cellVolume_,
            T stepSize_,
            unsigned int numStep_);
        RDF(const RDF&) = default;
        RDF(RDF&&) noexcept = default;
        ~RDF() = default;
        /* Operators */
        RDF& operator=(RDF obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Scalar U>
        void sample(const MDCell<U>& cell);
        VectorND<T> makeDists() const;
        VectorND<T> makeRDF() const;
        void swap(RDF& __restrict rdf) noexcept;
        [[nodiscard]] bool checkRadius(const PeriodicCell<T, 3>& cell) const;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return isFromParticle.getLength(); }
        [[nodiscard]] T getCellVolume() const noexcept { return cellVolume; }
        [[nodiscard]] T getMaxRadius() const noexcept { return stepSize * (numStep + 1); }
        /* Setters */
        void setCellVolume(T cellVolume_) { cellVolume = cellVolume_; }
    private:
        uint64_t sumFromParticle() const;
        uint64_t sumToParticle() const;
    };

    template<Scalar T>
    RDF<T>::RDF(Array<bool> isFromParticle_,
                         Array<bool> isToParticle_,
                         T cellVolume_,
                         T stepSize_,
                         unsigned int numStep_)
            : isFromParticle(std::move(isFromParticle_))
            , isToParticle(std::move(isToParticle_))
            , particleBucket(numStep_, 0)
            , cellVolume(cellVolume_)
            , stepSize(std::move(stepSize_))
            , numStep(numStep_)
            , numSample(0) {
        assert(isFromParticle.getLength() == isToParticle.getLength());
        assert(stepSize.isPositive());
    }
    /**
     * Optimize: If isFromParticle == isToParticle, we may make use of symmetry
     */
    template<Scalar T>
    template<Scalar U>
    void RDF<T>::sample(const MDCell<U>& cell) {
        using CellType = MDCell<U>;
        using VectorType = DenseVector<T, Dim>;

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
                    const T r = (from_pos - cell.getPos().row(to)).norm();
                    const uint64_t index = double(r / stepSize);
                    if (0 < index && index <= numStep)
                        ++particleBucket[index - 1];
                }
            }
        });
    }

    template<Scalar T>
    VectorND<T> RDF<T>::makeDists() const {
        VectorND<T> dists(numStep);
        for (size_t i = 0; i < dists.getLength(); ++i)
            dists[i] = stepSize * (i + 1.5);
        return dists;
    }

    template<Scalar T>
    VectorND<T> RDF<T>::makeRDF() const {
        VectorND<T> rdf(numStep);
        const uint64_t numFromParticle = sumFromParticle();
        const uint64_t numToParticle = sumToParticle();
        for (size_t i = 0; i < numStep; ++i) {
            const T temp = (stepSize * stepSize * stepSize) * ((i + 2) * (i + 2) * (i + 2) - (i + 1) * (i + 1) * (i + 1));
            const T numParticleIdealGas = T(4.0 / 3 * M_PI) * temp / cellVolume * numToParticle;
            const T factor = reciprocal(numParticleIdealGas * (numFromParticle * numSample));
            rdf[i] = factor * particleBucket[i];
        }
        return rdf;
    }

    template<Scalar T>
    void RDF<T>::swap(RDF& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        isFromParticle.swap(obj.isFromParticle);
        isToParticle.swap(obj.isToParticle);
        particleBucket.swap(obj.particleBucket);
        cellVolume.swap(obj.cellVolume);
        stepSize.swap(obj.stepSize);
        std::swap(numStep, obj.numStep);
        std::swap(numSample, obj.numSample);
    }

    template<Scalar T>
    bool RDF<T>::checkRadius(const PeriodicCell<T, 3>& cell) const {
        const T radius = getMaxRadius();
        const bool isBadRadius = cell.minDistVector(0, 0).norm() < radius;
        return isBadRadius;
    }

    template<Scalar T>
    uint64_t RDF<T>::sumFromParticle() const {
        uint64_t sum = 0;
        for (bool elem : isFromParticle)
            sum += elem;
        return sum;
    }

    template<Scalar T>
    uint64_t RDF<T>::sumToParticle() const {
        uint64_t sum = 0;
        for (bool elem : isToParticle)
            sum += elem;
        return sum;
    }
}
