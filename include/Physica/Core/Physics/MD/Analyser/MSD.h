/*
 * Copyright 2023 WeiBo He.
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

#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim>
    class MSD {
        using MDCellType = MDCell<ScalarType, Dim>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using VectorType = Vector<ScalarType, Dim>;

        MDCellType initCell;
        PositionMatrix buffer;
        size_t numSample;
    public:
        MSD(MDCellType initCell_);
        MSD(const MSD&) = default;
        MSD(MSD&&) noexcept = default;
        ~MSD() = default;
        /* Operators */
        MSD& operator=(MSD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void sample(const MDCellType& sample);
        template<size_t NumReplica, class ForceMatrixAllocator>
        inline void sample(const RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>& rpmd);

        [[nodiscard]] ScalarType calcMSD() const { return mean(calcAtomMSD()); }
        [[nodiscard]] ScalarType calcMSD2D() const { return mean(calcAtomMSD2D()); }
        [[nodiscard]] Vector<ScalarType> calcAtomMSD() const;
        [[nodiscard]] Vector<ScalarType> calcAtomMSD2D() const;
        [[nodiscard]] inline ScalarType calcFiniteSizeLimit(size_t atomId) const noexcept;
        void clear();
        void swap(MSD& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return initCell.getNumParticle(); }
        [[nodiscard]] size_t getNumSample() const noexcept { return numSample; }
    };

    template<class ScalarType, unsigned int Dim>
    MSD<ScalarType, Dim>::MSD(MDCellType initCell_) : initCell(std::move(initCell_)), numSample(0) {
        buffer.resize(getNumParticle(), 3);
        buffer = ScalarType(0);
    }

    template<class ScalarType, unsigned int Dim>
    void MSD<ScalarType, Dim>::sample(const MDCellType& sample) {
        const size_t numParticle = getNumParticle();
        const auto& pos = sample.getPos();
        for (size_t i = 0; i < numParticle; ++i) {
            const auto r = pos.row(i);
            const VectorType squaredDist = square(initCell.minDistVector(r, i));
            auto meanSquaredDist = buffer.row(i);
            toNextMean(meanSquaredDist, numSample, squaredDist);
        }
        numSample += 1;
    }

    template<class ScalarType, unsigned int Dim>
    template<size_t NumReplica, class ForceMatrixAllocator>
    inline void MSD<ScalarType, Dim>::sample(const RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>& rpmd) {
        for (size_t i = 0; i < rpmd.getNumReplica(); ++i)
            sample(rpmd.phaseToCell(i));
    }

    template<class ScalarType, unsigned int Dim>
    Vector<ScalarType> MSD<ScalarType, Dim>::calcAtomMSD() const {
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> result(numParticle);
        for (size_t i = 0; i < numParticle; ++i)
            result[i] = buffer.row(i).sum();
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    Vector<ScalarType> MSD<ScalarType, Dim>::calcAtomMSD2D() const {
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> result(numParticle);
        for (size_t i = 0; i < numParticle; ++i) {
            const auto row = buffer.row(i);
            result[i] = row[0] + row[1];
        }
        return result;
    }
    /**
     * \returns MSD cannot larger than this value because of finite size effect
     */
    template<class ScalarType, unsigned int Dim>
    inline ScalarType MSD<ScalarType, Dim>::calcFiniteSizeLimit(size_t atomId) const noexcept {
        return initCell.minDistVector(atomId, atomId).squaredNorm();
    }

    template<class ScalarType, unsigned int Dim>
    void MSD<ScalarType, Dim>::clear() {
        buffer = ScalarType(0);
        numSample = 0;
    }

    template<class ScalarType, unsigned int Dim>
    void MSD<ScalarType, Dim>::swap(MSD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        initCell.swap(obj.initCell);
        buffer.swap(obj.buffer);
        std::swap(numSample, obj.numSample);
    }
}
