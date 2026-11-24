/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"

namespace Physica {
    template<Scalar T, unsigned int Dim>
    class MSD {
        using MDCellType = MDCell<T, Dim>;
        using PositionMatrix = MDCellType::PositionMatrix;
        using VectorType = DenseVector<T, Dim>;

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
        void sample(const RPMD<T, Dim, NumReplica, ForceMatrixAllocator>& rpmd);

        [[nodiscard]] T calcMSD() const { return calcAtomMSD().mean(); }
        [[nodiscard]] T calcMSD2D() const { return calcAtomMSD2D().mean(); }
        [[nodiscard]] VectorND<T> calcAtomMSD() const;
        [[nodiscard]] VectorND<T> calcAtomMSD2D() const;
        [[nodiscard]] T calcFiniteSizeLimit(size_t atomId) const noexcept;
        void clear();
        void swap(MSD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return initCell.getNumParticle(); }
        [[nodiscard]] size_t getNumSample() const noexcept { return numSample; }
    };

    template<Scalar T, unsigned int Dim>
    MSD<T, Dim>::MSD(MDCellType initCell_) : initCell(std::move(initCell_)), numSample(0) {
        buffer.resize(getNumParticle(), 3);
        buffer.zeros();
    }

    template<Scalar T, unsigned int Dim>
    void MSD<T, Dim>::sample(const MDCellType& sample) {
        const size_t numParticle = getNumParticle();
        const auto& pos = sample.getPos();
        for (size_t i = 0; i < numParticle; ++i) {
            const auto r = pos.row(i);
            const VectorType squaredDist = square(initCell.minDistVector(r, i));
            auto meanSquaredDist = buffer.row(i);
            meanSquaredDist.toNextMean(numSample, squaredDist);
        }
        numSample += 1;
    }

    template<Scalar T, unsigned int Dim>
    template<size_t NumReplica, class ForceMatrixAllocator>
    void MSD<T, Dim>::sample(const RPMD<T, Dim, NumReplica, ForceMatrixAllocator>& rpmd) {
        for (size_t i = 0; i < rpmd.getNumReplica(); ++i)
            sample(rpmd.phaseToCell(i));
    }

    template<Scalar T, unsigned int Dim>
    VectorND<T> MSD<T, Dim>::calcAtomMSD() const {
        const size_t numParticle = getNumParticle();
        VectorND<T> result(numParticle);
        for (size_t i = 0; i < numParticle; ++i)
            result[i] = buffer.row(i).sum();
        return result;
    }

    template<Scalar T, unsigned int Dim>
    VectorND<T> MSD<T, Dim>::calcAtomMSD2D() const {
        const size_t numParticle = getNumParticle();
        VectorND<T> result(numParticle);
        for (size_t i = 0; i < numParticle; ++i) {
            const auto row = buffer.row(i);
            result[i] = row[0] + row[1];
        }
        return result;
    }
    /**
     * \returns MSD cannot larger than this value because of finite size effect
     */
    template<Scalar T, unsigned int Dim>
    T MSD<T, Dim>::calcFiniteSizeLimit(size_t atomId) const noexcept {
        return initCell.minDistVector(atomId, atomId).squaredNorm();
    }

    template<Scalar T, unsigned int Dim>
    void MSD<T, Dim>::clear() {
        buffer.zeros();
        numSample = 0;
    }

    template<Scalar T, unsigned int Dim>
    void MSD<T, Dim>::swap(MSD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        initCell.swap(obj.initCell);
        buffer.swap(obj.buffer);
        std::swap(numSample, obj.numSample);
    }
}
