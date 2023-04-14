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

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    class RingPolymer {
    public:
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using PhaseMatrix = DenseMatrix<PosScalarType, MatrixOption::Row | MatrixOption::Vector>;
        using BufferType = DenseMatrix<ComplexScalar<PosScalarType>, MatrixOption::Row | MatrixOption::Vector, 2>;
    private:
        PhaseMatrix phase;
        FFT<PosScalarType, 1> fft;
        BufferType buffer;
    public:
        RingPolymer() = default;
        RingPolymer(const MDCellType& cell, size_t numReplica);
        RingPolymer(const RingPolymer&) = default;
        RingPolymer(RingPolymer&&) noexcept = default;
        ~RingPolymer() = default;
        /* Operators */
        RingPolymer& operator=(RingPolymer obj) noexcept;
        /* Operations */
        void removeDrift();
        void toNormalRepr(size_t posID);
        void toBeadRepr(size_t posID);
        [[nodiscard]] PositionMatrix makeBeadPos(size_t replica) const;
        [[nodiscard]] PositionMatrix makeCentroidPos() const;
        [[nodiscard]] PositionMatrix makeBeadMomentum(size_t replica) const;
        [[nodiscard]] PositionMatrix makeCentroidMomentum() const;
        void swap(RingPolymer& obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] const PhaseMatrix& asMatrix() const noexcept { return phase; }
        [[nodiscard]] PhaseMatrix& asMatrix() noexcept { return phase; }
        [[nodiscard]] size_t getDOF() const noexcept { return phase.getRow() / 2U; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return getDOF() / Dim; }
        [[nodiscard]] size_t getNumReplica() const noexcept { return phase.getColumn(); }
        [[nodiscard]] FFT<PosScalarType, 1>& getCanonicalFFT() noexcept { return fft; }
        [[nodiscard]] const BufferType& getBuffer() const noexcept { return buffer; }
        [[nodiscard]] BufferType& getBuffer() noexcept { return buffer; }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return buffer.getColumn(); }
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    RingPolymer<ScalarType, PosScalarType, Dim>::RingPolymer(const MDCellType& cell, size_t numReplica)
            : phase(2 * cell.getDOF(), numReplica)
            , fft(numReplica, 1) {
        const size_t dof = getDOF();
        buffer.resize(2, fft.getKSpaceSize());

        auto momentum = phase.topRows(dof);
        momentum = ScalarType::Zero();
        /* Fill pos */ {
            size_t index = dof;
            for (auto elem : cell.getPos()) {
                phase(index, 0) = elem;
                ++index;
            }
            for (size_t i = 1; i < getNumReplica(); ++i) {
                auto col = phase.col(i);
                auto pos = col.tail(dof);
                pos = phase.col(0).tail(dof);
            }
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    RingPolymer<ScalarType, PosScalarType, Dim>&
    RingPolymer<ScalarType, PosScalarType, Dim>::operator=(RingPolymer<ScalarType, PosScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RingPolymer<ScalarType, PosScalarType, Dim>::removeDrift() {
        const size_t dof = getDOF();
        Vector<ScalarType, Dim> driftMomentum(Dim, 0);
        for (size_t i = 0; i < dof; ++i) {
            const size_t direction = i % Dim;
            for (size_t j = 0; j < getNumReplica(); ++j)
                driftMomentum[direction] += phase(i, j);
        }
        driftMomentum *= Core::reciprocal(ScalarType(getNumParticle() * getNumReplica()));

        for (size_t i = 0; i < dof; ++i) {
            auto row = phase.row(i);
            row -= driftMomentum[i % Dim];
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RingPolymer<ScalarType, PosScalarType, Dim>::toNormalRepr(size_t posID) {
        assert(posID < getDOF());
        fft.transform(phase.row(posID));
        auto momentum = buffer.row(0);
        momentum = fft.getKSpace();

        fft.transform(phase.row(posID + getDOF()));
        auto pos = buffer.row(1);
        pos = fft.getKSpace();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RingPolymer<ScalarType, PosScalarType, Dim>::toBeadRepr(size_t posID) {
        assert(posID < getDOF());
        fft.invTransform(buffer.row(0));
        auto momentum = phase.row(posID);
        momentum = fft.getRSpace();

        fft.invTransform(buffer.row(1));
        auto pos = phase.row(posID + getDOF());
        pos = fft.getRSpace();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RingPolymer<ScalarType, PosScalarType, Dim>::PositionMatrix
    RingPolymer<ScalarType, PosScalarType, Dim>::makeBeadPos(size_t replica) const {
        PositionMatrix result(getNumParticle(), Dim);
        auto col = phase.col(replica);
        size_t index = getDOF();
        for (auto& elem : result) {
            elem = PosScalarType(col[index]);
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RingPolymer<ScalarType, PosScalarType, Dim>::PositionMatrix
    RingPolymer<ScalarType, PosScalarType, Dim>::makeCentroidPos() const {
        PositionMatrix result(getNumParticle(), Dim);
        const size_t dof = getDOF();
        size_t index = dof;
        for (auto& elem : result) {
            elem = PosScalarType(mean(phase.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RingPolymer<ScalarType, PosScalarType, Dim>::PositionMatrix
    RingPolymer<ScalarType, PosScalarType, Dim>::makeBeadMomentum(size_t replica) const {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result) {
            elem = PosScalarType(phase(index, replica));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RingPolymer<ScalarType, PosScalarType, Dim>::PositionMatrix
    RingPolymer<ScalarType, PosScalarType, Dim>::makeCentroidMomentum() const {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result) {
            elem = PosScalarType(mean(phase.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RingPolymer<ScalarType, PosScalarType, Dim>::swap(RingPolymer& obj) noexcept {
        phase.swap(obj.phase);
        fft.swap(obj.fft);
        buffer.swap(obj.buffer);
    }
}
