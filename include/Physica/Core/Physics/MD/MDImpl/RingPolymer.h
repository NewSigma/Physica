/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica {
    template<Scalar T, unsigned int Dim, size_t NumReplica>
    class RingPolymer {
        using This = RingPolymer<T, Dim, NumReplica>;
        using Tv = T::ValueType;
        using Tcv = Tv::ComplexType;
        constexpr static int PhaseMatrixMajor = NumReplica == 1 ? MatrixOption::Col : MatrixOption::Row;
    public:
        using MDCellType = MDCell<T, Dim>;
        using MassVector = MDCellType::MassVector;
        using PositionMatrix = MDCellType::PositionMatrix;
        using PhaseMatrix = DenseMatrix<T, PhaseMatrixMajor, Dynamic, NumReplica>;
        using BufferType = DenseMatrix<Tcv, MatrixOption::Row, 2>;
        using FFTType = FFT<T, 1>;
    private:
        PhaseMatrix phase;
        MassVector massVec;
        FFTType fft;
        BufferType buffer;
    public:
        RingPolymer() = default;
        RingPolymer(const MDCellType& cell, size_t numReplica);
        RingPolymer(const This&) = default;
        RingPolymer(This&&) noexcept = default;
        ~RingPolymer() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class KineticModel, RNG R> void initMomentum(T temperatureT);
        template<class KineticModel> void scaleVelocity(T temperatureT);
        [[nodiscard]] DenseVector<T, Dim> makeDriftMomentum() const;
        void removeDrift();
        void toNormalRepr(size_t posID);
        void toBeadRepr(size_t posID);
        void toNormalRepr(size_t posID, const PhaseMatrix& outer_phase);
        void toBeadRepr(size_t posID, PhaseMatrix& outer_phase);
        void toNormalRepr(size_t posID, const PhaseMatrix& outer_phase, BufferType& outer_buffer, FFTType& outer_fft) const;
        void toBeadRepr(size_t posID, PhaseMatrix& outer_phase, const BufferType& outer_buffer, FFTType& outer_fft) const;
        [[nodiscard]] PositionMatrix makeBeadPos(size_t replica) const;
        [[nodiscard]] PositionMatrix makeCentroidPos() const;
        [[nodiscard]] PositionMatrix makeBeadMomentum(size_t replica) const;
        [[nodiscard]] PositionMatrix makeCentroidMomentum() const;

        [[nodiscard]] T calcKineticClassical() const;
        template<class KineticModel>
        [[nodiscard]] T calcTemperature() const;

        [[nodiscard]] T calcRepBeta(T temperatureT) const noexcept;
        [[nodiscard]] T calcOmegaW(T temperatureT) const noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] auto&& asMatrix(this auto&&) noexcept;
        [[nodiscard]] size_t getDOF() const noexcept { return phase.getRow() / 2U; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return getDOF() / Dim; }
        [[nodiscard]] size_t getNumReplica() const noexcept { return phase.getCol(); }

        [[nodiscard]] const MassVector& getMassVec() const noexcept { return massVec; }
        [[nodiscard]] FFT<T, 1>& getCanonicalFFT() noexcept { return fft; }

        [[nodiscard]] const FFTType& getFFT() const noexcept { return fft; }
        [[nodiscard]] auto&& getBuffer(this auto&&) noexcept;
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return buffer.getCol(); }
        /* Static members */
        [[nodiscard]] static T calcRepBeta(T temperatureT, size_t numReplica) noexcept;
        [[nodiscard]] static T calcOmegaW(T temperatureT, size_t numReplica) noexcept;
    };

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    RingPolymer<T, Dim, NumReplica>::RingPolymer(const MDCellType& cell, size_t numReplica)
            : phase(2 * cell.getDOF(), numReplica)
            , massVec(cell.getMassVec())
            , fft(numReplica, PlanFlag::Estimate) {
        assert(NumReplica == Dynamic || NumReplica == numReplica);
        const size_t dof = getDOF();
        buffer.resize(2, fft.getKSpaceSize());

        auto momentum = phase.topRows(dof);
        momentum = T(0);
        /* Fill pos */ {
            size_t index = dof;
            for (auto elem : cell.getPos().asArray()) {
                phase[index, 0] = elem;
                ++index;
            }
            for (size_t i = 1; i < getNumReplica(); ++i) {
                auto col = phase.col(i);
                auto pos = col.tail(dof);
                pos = phase.col(0).tail(dof);
            }
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    template<class KineticModel, RNG R>
    void RingPolymer<T, Dim, NumReplica>::initMomentum(T temperatureT) {
        const size_t dof = getDOF();
        if (temperatureT.isZero()) [[unlikely]] {
            auto momentum = phase.topRows(dof);
            momentum = Tv(0);
            return;
        }

        constexpr bool IsPeriodBoundary = Traits<KineticModel>::IsPeriodBoundary;
        const T repBeta = calcRepBeta(temperatureT);
        phase.topRows(dof).template random_normal<R>();
        if constexpr (IsPeriodBoundary) {
            DenseVector<T, Dim> driftMomentum(Dim, 0);
            for (size_t i = 0; i < getNumParticle(); ++i) {
                phase.rows(i * Dim, Dim) *= sqrt(repBeta * massVec[i]);
                for (int dim = 0; dim < int(Dim); ++dim)
                    driftMomentum[dim] += phase.row(i * Dim + dim).sum();
            }
            driftMomentum *= reciprocal(T(getNumParticle() * getNumReplica()));

            for (size_t i = 0; i < dof; ++i)
                phase.row(i) -= driftMomentum[i % Dim];
        }
        else {
            for (size_t i = 0; i < getNumParticle(); ++i)
                phase.rows(i * Dim, Dim) *= sqrt(repBeta * massVec[i]);
        }
        scaleVelocity<KineticModel>(temperatureT);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    template<class KineticModel>
    void RingPolymer<T, Dim, NumReplica>::scaleVelocity(T temperatureT) {
        const T temperatureT0 = calcTemperature<KineticModel>();
        assert(temperatureT0.isPositive() && "[Error]: Invalid temperature");
        phase.topRows(getDOF()) *= sqrt(temperatureT / temperatureT0);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    DenseVector<T, Dim> RingPolymer<T, Dim, NumReplica>::makeDriftMomentum() const {
        const size_t dof = getDOF();
        DenseVector<T, Dim> result(Dim, 0);
        for (size_t i = 0; i < dof; ++i) {
            const size_t direction = i % Dim;
            result[direction] += phase.row(i).sum();
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::removeDrift() {
        const DenseVector<T, Dim> driftVelocity = makeDriftMomentum() * reciprocal(T(getNumReplica()) * massVec.sum());
        for (size_t i = 0; i < getDOF(); ++i) {
            const auto mass = massVec[i / Dim];
            auto row = phase.row(i);
            row -= mass * driftVelocity[i % Dim];
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::toNormalRepr(size_t posID) {
        toNormalRepr(posID, phase);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::toBeadRepr(size_t posID) {
        toBeadRepr(posID, phase);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::toNormalRepr(size_t posID, const PhaseMatrix& outer_phase) {
        assert(posID < outer_phase.getRow() / 2);
        fft.transform(outer_phase.row(posID));
        auto momentum = buffer.row(0);
        momentum = fft.getKSpace();

        fft.transform(outer_phase.row(posID + getDOF()));
        auto pos = buffer.row(1);
        pos = fft.getKSpace();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::toBeadRepr(size_t posID, PhaseMatrix& outer_phase) {
        assert(posID < outer_phase.getRow() / 2);
        fft.invTransform(buffer.row(0));
        auto momentum = outer_phase.row(posID);
        momentum = fft.getRSpace();

        fft.invTransform(buffer.row(1));
        auto pos = outer_phase.row(posID + getDOF());
        pos = fft.getRSpace();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::toNormalRepr(
            size_t posID, const PhaseMatrix& outer_phase, BufferType& outer_buffer, FFTType& outer_fft) const {
        assert(posID < getDOF());
        assert(outer_fft.getRSpaceSize() == getNumReplica());
        assert(outer_fft.getKSpaceSize() == outer_buffer.getCol());
        outer_fft.getRSpace() = outer_phase.row(posID);
        FFTType::transform(fft, outer_fft);
        auto momentum = outer_buffer.row(0);
        momentum = outer_fft.getKSpace();

        outer_fft.getRSpace() = outer_phase.row(posID + getDOF());
        FFTType::transform(fft, outer_fft);
        auto pos = outer_buffer.row(1);
        pos = outer_fft.getKSpace();
    }
    
    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::toBeadRepr(
            size_t posID, PhaseMatrix& outer_phase, const BufferType& outer_buffer, FFTType& outer_fft) const {
        assert(posID < getDOF());
        assert(outer_fft.getRSpaceSize() == getNumReplica());
        assert(outer_fft.getKSpaceSize() == outer_buffer.getCol());
        outer_fft.getKSpace() = outer_buffer.row(0);
        FFTType::invTransform(fft, outer_fft);
        auto momentum = outer_phase.row(posID);
        momentum = outer_fft.getRSpace();

        outer_fft.getKSpace() = outer_buffer.row(1);
        FFTType::invTransform(fft, outer_fft);
        auto pos = outer_phase.row(posID + getDOF());
        pos = outer_fft.getRSpace();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    auto RingPolymer<T, Dim, NumReplica>::makeBeadPos(size_t replica) const -> PositionMatrix {
        PositionMatrix result(getNumParticle(), Dim);
        auto col = phase.col(replica);
        size_t index = getDOF();
        for (auto& elem : result.asArray()) {
            elem = T(col[index]);
            ++index;
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    auto RingPolymer<T, Dim, NumReplica>::makeCentroidPos() const -> PositionMatrix {
        PositionMatrix result(getNumParticle(), Dim);
        const size_t dof = getDOF();
        size_t index = dof;
        for (auto& elem : result.asArray()) {
            elem = T(phase.row(index).mean());
            ++index;
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    auto RingPolymer<T, Dim, NumReplica>::makeBeadMomentum(size_t replica) const -> PositionMatrix {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result.asArray()) {
            elem = T(phase(index, replica));
            ++index;
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    auto RingPolymer<T, Dim, NumReplica>::makeCentroidMomentum() const -> PositionMatrix {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result.asArray()) {
            elem = T(phase.row(index).mean());
            ++index;
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    T RingPolymer<T, Dim, NumReplica>::calcKineticClassical() const {
        const size_t dof = getDOF();
        T result = 0;
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = massVec[i / Dim];
            auto p = phase.row(i);
            result += square(p).sum() / (mass * Tv(2));
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    template<class KineticModel>
    T RingPolymer<T, Dim, NumReplica>::calcTemperature() const {
        constexpr bool IsPeriodBoundary = Traits<KineticModel>::IsPeriodBoundary;
        constexpr size_t NumConstraint = IsPeriodBoundary ? 1 : 0;
        assert(getNumParticle() * getNumReplica() > NumConstraint && "[Error]: Temperature of a periodic boundary system with only one atom is not well defined");
        const auto factor1 = Tv(2 / (Dim * PhyConst<AU>::boltzmannK));
        const auto factor2 = factor1 / Tv((getNumParticle() * getNumReplica() - NumConstraint) * getNumReplica());
        return calcKineticClassical() * factor2;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    T RingPolymer<T, Dim, NumReplica>::calcRepBeta(T temperatureT) const noexcept {
        return calcRepBeta(temperatureT, getNumReplica());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    T RingPolymer<T, Dim, NumReplica>::calcOmegaW(T temperatureT) const noexcept {
        return calcOmegaW(temperatureT, getNumReplica());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    void RingPolymer<T, Dim, NumReplica>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        phase.swap(obj.phase);
        massVec.swap(obj.massVec);
        fft.swap(obj.fft);
        buffer.swap(obj.buffer);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    auto&& RingPolymer<T, Dim, NumReplica>::asMatrix(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).phase;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    auto&& RingPolymer<T, Dim, NumReplica>::getBuffer(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).buffer;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    T RingPolymer<T, Dim, NumReplica>::calcRepBeta(T temperatureT, size_t numReplica) noexcept {
        return temperatureT * Tv(PhyConst<AU>::boltzmannK * numReplica);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica>
    T RingPolymer<T, Dim, NumReplica>::calcOmegaW(T temperatureT, size_t numReplica) noexcept {
        return calcRepBeta(temperatureT, numReplica) / Tv(PhyConst<AU>::reducedPlanck);
    }
}
