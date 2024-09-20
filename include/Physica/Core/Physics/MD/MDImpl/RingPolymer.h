/*
 * Copyright 2023-2024 Weibo He.
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
    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    class RingPolymer {
        using PlainScalar = typename ScalarType::PlainScalar;
        using BufferScalarType = typename PlainScalar::ComplexType;
        constexpr static int PhaseMatrixMajor = NumReplica == 1 ? MatrixOption::Column : MatrixOption::Row;
    public:
        using MDCellType = MDCell<ScalarType, Dim>;
        using MassVector = typename MDCellType::MassVector;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using PhaseMatrix = DenseMatrix<ScalarType, PhaseMatrixMajor | MatrixOption::Vector, Dynamic, NumReplica>;
        using BufferType = DenseMatrix<BufferScalarType, MatrixOption::Row | MatrixOption::Vector, 2>;
        using FFTType = FFT<ScalarType, 1>;
    private:
        PhaseMatrix phase;
        MassVector massVec;
        FFTType fft;
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
        template<class KineticModel, class RandomGenerator> void initMomentum(ScalarType temperatureT, RandomGenerator& gen);
        template<class KineticModel> void scaleVelocity(ScalarType temperatureT);
        [[nodiscard]] Vector<ScalarType, Dim> makeDriftMomentum() const;
        void removeDrift();
        inline void toNormalRepr(size_t posID);
        inline void toBeadRepr(size_t posID);
        void toNormalRepr(size_t posID, const PhaseMatrix& outer_phase);
        void toBeadRepr(size_t posID, PhaseMatrix& outer_phase);
        void toNormalRepr(size_t posID, const PhaseMatrix& outer_phase, BufferType& outer_buffer, FFTType& outer_fft) const;
        void toBeadRepr(size_t posID, PhaseMatrix& outer_phase, const BufferType& outer_buffer, FFTType& outer_fft) const;
        [[nodiscard]] PositionMatrix makeBeadPos(size_t replica) const;
        [[nodiscard]] PositionMatrix makeCentroidPos() const;
        [[nodiscard]] PositionMatrix makeBeadMomentum(size_t replica) const;
        [[nodiscard]] PositionMatrix makeCentroidMomentum() const;

        [[nodiscard]] ScalarType calcKineticClassical() const;
        template<class KineticModel>
        [[nodiscard]] ScalarType calcTemperature() const;

        [[nodiscard]] ScalarType calcRepBeta(ScalarType temperatureT) const noexcept;
        [[nodiscard]] ScalarType calcOmegaW(ScalarType temperatureT) const noexcept;
        void swap(RingPolymer& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] const PhaseMatrix& asMatrix() const noexcept { return phase; }
        [[nodiscard]] PhaseMatrix& asMatrix() noexcept { return phase; }
        [[nodiscard]] size_t getDOF() const noexcept { return phase.getRow() / 2U; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return getDOF() / Dim; }
        [[nodiscard]] size_t getNumReplica() const noexcept { return phase.getColumn(); }

        [[nodiscard]] const MassVector& getMassVec() const noexcept { return massVec; }
        [[nodiscard]] FFT<ScalarType, 1>& getCanonicalFFT() noexcept { return fft; }

        [[nodiscard]] const FFTType& getFFT() const noexcept { return fft; }
        [[nodiscard]] const BufferType& getBuffer() const noexcept { return buffer; }
        [[nodiscard]] BufferType& getBuffer() noexcept { return buffer; }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return buffer.getColumn(); }
        /* Static members */
        [[nodiscard]] static ScalarType calcRepBeta(ScalarType temperatureT, size_t numReplica) noexcept;
        [[nodiscard]] static ScalarType calcOmegaW(ScalarType temperatureT, size_t numReplica) noexcept;
    };

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    RingPolymer<ScalarType, Dim, NumReplica>::RingPolymer(const MDCellType& cell, size_t numReplica)
            : phase(2 * cell.getDOF(), numReplica)
            , massVec(cell.getMassVec())
            , fft(numReplica, PlanFlag::Estimate) {
        assert(NumReplica == Dynamic || NumReplica == numReplica);
        const size_t dof = getDOF();
        buffer.resize(2, fft.getKSpaceSize());

        auto momentum = phase.topRows(dof);
        momentum = ScalarType(0);
        /* Fill pos */ {
            size_t index = dof;
            for (auto elem : cell.getPos().asArray()) {
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

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    RingPolymer<ScalarType, Dim, NumReplica>&
    RingPolymer<ScalarType, Dim, NumReplica>::operator=(RingPolymer<ScalarType, Dim, NumReplica> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    template<class KineticModel, class RandomGenerator>
    void RingPolymer<ScalarType, Dim, NumReplica>::initMomentum(ScalarType temperatureT, RandomGenerator& gen) {
        const size_t dof = getDOF();
        [[unlikely]] if (temperatureT.isZero()) {
            auto momentum = phase.topRows(dof);
            momentum = PlainScalar(0);
            return;
        }

        std::normal_distribution<> dist{};
        const ScalarType repBeta = calcRepBeta(temperatureT);
        Vector<ScalarType, Dim> driftMomentum(Dim, 0);
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = massVec[i / Dim];
            const size_t direction = i % Dim;
            const ScalarType factor = sqrt(repBeta * mass);
            for (size_t j = 0; j < getNumReplica(); ++j) {
                const ScalarType temp = factor * PlainScalar(dist(gen));
                phase(i, j) = temp;
                driftMomentum[direction] += temp;
            }
        }
        driftMomentum *= Core::reciprocal(ScalarType(getNumParticle() * getNumReplica()));

        for (size_t i = 0; i < dof; ++i) {
            auto row = phase.row(i);
            row.asVector() -= driftMomentum[i % Dim];
        }
        scaleVelocity<KineticModel>(temperatureT);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    template<class KineticModel>
    void RingPolymer<ScalarType, Dim, NumReplica>::scaleVelocity(ScalarType temperatureT) {
        const ScalarType temperatureNow = calcTemperature<KineticModel>();
        assert(!temperatureNow.isZero() && "[Error]: Velocity-scaling fails at 0K");
        const size_t dof = getDOF();
        const ScalarType factor = sqrt(temperatureT / temperatureNow);
        auto momentum = phase.topRows(dof);
        momentum.asMatrix() *= factor;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    Vector<ScalarType, Dim> RingPolymer<ScalarType, Dim, NumReplica>::makeDriftMomentum() const {
        const size_t dof = getDOF();
        Vector<ScalarType, Dim> result(Dim, 0);
        for (size_t i = 0; i < dof; ++i) {
            const size_t direction = i % Dim;
            result[direction] += phase.row(i).asVector().sum();
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void RingPolymer<ScalarType, Dim, NumReplica>::removeDrift() {
        const Vector<ScalarType, Dim> driftVelocity = makeDriftMomentum() * reciprocal(ScalarType(getNumReplica()) * massVec.sum());
        for (size_t i = 0; i < getDOF(); ++i) {
            const auto mass = massVec[i / Dim];
            auto row = phase.row(i);
            row.asVector() -= mass * driftVelocity[i % Dim];
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    inline void RingPolymer<ScalarType, Dim, NumReplica>::toNormalRepr(size_t posID) {
        toNormalRepr(posID, phase);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    inline void RingPolymer<ScalarType, Dim, NumReplica>::toBeadRepr(size_t posID) {
        toBeadRepr(posID, phase);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void RingPolymer<ScalarType, Dim, NumReplica>::toNormalRepr(size_t posID, const PhaseMatrix& outer_phase) {
        assert(posID < outer_phase.getRow() / 2);
        fft.transform(outer_phase.row(posID));
        auto momentum = buffer.row(0);
        momentum = fft.getKSpace();

        fft.transform(outer_phase.row(posID + getDOF()));
        auto pos = buffer.row(1);
        pos = fft.getKSpace();
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void RingPolymer<ScalarType, Dim, NumReplica>::toBeadRepr(size_t posID, PhaseMatrix& outer_phase) {
        assert(posID < outer_phase.getRow() / 2);
        fft.invTransform(buffer.row(0));
        auto momentum = outer_phase.row(posID);
        momentum = fft.getRSpace();

        fft.invTransform(buffer.row(1));
        auto pos = outer_phase.row(posID + getDOF());
        pos = fft.getRSpace();
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void RingPolymer<ScalarType, Dim, NumReplica>::toNormalRepr(
            size_t posID, const PhaseMatrix& outer_phase, BufferType& outer_buffer, FFTType& outer_fft) const {
        assert(posID < getDOF());
        assert(outer_fft.getRSpaceSize() == getNumReplica());
        assert(outer_fft.getKSpaceSize() == outer_buffer.getColumn());
        outer_fft.getRSpace() = outer_phase.row(posID);
        FFTType::transform(fft, outer_fft);
        auto momentum = outer_buffer.row(0);
        momentum = outer_fft.getKSpace();

        outer_fft.getRSpace() = outer_phase.row(posID + getDOF());
        FFTType::transform(fft, outer_fft);
        auto pos = outer_buffer.row(1);
        pos = outer_fft.getKSpace();
    }
    
    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void RingPolymer<ScalarType, Dim, NumReplica>::toBeadRepr(
            size_t posID, PhaseMatrix& outer_phase, const BufferType& outer_buffer, FFTType& outer_fft) const {
        assert(posID < getDOF());
        assert(outer_fft.getRSpaceSize() == getNumReplica());
        assert(outer_fft.getKSpaceSize() == outer_buffer.getColumn());
        outer_fft.getKSpace() = outer_buffer.row(0);
        FFTType::invTransform(fft, outer_fft);
        auto momentum = outer_phase.row(posID);
        momentum = outer_fft.getRSpace();

        outer_fft.getKSpace() = outer_buffer.row(1);
        FFTType::invTransform(fft, outer_fft);
        auto pos = outer_phase.row(posID + getDOF());
        pos = outer_fft.getRSpace();
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    typename RingPolymer<ScalarType, Dim, NumReplica>::PositionMatrix
    RingPolymer<ScalarType, Dim, NumReplica>::makeBeadPos(size_t replica) const {
        PositionMatrix result(getNumParticle(), Dim);
        auto col = phase.col(replica);
        size_t index = getDOF();
        for (auto& elem : result.asArray()) {
            elem = ScalarType(col[index]);
            ++index;
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    typename RingPolymer<ScalarType, Dim, NumReplica>::PositionMatrix
    RingPolymer<ScalarType, Dim, NumReplica>::makeCentroidPos() const {
        PositionMatrix result(getNumParticle(), Dim);
        const size_t dof = getDOF();
        size_t index = dof;
        for (auto& elem : result.asArray()) {
            elem = ScalarType(mean(phase.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    typename RingPolymer<ScalarType, Dim, NumReplica>::PositionMatrix
    RingPolymer<ScalarType, Dim, NumReplica>::makeBeadMomentum(size_t replica) const {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result.asArray()) {
            elem = ScalarType(phase(index, replica));
            ++index;
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    typename RingPolymer<ScalarType, Dim, NumReplica>::PositionMatrix
    RingPolymer<ScalarType, Dim, NumReplica>::makeCentroidMomentum() const {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result.asArray()) {
            elem = ScalarType(mean(phase.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    ScalarType RingPolymer<ScalarType, Dim, NumReplica>::calcKineticClassical() const {
        const size_t dof = getDOF();
        ScalarType result = 0;
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = massVec[i / Dim];
            auto p = phase.row(i);
            result += square(p.asVector()).sum() / (mass * PlainScalar(2));
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    template<class KineticModel>
    ScalarType RingPolymer<ScalarType, Dim, NumReplica>::calcTemperature() const {
        constexpr bool IsPeriodBoundary = Traits<KineticModel>::IsPeriodBoundary;
        constexpr size_t NumConstraint = IsPeriodBoundary ? 1 : 0;
        const auto factor1 = PlainScalar(2 / (Dim * PhyConst<AU>::boltzmannK));
        const auto factor2 = factor1 / PlainScalar((getNumParticle() * getNumReplica() - NumConstraint) * getNumReplica());
        return calcKineticClassical() * factor2;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    ScalarType RingPolymer<ScalarType, Dim, NumReplica>::calcRepBeta(ScalarType temperatureT) const noexcept {
        return calcRepBeta(temperatureT, getNumReplica());
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    ScalarType RingPolymer<ScalarType, Dim, NumReplica>::calcOmegaW(ScalarType temperatureT) const noexcept {
        return calcOmegaW(temperatureT, getNumReplica());
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    void RingPolymer<ScalarType, Dim, NumReplica>::swap(RingPolymer& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        phase.swap(obj.phase);
        massVec.swap(obj.massVec);
        fft.swap(obj.fft);
        buffer.swap(obj.buffer);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    ScalarType RingPolymer<ScalarType, Dim, NumReplica>::calcRepBeta(ScalarType temperatureT, size_t numReplica) noexcept {
        return temperatureT * PlainScalar(PhyConst<AU>::boltzmannK * numReplica);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica>
    ScalarType RingPolymer<ScalarType, Dim, NumReplica>::calcOmegaW(ScalarType temperatureT, size_t numReplica) noexcept {
        return calcRepBeta(temperatureT, numReplica) / PlainScalar(PhyConst<AU>::reducedPlanck);
    }
}
