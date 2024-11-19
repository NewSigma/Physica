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

#include <Physica/Core/Exception/BadConvergenceException.h>
#include "FreeModel.h"

namespace Physica::Core {
    template<Scalar T, unsigned int Dim, size_t NumReplica> class RingPolymer;

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor = SequentialExecutor>
    class HardCore : private FreeModel<T, 1, NumReplica, Integrator> {
        using Base = FreeModel<T, 1, NumReplica, Integrator>;
        using ValueType = typename T::ValueType;
    public:
        using RingPolymerType = RingPolymer<T, 1, NumReplica>;
        using PhaseMatrix = typename RingPolymerType::PhaseMatrix;
    private:
        T latticeSize;
        ValueType collideFactor;
        T temperatureT;
        VectorND<T> repMass;
        PhaseMatrix buffer;
        size_t maxHandleNum;
        size_t handleNum;
    public:
        HardCore(T latticeSize_,
                 ValueType collideFactor_,
                 T temperatureT_,
                 size_t numParticle,
                 size_t numReplica,
                 size_t maxHandleNum_);
        HardCore(const HardCore&) = default;
        HardCore(HardCore&&) noexcept = default;
        ~HardCore() = default;
        /* Operators */
        HardCore& operator=(HardCore obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, T deltaT);
        void updateMass(const RingPolymerType& ringPolymer);
        void swap(HardCore& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return repMass.getLength(); }
        [[nodiscard]] size_t getNumReplica() const noexcept { return buffer.getCol(); }
        [[nodiscard]] const VectorND<T>& getRepMass() const noexcept { return repMass; }
        [[nodiscard]] size_t getHandleNum() { return handleNum; }
        /* Static members */
        static void checkParam(ValueType collideFactor, size_t numReplica);
        [[nodiscard]] __host__ __device__ static bool checkStepSize(
            T latticeSize,
            T temperatureT,
            T collideStep,
            T maxMass);
    private:
        bool checkCollision([[maybe_unused]] size_t id_dof, const RingPolymerType& ringPolymer) const;
        void handleCollision(const RingPolymerType& ringPolymer);
        bool checkRepMass() const;
    };

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::HardCore(
            T latticeSize_, ValueType collideFactor_, T temperatureT_, size_t numParticle, size_t numReplica, size_t maxHandleNum_)
            : Base(temperatureT_, numReplica)
            , latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , temperatureT(temperatureT_)
            , repMass(numParticle, 0)
            , buffer(numParticle * 2, numReplica)
            , maxHandleNum(maxHandleNum_) {
        checkParam(collideFactor, numReplica);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::nve_step(RingPolymerType& ringPolymer, T deltaT) {
        const size_t numParticle = getNumParticle();
        const ValueType collideStep = collideFactor * deltaT.getValue();
        auto& phase = ringPolymer.asMatrix();
        assert(numParticle == ringPolymer.getNumParticle());
        assert(getNumReplica() == ringPolymer.getNumReplica());
        assert(checkRepMass() && "[Error]: Mass is not initialized");
        assert(checkStepSize(latticeSize, temperatureT, collideStep, ringPolymer.getMassVec().max())
            && "[Warn]: Stepsize is too small, add precision or adjust params is recommanded");

        buffer = phase;
        T lStep = 0;
        T rStep = deltaT;
        T from = 0;
        T to = deltaT;
        handleNum = 0;
        while (true) {
            const T step = to - from;
            bool isCollided = false;
            if constexpr (NumReplica != 1) {
                Base::pre_nve_step_impl(ringPolymer, step);
                if constexpr (IsFixedBoundary) {
                    for (size_t i = 0; i < numParticle; ++i) {
                        Base::do_nve_step_impl(i, ringPolymer, buffer, phase);
                        isCollided = checkCollision(i, ringPolymer);
                        if (isCollided) [[unlikely]]
                            break;
                    }
                }
                else {
                    Base::do_nve_step_impl(0, ringPolymer, buffer, phase);
                    for (size_t i = 1; i < numParticle; ++i) {
                        Base::do_nve_step_impl(i, ringPolymer, buffer, phase);
                        isCollided = checkCollision(i, ringPolymer);
                        if (isCollided) [[unlikely]]
                            break;
                    }
                }
            }
            else {
                auto col = phase.col(0);
                auto pos = col.tail(numParticle);
                auto col1 = buffer.col(0);
                auto momentum_buffer = col1.head(numParticle);
                auto pos_buffer = col1.tail(numParticle);
                pos = pos_buffer + hadamard(momentum_buffer, repMass) * step;
                constexpr int unused = 0;
                isCollided = checkCollision(unused, ringPolymer);
            }

            if (isCollided)
                rStep = to;
            else {
                lStep = to;
                const bool isDone = lStep == deltaT;
                if (isDone)
                    break;
            }

            const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
            if (isDeltaSmallEnough) {
                if (handleNum == maxHandleNum) [[unlikely]]
                    throw BadConvergenceException("[Error]: Too many collision within a step, possibly numerical error is large");

                if constexpr (NumReplica != 1) {
                    Base::nve_step_impl(ringPolymer, buffer, phase, rStep - from);
                    Base::nve_step_impl(ringPolymer, buffer, buffer, lStep - from);
                }
                else {
                    auto col = phase.col(0);
                    auto pos = col.tail(numParticle);
                    auto col1 = buffer.col(0);
                    auto momentum_buffer = col1.head(numParticle);
                    auto pos_buffer = col1.tail(numParticle);
                    pos = pos_buffer + hadamard(momentum_buffer, repMass) * (rStep - from);
                    pos_buffer += hadamard(momentum_buffer, repMass) * (lStep - from);
                }
                from = lStep;
                to = deltaT;
                handleCollision(ringPolymer);
                handleNum += 1;
            }
            else
                to = (lStep + rStep) * ValueType(0.5);
        }

        if constexpr (NumReplica == 1) {
            auto momentum = phase.topRows(numParticle);
            momentum = buffer.topRows(numParticle);
        }
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::updateMass(const RingPolymerType& ringPolymer) {
        repMass = reciprocal(ringPolymer.getMassVec());
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::swap(HardCore& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        temperatureT.swap(obj.temperatureT);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        std::swap(maxHandleNum, obj.maxHandleNum);
        std::swap(handleNum, obj.handleNum);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::checkParam(ValueType collideFactor, size_t numReplica) {
        if (!(collideFactor < ValueType(1.0) && collideFactor.isPositive())) [[unlikely]]
            throw std::invalid_argument("[Error]: Collide factor must be in (0, 1)");
        if (NumReplica != Dynamic && NumReplica != numReplica) [[unlikely]]
            throw std::invalid_argument("[Error]: Number of replica is not consistent");
        if (collideFactor <= ValueType(std::numeric_limits<T>::epsilon())) [[unlikely]]
            throw std::invalid_argument("[Error]: Collide factor is too small, numerical error will be large");
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    bool HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::checkStepSize(
            T latticeSize,
            T temperatureT,
            T collideStep,
            T maxMass) {
        const ValueType epsilonStep = latticeSize * std::numeric_limits<ValueType>::epsilon();
        const ValueType meanVelocity = sqrt(ValueType(PhyConst<AU>::boltzmannK) * temperatureT / maxMass);
        return collideStep * meanVelocity > epsilonStep;
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    bool HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::checkCollision(
            [[maybe_unused]] size_t id_dof, const RingPolymerType& ringPolymer) const {
        using PacketType = typename BestPacket<T, Dynamic>::Type;
        using BoolPacketType = typename Traits<PacketType>::BoolSIMDType;
        const size_t numParticle = getNumParticle();
        if constexpr (NumReplica == 1) {
            auto phase = ringPolymer.asMatrix().col(0);
            auto pos = phase.tail(numParticle);
            const size_t length = numParticle - 1;
            {
                auto head = pos.head(length);
                auto tail = pos.tail(1);
                const size_t to = length / PacketType::size() * PacketType::size();
                size_t i = 0;
                for (; i < to; i += PacketType::size()) {
                    const auto boolPacket = BoolPacketType(head.template packet<PacketType>(i) > tail.template packet<PacketType>(i));
                    if (boolPacket.horizontal_or()) [[unlikely]]
                        return true;
                }

                if (to != length) {
                    for (size_t i = to; i < length; ++i)
                        if (head[i] > tail[i]) [[unlikely]]
                            return true;
                }
            }
            if constexpr (IsFixedBoundary) {
                if (pos[0].isNegative() || pos[length] > latticeSize) [[unlikely]]
                    return true;
            }
            else {
                if (pos[length] - latticeSize > pos[0]) [[unlikely]]
                    return true;
            }
        }
        else {
            const size_t numReplica = getNumReplica();
            const size_t to = numReplica / PacketType::size() * PacketType::size();
            auto pos = ringPolymer.asMatrix().row(numParticle + id_dof);
            [[unlikely]] if (id_dof == 0) {
                if constexpr (IsFixedBoundary) {
                    const PacketType zeros(0);
                    size_t i = 0;
                    for (; i < to; i += PacketType::size()) {
                        const auto boolPacket = BoolPacketType(pos.template packet<PacketType>(i) < zeros);
                        if (boolPacket.horizontal_or()) [[unlikely]]
                            return true;
                    }

                    if (to != numReplica) {
                        for (size_t i = to; i < numReplica; ++i)
                            if (pos[i].isNegative()) [[unlikely]]
                                return true;
                    }
                }
                else {
                    auto pos_end = ringPolymer.asMatrix().row(numParticle * 2 - 1);
                    const PacketType latticeSizes(latticeSize.toMachine());
                    size_t i = 0;
                    for (; i < to; i += PacketType::size()) {
                        const PacketType pack1 = pos_end.template packet<PacketType>(i) - latticeSizes;
                        const PacketType pack2 = pos.template packet<PacketType>(i);
                        const auto boolPacket = BoolPacketType(pack1 > pack2);
                        if (boolPacket.horizontal_or()) [[unlikely]]
                            return true;
                    }

                    if (to != numReplica) {
                        for (size_t i = to; i < numReplica; ++i)
                            if (pos_end[i] > pos[i] + latticeSize) [[unlikely]]
                                return true;
                    }
                }
            }
            else {
                /* Particles far from the wall */ {
                    auto pos0 = ringPolymer.asMatrix().row(numParticle + id_dof - 1);
                    size_t i = 0;
                    for (; i < to; i += PacketType::size()) {
                        const auto boolPacket = BoolPacketType(pos0.template packet<PacketType>(i) > pos.template packet<PacketType>(i));
                        if (boolPacket.horizontal_or()) [[unlikely]]
                            return true;
                    }

                    if (to != numReplica) {
                        for (size_t i = to; i < numReplica; ++i)
                            if (pos0[i] > pos[i]) [[unlikely]]
                                return true;
                    }
                }
                if constexpr (IsFixedBoundary) {
                    [[unlikely]] if (id_dof == numParticle - 1) {
                        const PacketType latticeSizes(latticeSize.getValue().toMachine());
                        size_t i = 0;
                        for (; i < to; i += PacketType::size()) {
                            const auto boolPacket = BoolPacketType(pos.template packet<PacketType>(i) > latticeSizes);
                            if (boolPacket.horizontal_or()) [[unlikely]]
                                return true;
                        }

                        if (to != numReplica) {
                        for (size_t i = to; i < numReplica; ++i)
                            if (pos[i] > latticeSize) [[unlikely]]
                                return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::handleCollision(const RingPolymerType& ringPolymer) {
        const size_t numReplica = getNumReplica();
        const size_t numParticle = getNumParticle();
        const auto& mass = ringPolymer.getMassVec();
        auto momentumMatrix = buffer.topRows(numParticle);
        auto posMatrix = ringPolymer.asMatrix().bottomRows(numParticle);
        for (size_t replica = 0; replica < numReplica; ++replica) {
            auto momentum = momentumMatrix.col(replica);
            auto pos = posMatrix.col(replica);
            size_t i = 0;
            for (; i < numParticle - 1; ++i) {
                if (pos[i] > pos[i + 1]) {
                    const T m1 = mass[i];
                    const T m2 = mass[i + 1];
                    const T p1 = momentum[i];
                    const T p2 = momentum[i + 1];
                    const T next_p1 = ((m1 - m2) * p1 + T(2) * m1 * p2) * reciprocal(m1 + m2);
                    const T next_p2 = p1 + p2 - next_p1;
                    momentum[i] = next_p1;
                    momentum[i + 1] = next_p2;
                }
            }
            if constexpr (IsFixedBoundary) {
                if (!pos[0].isPositive())
                    momentum[0] = abs(momentum[0]);

                if (pos[i] > latticeSize)
                    momentum[i] = -abs(momentum[i]);
            }
            else {
                if (pos[i] - latticeSize > pos[0]) {
                    const T m1 = mass[i];
                    const T m2 = mass[0];
                    const T p1 = momentum[i];
                    const T p2 = momentum[0];
                    const T next_p1 = ((m1 - m2) * p1 + T(2) * m1 * p2) * reciprocal(m1 + m2);
                    const T next_p2 = p1 + p2 - next_p1;
                    momentum[i] = next_p1;
                    momentum[0] = next_p2;
                }
            }
        }
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    bool HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>::checkRepMass() const {
        bool isGood = true;
        for (T elem : repMass)
            isGood &= elem.isPositive();
        return isGood;
    }
}

namespace Physica {
    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator, class Executor>
    class Traits<Core::HardCore<T, IsFixedBoundary, NumReplica, Integrator, Executor>> {
        static_assert(!T::isComplex);
    public:
        constexpr static bool IsPeriodBoundary = !IsFixedBoundary;
    };
}
