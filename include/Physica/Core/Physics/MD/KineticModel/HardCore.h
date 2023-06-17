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

#include "FreeModel.h"
#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica> class RingPolymer;

    template<class ScalarType, size_t NumReplica = Dynamic, class Executor = SequentialExecutor>
    class HardCore : private FreeModel<ScalarType, ScalarType, 1, NumReplica> {
        using Base = FreeModel<ScalarType, ScalarType, 1, NumReplica>;
    public:
        using RingPolymerType = RingPolymer<ScalarType, ScalarType, 1, NumReplica>;
        using PhaseMatrix = typename RingPolymerType::PhaseMatrix;
    private:
        ScalarType latticeSize;
        ScalarType collideFactor;
        ScalarType temperatureT;
        Vector<ScalarType> repMass;
        PhaseMatrix buffer;
        size_t maxHandleNum;
    public:
        HardCore(ScalarType latticeSize_,
                 ScalarType collideFactor_,
                 ScalarType temperatureT_,
                 size_t numParticle,
                 size_t numReplica,
                 size_t maxHandleNum_);
        HardCore(const HardCore&) = default;
        HardCore(HardCore&&) noexcept = default;
        ~HardCore() = default;
        /* Operators */
        HardCore& operator=(HardCore obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT);
        void updateMass(RingPolymerType& ringPolymer);
        void swap(HardCore& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return repMass.getLength(); }
        [[nodiscard]] size_t getNumReplica() const noexcept { return buffer.getColumn(); }
        [[nodiscard]] const Vector<ScalarType>& getRepMass() const noexcept { return repMass; }
    private:
        bool checkCollision([[maybe_unused]] size_t id_dof, const RingPolymerType& ringPolymer) const;
        bool handleCollision(const RingPolymerType& ringPolymer);
        bool checkRepMass() const;
    };

    template<class ScalarType, size_t NumReplica, class Executor>
    HardCore<ScalarType, NumReplica, Executor>::HardCore(
            ScalarType latticeSize_, ScalarType collideFactor_, ScalarType temperatureT_, size_t numParticle, size_t numReplica, size_t maxHandleNum_)
            : Base(temperatureT_, numReplica)
            , latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , repMass(numParticle, 0)
            , buffer(numParticle * 2, numReplica)
            , maxHandleNum(maxHandleNum_) {
        assert(collideFactor < ScalarType(1.0) && collideFactor.isPositive());
        assert(NumReplica == Dynamic || NumReplica == numReplica);
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    HardCore<ScalarType, NumReplica, Executor>&
    HardCore<ScalarType, NumReplica, Executor>::operator=(HardCore<ScalarType, NumReplica, Executor> obj) noexcept {
        swap(*this);
        return *this;
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    void HardCore<ScalarType, NumReplica, Executor>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT) {
        const size_t numParticle = getNumParticle();
        const ScalarType collideStep = collideFactor * deltaT;
        auto& phase = ringPolymer.asMatrix();
        assert(numParticle == ringPolymer.getNumParticle());
        assert(getNumReplica() == ringPolymer.getNumReplica());
        assert(checkRepMass());

        buffer = phase;
        ScalarType lStep = 0;
        ScalarType rStep = deltaT;
        ScalarType from = 0;
        ScalarType to = deltaT;
        size_t handleNum = 0;
        bool isDrifted = false;
        while (true) {
            const ScalarType step = to - from;
            if constexpr (NumReplica != 1) {
                Base::pre_nve_step_impl(ringPolymer, step);
                for (size_t i = 0; i < numParticle; ++i) {
                    Base::do_nve_step_impl(i, ringPolymer, buffer, phase, step);
                    [[unlikely]] if (checkCollision(i, ringPolymer)) {
                        rStep = to;
                        goto done;
                    }
                }
                lStep = to;
            done:;
            }
            else {
                auto col = phase.col(0);
                auto pos = col.tail(numParticle);
                auto col1 = buffer.col(0);
                auto momentum_buffer = col1.head(numParticle);
                auto pos_buffer = col1.tail(numParticle);
                pos = pos_buffer + hadamard(momentum_buffer, repMass) * step;
                constexpr int unused = 0;
                if (checkCollision(unused, ringPolymer))
                    rStep = to;
                else
                    lStep = to;
            }
            to = (lStep + rStep) * 0.5;

            const bool isDone = lStep == deltaT;
            if (isDone) {
                if constexpr (NumReplica == 1) {
                    auto momentum = phase.topRows(numParticle);
                    momentum = buffer.topRows(numParticle);
                }
                break;
            }

            const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
            if (isDeltaSmallEnough) {
                if (handleNum == maxHandleNum) [[unlikely]]
                    throw BadConvergenceException("[Error]: Too many collision with in a step");

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
                handleNum += 1;
                from = lStep;
                to = deltaT;
                isDrifted |= handleCollision(ringPolymer);
            }
        }

        if (isDrifted) {
            const ScalarType temperatureT = ringPolymer.calcTemperature();
            ringPolymer.removeDrift();
            ringPolymer.scaleVelocity(temperatureT);
        }
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    void HardCore<ScalarType, NumReplica, Executor>::updateMass(RingPolymerType& ringPolymer) {
        repMass = reciprocal(ringPolymer.getMassVec());
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    void HardCore<ScalarType, NumReplica, Executor>::swap(HardCore& obj) noexcept {
        Base::swap(obj);
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        std::swap(maxHandleNum, obj.maxHandleNum);
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    bool HardCore<ScalarType, NumReplica, Executor>::checkCollision(
            [[maybe_unused]] size_t id_dof, const RingPolymerType& ringPolymer) const {
        using PacketType = typename Internal::BestPacket<ScalarType, Dynamic>::Type;
        const size_t numParticle = getNumParticle();
        if constexpr (NumReplica == 1) {
            const size_t length = numParticle - 1;
            const size_t to = length / PacketType::size() * PacketType::size();
            auto phase = ringPolymer.asMatrix().col(0);
            auto pos = phase.tail(numParticle);
            if (pos[0].isNegative()) [[unlikely]]
                return true;

            auto head = pos.head(length);
            auto tail = pos.tail(1);
            {
                size_t i = 0;
                for (; i < to; i += PacketType::size()) {
                    const auto boolPacket = head.template packet<PacketType>(i) > tail.template packet<PacketType>(i);
                    if (horizontal_or(boolPacket)) [[unlikely]]
                        return true;
                }

                if (to != length) {
                    const size_t count = length - i;
                    const auto boolPacket = head.template packetPartial<PacketType>(i, count) > tail.template packetPartial<PacketType>(i, count);
                    if (horizontal_or(boolPacket)) [[unlikely]]
                        return true;
                }
            }
            if (pos[length] > latticeSize) [[unlikely]]
                return true;
        }
        else {
            const size_t numReplica = getNumReplica();
            const size_t to = numReplica / PacketType::size() * PacketType::size();
            auto pos = ringPolymer.asMatrix().row(numParticle + id_dof);
            [[unlikely]] if (id_dof == 0) {
                const PacketType zeros(0);
                size_t i = 0;
                for (; i < to; i += PacketType::size()) {
                    const auto boolPacket = pos.template packet<PacketType>(i) < zeros;
                    if (horizontal_or(boolPacket)) [[unlikely]]
                        return true;
                }

                if (to != numReplica) {
                    const size_t count = numReplica - i;
                    const auto boolPacket = pos.template packetPartial<PacketType>(i, count) < zeros;
                    if (horizontal_or(boolPacket)) [[unlikely]]
                        return true;
                }
            }
            else {
                /* Particles far from the wall */ {
                    auto pos0 = ringPolymer.asMatrix().row(numParticle + id_dof - 1);
                    size_t i = 0;
                    for (; i < to; i += PacketType::size()) {
                        const auto boolPacket = pos0.template packet<PacketType>(i) > pos.template packet<PacketType>(i);
                        if (horizontal_or(boolPacket)) [[unlikely]]
                            return true;
                    }

                    if (to != numReplica) {
                        const size_t count = numReplica - i;
                        const auto boolPacket = pos0.template packetPartial<PacketType>(i, count) > pos.template packetPartial<PacketType>(i, count);
                        if (horizontal_or(boolPacket)) [[unlikely]]
                            return true;
                    }
                }
                [[unlikely]] if (id_dof == numParticle - 1) {
                    const PacketType latticeSizes(latticeSize.getTrivial());
                    size_t i = 0;
                    for (; i < to; i += PacketType::size()) {
                        const auto boolPacket = pos.template packet<PacketType>(i) > latticeSizes;
                        if (horizontal_or(boolPacket)) [[unlikely]]
                            return true;
                    }

                    if (to != numReplica) {
                        const size_t count = numReplica - i;
                        const auto boolPacket = pos.template packetPartial<PacketType>(i, count) > latticeSizes;
                        if (horizontal_or(boolPacket)) [[unlikely]]
                            return true;
                    }
                }
            }
        }
        return false;
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    bool HardCore<ScalarType, NumReplica, Executor>::handleCollision(const RingPolymerType& ringPolymer) {
        const size_t numReplica = getNumReplica();
        const size_t numParticle = getNumParticle();
        const auto& mass = ringPolymer.getMassVec();
        auto momentumMatrix = buffer.topRows(numParticle);
        auto posMatrix = ringPolymer.asMatrix().bottomRows(numParticle);
        bool isDrifted = false;
        for (size_t replica = 0; replica < numReplica; ++replica) {
            auto momentum = momentumMatrix.col(replica);
            auto pos = posMatrix.col(replica);
            size_t i = 0;
            for (; i < numParticle - 1; ++i) {
                if (pos[i] > pos[i + 1]) {
                    const ScalarType m1 = mass[i];
                    const ScalarType m2 = mass[i + 1];
                    const ScalarType p1 = momentum[i];
                    const ScalarType p2 = momentum[i + 1];
                    const ScalarType next_p1 = ((m1 - m2) * p1 + ScalarType(2) * m1 * p2) * reciprocal(m1 + m2);
                    const ScalarType next_p2 = p1 + p2 - next_p1;
                    momentum[i] = next_p1;
                    momentum[i + 1] = next_p2;
                }
            }
            if (!pos[0].isPositive()) {
                momentum[0] = abs(momentum[0]);
                isDrifted = true;
            }
            if (pos[i] > latticeSize) {
                momentum[i] = -abs(momentum[i]);
                isDrifted = true;
            }
        }
        return isDrifted;
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    bool HardCore<ScalarType, NumReplica, Executor>::checkRepMass() const {
        bool isGood = true;
        for (ScalarType elem : repMass)
            isGood &= elem.isPositive();
        return isGood;
    }
}
