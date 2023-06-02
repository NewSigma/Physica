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

#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica> class RingPolymer;

    template<class ScalarType, size_t NumReplica = Dynamic, class Executor = SequentialExecutor>
    class HardCore {
    public:
        using RingPolymerType = RingPolymer<ScalarType, ScalarType, 1, NumReplica>;
        using PhaseMatrix = typename RingPolymerType::PhaseMatrix;
    private:
        ScalarType latticeSize;
        ScalarType collideFactor;
        ScalarType temperatureT;
        Vector<ScalarType> repMass;
        PhaseMatrix buffer;
        PhaseMatrix velocity;
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
        [[nodiscard]] const PhaseMatrix& getVelocity() const noexcept { return velocity; }
    private:
        void updateVelocity(const RingPolymerType& ringPolymer);
        bool checkCollision(const RingPolymerType& ringPolymer) const;
        bool handleCollision(RingPolymerType& ringPolymer);
        bool checkRepMass() const;
    };

    template<class ScalarType, size_t NumReplica, class Executor>
    HardCore<ScalarType, NumReplica, Executor>::HardCore(
            ScalarType latticeSize_, ScalarType collideFactor_, ScalarType temperatureT_, size_t numParticle, size_t numReplica, size_t maxHandleNum_)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , temperatureT(temperatureT_)
            , repMass(numParticle, 0)
            , buffer(numParticle, numReplica)
            , velocity(numParticle, numReplica)
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
        auto pos = phase.bottomRows(numParticle);
        assert(numParticle == ringPolymer.getNumParticle());
        assert(numReplica == ringPolymer.getNumReplica());
        assert(checkRepMass());

        buffer = pos;
        updateVelocity(ringPolymer);
        ScalarType lStep = 0;
        ScalarType rStep = deltaT;
        ScalarType from = 0;
        ScalarType to = deltaT;
        size_t handleNum = 0;
        bool isDrifted = false;
        while (true) {
            const ScalarType step = to - from;
            pos = buffer + velocity * step;
            if (checkCollision(ringPolymer))
                rStep = to;
            else
                lStep = to;
            to = (lStep + rStep) * 0.5;

            const bool isDone = lStep == deltaT;
            if (isDone)
                break;

            const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
            if (isDeltaSmallEnough) {
                if (handleNum == maxHandleNum) [[unlikely]]
                    throw BadConvergenceException("[Error]: Too many collision with in a step");
                handleNum += 1;
                pos = buffer + velocity * (rStep - from);
                buffer += velocity * (lStep - from);
                from = lStep;
                to = deltaT;
                rStep = deltaT;
                isDrifted |= handleCollision(ringPolymer);
                updateVelocity(ringPolymer);
            }
        }

        if (isDrifted) {
            ringPolymer.removeDrift();
            const size_t numParticle = ringPolymer.getNumParticle();
            const ScalarType temperatureNow = ringPolymer.calcClassicalKinetic() * 2.0 / ScalarType(numParticle);
            auto phase = ringPolymer.asMatrix().col(0);
            auto momentum = phase.head(numParticle);
            momentum *= sqrt(temperatureT / temperatureNow);
        }
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    void HardCore<ScalarType, NumReplica, Executor>::updateMass(RingPolymerType& ringPolymer) {
        repMass = reciprocal(ringPolymer.getMassVec());
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    void HardCore<ScalarType, NumReplica, Executor>::swap(HardCore& obj) noexcept {
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        temperatureT.swap(temperatureT);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        velocity.swap(obj.velocity);
        std::swap(maxHandleNum, obj.maxHandleNum);
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    void HardCore<ScalarType, NumReplica, Executor>::updateVelocity(const RingPolymerType& ringPolymer) {
        const auto phase = ringPolymer.asMatrix();
        const auto momentum = phase.topRows(getNumParticle());
        const size_t numReplica = getNumReplica();
        for (size_t replica = 0; replica < numReplica; ++replica) {
            auto v = velocity.col(replica);
            v = hadamard(momentum.col(replica), repMass);
        }
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    bool HardCore<ScalarType, NumReplica, Executor>::checkCollision(const RingPolymerType& ringPolymer) const {
        const size_t numReplica = getNumReplica();
        const size_t numParticle = getNumParticle();
        for (size_t replica = 0; replica < numReplica; ++replica) {
            auto phase = ringPolymer.asMatrix().col(replica);
            auto pos = phase.tail(numParticle);
            if (pos[0].isNegative()) [[unlikely]]
                return true;

            const size_t length = numParticle - 1;
            auto head = pos.head(length);
            auto tail = pos.tail(1);
            {
                using PacketType = typename Internal::BestPacket<ScalarType, Dynamic>::Type;
                size_t i = 0;
                const size_t to = length / PacketType::size() * PacketType::size();

                for (; i < to; i += PacketType::size()) {
                    const auto boolPacket = head.template packet<PacketType>(i) > tail.template packet<PacketType>(i);
                    if (horizontal_or(boolPacket)) [[unlikely]]
                        return true;
                }

                [[likely]] if (to != length) {
                    const size_t count = length - i;
                    const auto boolPacket = head.template packetPartial<PacketType>(i, count) > tail.template packetPartial<PacketType>(i, count);
                    if (horizontal_or(boolPacket)) [[unlikely]]
                        return true;
                }
            }
            if (pos[length] > latticeSize) [[unlikely]]
                return true;
        }
        return false;
    }

    template<class ScalarType, size_t NumReplica, class Executor>
    bool HardCore<ScalarType, NumReplica, Executor>::handleCollision(RingPolymerType& ringPolymer) {
        const size_t numReplica = getNumReplica();
        const size_t numParticle = getNumParticle();
        const auto& mass = ringPolymer.getMassVec();
        bool isDrifted = false;
        for (size_t replica = 0; replica < numReplica; ++replica) {
            auto phase = ringPolymer.asMatrix().col(replica);
            auto pos = phase.tail(numParticle);
            size_t i = 0;
            for (; i < numParticle - 1; ++i) {
                if (pos[i] > pos[i + 1]) {
                    const ScalarType m1 = mass[i];
                    const ScalarType m2 = mass[i + 1];
                    const ScalarType p1 = phase[i];
                    const ScalarType p2 = phase[i + 1];
                    const ScalarType next_p1 = ((m1 - m2) * p1 + ScalarType(2) * m1 * p2) * reciprocal(m1 + m2);
                    const ScalarType next_p2 = p1 + p2 - next_p1;
                    phase[i] = next_p1;
                    phase[i + 1] = next_p2;
                }
            }
            if (!pos[0].isPositive()) {
                phase[0] = abs(phase[0]);
                isDrifted = true;
            }
            if (pos[i] > latticeSize) {
                phase[i] = -abs(phase[i]);
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
