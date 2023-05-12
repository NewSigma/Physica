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
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica> class RingPolymer;

    template<class ScalarType, class Executor = SequentialExecutor>
    class HardCore {
    public:
        using RingPolymerType = RingPolymer<ScalarType, ScalarType, 1, 1>;
    private:
        ScalarType latticeSize;
        ScalarType collideFactor;
        Vector<ScalarType> repMass;
        Vector<ScalarType> buffer;
        Vector<ScalarType> velocity;
    public:
        HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle);
        HardCore(const HardCore&) = default;
        HardCore(HardCore&&) noexcept = default;
        ~HardCore() = default;
        /* Operators */
        HardCore& operator=(HardCore obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT);
        void updateMass(RingPolymerType& ringPolymer);
        void swap(HardCore& obj) noexcept;
        /* Static members */
        static void handleCollision(ScalarType latticeSize, RingPolymerType& ringPolymer);
    private:
        bool checkCollision(const RingPolymerType& ringPolymer) const;
    };

    template<class ScalarType, class Executor>
    HardCore<ScalarType, Executor>::HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , repMass(numParticle, 0)
            , buffer(numParticle)
            , velocity(numParticle) {
        assert(collideFactor < ScalarType(1.0) && collideFactor.isPositive());
    }

    template<class ScalarType, class Executor>
    HardCore<ScalarType, Executor>& HardCore<ScalarType, Executor>::operator=(HardCore<ScalarType, Executor> obj) noexcept {
        swap(*this);
        return *this;
    }

    template<class ScalarType, class Executor>
    void HardCore<ScalarType, Executor>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT) {
        const size_t numParticle = ringPolymer.getNumParticle();
        const ScalarType collideStep = collideFactor * deltaT;
        auto phase = ringPolymer.asMatrix().col(0);
        auto momentum = phase.head(numParticle);
        auto pos = phase.tail(numParticle);
        assert(numParticle == repMass.getLength());

        ScalarType lStep = 0;
        ScalarType rStep = deltaT;
        ScalarType from = 0;
        ScalarType to = deltaT;
        buffer = pos;
        velocity = hadamard(momentum, repMass);
        while (lStep != deltaT) {
            const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
            if (isDeltaSmallEnough) {
                pos = buffer + velocity * (rStep - from);
                buffer += velocity * (lStep - from);
                from = lStep;
                to = deltaT;
                rStep = deltaT;
                handleCollision(latticeSize, ringPolymer);
                velocity = hadamard(momentum, repMass);
            }

            const ScalarType step = to - from;
            pos = buffer + velocity * step;
            if (checkCollision(ringPolymer))
                rStep = to;
            else
                lStep = to;
            to = (lStep + rStep) * 0.5;
        }
    }

    template<class ScalarType, class Executor>
    void HardCore<ScalarType, Executor>::updateMass(RingPolymerType& ringPolymer) {
        repMass = reciprocal(ringPolymer.getMassVec());
    }

    template<class ScalarType, class Executor>
    void HardCore<ScalarType, Executor>::swap(HardCore& obj) noexcept {
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        velocity.swap(obj.velocity);
    }

    template<class ScalarType, class Executor>
    void HardCore<ScalarType, Executor>::handleCollision(ScalarType latticeSize, RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto pos = phase.tail(numParticle);
        if (!pos[0].isPositive())
            phase[0].toOpposite();

        const auto& mass = ringPolymer.getMassVec();
        size_t i = 0;
        for (; i < numParticle - 1; ++i) {
            if (pos[i] >= pos[i + 1]) {
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

        if (pos[i] >= latticeSize)
            phase[i].toOpposite();
    }

    template<class ScalarType, class Executor>
    bool HardCore<ScalarType, Executor>::checkCollision(const RingPolymerType& ringPolymer) const {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto pos = phase.tail(numParticle);
        bool isGood = pos[0].isPositive();

        size_t i = 0;
        for (; i < numParticle - 1; ++i)
            isGood &= pos[i] < pos[i + 1];
        isGood &= pos[i] < latticeSize;
        return !isGood;
    }
}
