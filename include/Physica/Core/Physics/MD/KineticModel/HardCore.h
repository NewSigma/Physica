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
    template<class ScalarType, class PosScalarType, unsigned int Dim> class RingPolymer;

    template<class ScalarType>
    class HardCore {
        using RingPolymerType = RingPolymer<ScalarType, ScalarType, 1>;

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
    private:
        bool checkCollision(const RingPolymerType& ringPolymer) const;
        void handleCollision(RingPolymerType& ringPolymer);
    };

    template<class ScalarType>
    HardCore<ScalarType>::HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , repMass(numParticle, 0)
            , buffer(numParticle)
            , velocity(numParticle) {
        assert(collideFactor < ScalarType(1.0) && collideFactor.isPositive());
    }

    template<class ScalarType>
    HardCore<ScalarType>& HardCore<ScalarType>::operator=(HardCore<ScalarType> obj) noexcept {
        swap(*this);
        return *this;
    }

    template<class ScalarType>
    void HardCore<ScalarType>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT) {
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
                handleCollision(ringPolymer);
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

    template<class ScalarType>
    void HardCore<ScalarType>::updateMass(RingPolymerType& ringPolymer) {
        repMass = reciprocal(ringPolymer.getMassVec());
    }

    template<class ScalarType>
    void HardCore<ScalarType>::swap(HardCore& obj) noexcept {
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        velocity.swap(obj.velocity);
    }

    template<class ScalarType>
    bool HardCore<ScalarType>::checkCollision(const RingPolymerType& ringPolymer) const {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto pos = phase.tail(numParticle);
        bool isHappened = pos[0].isPositive();

        size_t i = 0;
        for (; i < numParticle - 1; ++i)
            isHappened &= pos[i] < pos[i + 1];
        isHappened &= pos[i] < latticeSize;
        return !isHappened;
    }

    template<class ScalarType>
    void HardCore<ScalarType>::handleCollision(RingPolymerType& ringPolymer) {
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
}
