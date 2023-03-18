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

#include <cstdlib>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    template<class ScalarType>
    class MD1D {
        using VectorType = Vector<ScalarType>;

        VectorType pos;
        VectorType velocity;
        VectorType mass;
        ScalarType latticeSize;
        ScalarType timeStep;
        ScalarType collideStep;

        VectorType buffer;
    public:
        MD1D(size_t numParticle, ScalarType latticeSize_, ScalarType timeStep_, ScalarType collideStep_);
        MD1D(const MD1D&) = default;
        MD1D(MD1D&&) noexcept = default;
        ~MD1D() = default;
        /* Operators */
        MD1D& operator=(MD1D obj) noexcept;
        /* Operations */
        template<class RandomGenerator>
        void init(RandomGenerator& gen);
        void nve_step();
        void swap(MD1D& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return pos.getLength(); }
        [[nodiscard]] const VectorType& getPos() const noexcept { return pos; }
        [[nodiscard]] const VectorType& getVelocity() const noexcept { return velocity; }
    private:
        bool checkCollision() const;
        void handleCollision();
    };

    template<class ScalarType>
    MD1D<ScalarType>::MD1D(size_t numParticle, ScalarType latticeSize_, ScalarType timeStep_, ScalarType collideStep_)
            : pos(numParticle)
            , velocity(numParticle)
            , mass(numParticle)
            , latticeSize(latticeSize_)
            , timeStep(timeStep_)
            , collideStep(collideStep_)
            , buffer(numParticle) {
        assert(timeStep > collideStep);
        assert(collideStep.isPositive());
    }

    template<class ScalarType>
    MD1D<ScalarType>& MD1D<ScalarType>::operator=(MD1D obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    template<class RandomGenerator>
    void MD1D<ScalarType>::init(RandomGenerator& gen) {
        using TrivialType = typename ScalarType::TrivialType;
        std::uniform_real_distribution<TrivialType> uniform_dist{};
        std::normal_distribution<TrivialType> normal_dist{};
        for (size_t i = 0; i < getNumParticle(); ++i) {
            pos[i] = uniform_dist(gen);
            velocity[i] = normal_dist(gen);
            mass[i] = uniform_dist(gen);
        }
        std::sort(pos.begin(), pos.end());
    }

    template<class ScalarType>
    void MD1D<ScalarType>::nve_step() {
        ScalarType from = 0;
        ScalarType lStep = 0;
        ScalarType rStep = timeStep;
        ScalarType to = timeStep;
        buffer = pos;
        while (lStep != timeStep) {
            const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
            if (isDeltaSmallEnough) {
                pos = buffer + velocity * (rStep - from);
                buffer += velocity * (lStep - from);
                from = lStep;
                to = timeStep;
                rStep = timeStep;
                handleCollision();
            }

            const ScalarType step = to - from;
            pos = buffer + velocity * step;
            if (checkCollision())
                rStep = to;
            else
                lStep = to;
            to = (lStep + rStep) * 0.5;
        }
    }

    template<class ScalarType>
    void MD1D<ScalarType>::swap(MD1D& obj) noexcept {
        pos.swap(obj.pos);
        velocity.swap(obj.velocity);
        mass.swap(obj.mass);
        latticeSize.swap(obj.latticeSize);
        timeStep.swap(obj.timeStep);
        collideStep.swap(obj.collideStep);
        buffer.swap(obj.buffer);
    }

    template<class ScalarType>
    bool MD1D<ScalarType>::checkCollision() const {
        bool isHappened = pos[0].isPositive();

        size_t i = 0;
        for (; i < getNumParticle() - 1; ++i)
            isHappened &= pos[i] < pos[i + 1];
        isHappened &= pos[i] < latticeSize;
        return !isHappened;
    }

    template<class ScalarType>
    void MD1D<ScalarType>::handleCollision() {
        if (!pos[0].isPositive())
            velocity[0].toOpposite();

        size_t i = 0;
        for (; i < getNumParticle() - 1; ++i) {
            if (pos[i] >= pos[i + 1]) {
                const ScalarType m1 = mass[i];
                const ScalarType m2 = mass[i + 1];
                const ScalarType v1 = velocity[i];
                const ScalarType v2 = velocity[i + 1];
                const ScalarType next_v1 = ((m1 - m2) * v1 + ScalarType(2) * m2 * v2) * reciprocal(m1 + m2);
                const ScalarType next_v2 = next_v1 + v1 - v2;
                velocity[i] = next_v1;
                velocity[i + 1] = next_v2;
            }
        }

        if (pos[i] >= latticeSize)
            velocity[i].toOpposite();
    }
}
