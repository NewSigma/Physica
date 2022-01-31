/*
 * Copyright 2022 WeiBo He.
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
    /**
     * \class KPointBase contains members that do not depends on spin polarization
     */
    template<class ScalarType>
    class KPointBase {
    public:
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using OccupacyVector = Vector<Scalar<Float, false>>;
    private:
        Vector3D pos;
        ScalarType weight;
    public:
        KPointBase() = default;
        KPointBase(Vector3D pos_, ScalarType weight_);
        KPointBase(const KPointBase&) = default;
        KPointBase(KPointBase&&) noexcept = default;
        ~KPointBase() = default;
        /* Operators */
        KPointBase& operator=(const KPointBase&) = delete;
        KPointBase& operator=(KPointBase&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] const Vector3D& getPos() const noexcept { return pos; }
        [[nodiscard]] const ScalarType& getWeight() const noexcept { return weight; }
        /* Helpers */
        void swap(KPointBase& base) noexcept;
    };

    template<class ScalarType>
    KPointBase<ScalarType>::KPointBase(Vector3D pos_, ScalarType weight_) : pos(std::move(pos_)), weight(std::move(weight_)) {}

    template<class ScalarType>
    void KPointBase<ScalarType>::swap(KPointBase& base) noexcept {
        pos.swap(base.pos);
        weight.swap(base.weight);
    }
}
