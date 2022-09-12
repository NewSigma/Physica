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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
     */
    template<class ScalarType>
    class Q_TIP4P {
        using HytrogenListType = Utils::Array<std::pair<size_t, size_t>>;
        HytrogenListType hytrogenList;
        ScalarType cutoff;
    public:
        Q_TIP4P(const MDCell& refer_cell, ScalarType cutoff_);
        Q_TIP4P(const Q_TIP4P&) = default;
        Q_TIP4P(Q_TIP4P&&) noexcept = default;
        ~Q_TIP4P() = default;
        /* Operators */
        Q_TIP4P& operator=(Q_TIP4P model) noexcept;
        /* Getters */
        [[nodiscard]] const HytrogenListType& getHytrogenList() const noexcept { return hytrogenList; }
        /* Helpers */
        bool checkHytrogenList() const;
        void swap(Q_TIP4P& model) noexcept;
    private:
        void makeHytrogenList(const MDCell& refer_cell);
        ScalarType minSquaredDist(const MDCell& cell, size_t from_id, size_t to_id);
    };

    template<class ScalarType>
    Q_TIP4P<ScalarType>::Q_TIP4P(const MDCell& refer_cell, ScalarType cutoff_) : cutoff(std::move(cutoff_)) {
        assert(refer_cell.getNumParticle() % 3 == 0);
        makeHytrogenList(refer_cell);
    }

    template<class ScalarType>
    Q_TIP4P<ScalarType>& Q_TIP4P<ScalarType>::operator=(Q_TIP4P model) noexcept {
        swap(model);
        return *this;
    }

    template<class ScalarType>
    bool Q_TIP4P<ScalarType>::checkHytrogenList() const {
        const size_t maxHytrogenId = hytrogenList.getLength() * 2;
        for (size_t i = 0; i < hytrogenList.getLength(); ++i) {
            const auto& pair = hytrogenList[i];
            if (pair.first == pair.second)
                return false;

            if (pair.first >= maxHytrogenId || pair.second >= maxHytrogenId)
                return false;

            for (size_t j = i + 1; j < hytrogenList.getLength(); ++j) {
                const auto& pair2 = hytrogenList[j];
                if (pair.first == pair2.first
                 || pair.first == pair2.second
                 || pair.second == pair2.first
                 || pair.second == pair2.second)
                    return false;
            }
        }
        return true;
    }

    template<class ScalarType>
    void Q_TIP4P<ScalarType>::swap(Q_TIP4P& model) noexcept {
        hytrogenList.swap(model.hytrogenList);
        cutoff.swap(model.cutoff);
    }

    template<class ScalarType>
    void Q_TIP4P<ScalarType>::makeHytrogenList(const MDCell& refer_cell) {
        const size_t numMolecule = refer_cell.getNumParticle() / 3;
        hytrogenList.resize(numMolecule);

        const auto& pos = refer_cell.getPos();
        const size_t maxHytrogenId = refer_cell.getNumParticle() - numMolecule;
        ScalarType bondLength1, bondLength2;
        size_t hytrogenId1, hytrogenId2;
        for (size_t i = 0; i < numMolecule; ++i) {
            bondLength1 = bondLength2 = std::numeric_limits<ScalarType>::max();

            auto pos_O = pos.row(maxHytrogenId + i);
            for (size_t j = 0; j < maxHytrogenId; ++j) {
                auto pos_H = pos.row(j);
                const ScalarType squared_dist = minSquaredDist(refer_cell, maxHytrogenId + i, j);
                if (squared_dist < bondLength1) {
                    if (bondLength2 > bondLength1) {
                        bondLength2 = squared_dist;
                        hytrogenId2 = j;
                    }
                    else {
                        bondLength1 = squared_dist;
                        hytrogenId1 = j;
                    }
                }
                else if (squared_dist < bondLength2) {
                    bondLength2 = squared_dist;
                    hytrogenId2 = j;
                }
            }
            auto& hytrogenPair = hytrogenList[i];
            hytrogenPair.first = hytrogenId1;
            hytrogenPair.second = hytrogenId2;
        }
    }

    template<class ScalarType>
    ScalarType Q_TIP4P<ScalarType>::minSquaredDist(const MDCell& cell, size_t from_id, size_t to_id) {
        using Vector3D = Vector<ScalarType, 3>;

        const auto& pos = cell.getPos();
        const auto& lattice = cell.getLattice();
        auto pos1 = pos.row(from_id);
        auto pos2 = pos.row(to_id);

        ScalarType result = std::numeric_limits<ScalarType>::max();
        for (int x = -1; x <= 1; ++x) {
            const Vector3D v1 = pos2.asVector() + lattice.row(0).asVector() * ScalarType(x);
            for (int y = -1; y <= 1; ++y) {
                const Vector3D v2 = v1 + lattice.row(1).asVector() * ScalarType(y);
                for (int z = -1; z <= 1; ++z) {
                    const Vector3D v3 = v2 + lattice.row(2).asVector() * ScalarType(z);
                    result = std::min(result, (v3 - pos1.asVector()).squaredNorm());
                }
            }
        }
        return result;
    }
}
