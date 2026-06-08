/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Physics/SolidState/CrystalCell.h"

namespace Physica {
    template<Scalar T>
    class DeepModelDataset {
        using CellType = CrystalCell<T>;
        using CellArray = Array<CellType>;
        using VectorType = VectorND<T>;
        using ForceArray = Array<VectorType>;

        CellArray cells;
        VectorType energys;
        ForceArray forces;
    public:
        DeepModelDataset() = default;
        DeepModelDataset(const DeepModelDataset&) = default;
        DeepModelDataset(DeepModelDataset&&) noexcept = default;
        ~DeepModelDataset() = default;
        /* Operators */
        DeepModelDataset& operator=(DeepModelDataset obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void append(CellType cell, T energy, VectorType force);
        H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(DeepModelDataset& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumSample() const noexcept { return cells.getLength(); }
        [[nodiscard]] const CellArray& getCells() const noexcept { return cells; }
        [[nodiscard]] const VectorType& getEnergys() const noexcept { return energys; }
        [[nodiscard]] const ForceArray& getForces() const noexcept { return forces; }
    };

    template<Scalar T>
    void DeepModelDataset<T>::append(CellType cell, T energy, VectorType force) {
        assert(cell.getNumParticle() * CellType::Dim == force.getLength());
        cells.append(std::move(cell));
        energys.append(std::move(energy));
        forces.append(std::move(force));
    }

    template<Scalar T>
    H5Group DeepModelDataset<T>::read(const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        energys.read(group, "Energys");
        
        const size_t numSample = energys.getLength();
        cells.resize(numSample);
        forces.resize(numSample);

        auto cellGroup = group.openGroup("Cells");
        auto forceGroup = group.openGroup("Forces");
        std::array<char, 32> buffer{}; //32 should be enough to hold size_t
        for (size_t i = 0; i < numSample; ++i) {
            [[maybe_unused]] const auto count = std::format_to_n(buffer.data(), buffer.size(), "%zu", i).size;
            cells[i].read(cellGroup, buffer);
            forces[i].read(forceGroup, buffer);
        }
        return group;
    }

    template<Scalar T>
    H5Group DeepModelDataset<T>::write(H5Loc& loc, const char* name) const {
        auto group = loc.openGroup(name);
        energys.write(group, "Energys");

        auto cellGroup = group.openGroup("Cells");
        auto forceGroup = group.openGroup("Forces");
        std::array<char, 32> buffer{}; //32 should be enough to hold size_t
        for (size_t i = 0; i < getNumSample(); ++i) {
            [[maybe_unused]] const auto count = std::format_to_n(buffer.data(), buffer.size(), "%zu", i).size;
            cells[i].write(cellGroup, buffer);
            forces[i].write(forceGroup, buffer);
        }
        return H5Group(group);
    }

    template<Scalar T>
    void DeepModelDataset<T>::swap(DeepModelDataset& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        cells.swap(obj.cells);
        energys.swap(obj.energys);
        forces.swap(obj.forces);
    }
}
