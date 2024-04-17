/*
 * Copyright 2024 WeiBo He.
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

#include "SpinlessElectron.h"

namespace Physica::Core {
    class SpinElectron {
        SpinlessElectron spinUp;
        SpinlessElectron spinDown;
    public:
        SpinElectron() = default;
        SpinElectron(SpinlessElectron spinUp_, SpinlessElectron spinDown_);
        SpinElectron(const SpinElectron&) = default;
        SpinElectron(SpinElectron&&) noexcept = default;
        ~SpinElectron() = default;
        /* Operators */
        SpinElectron& operator=(SpinElectron obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(SpinElectron other) const noexcept { return spinUp == other.spinUp && spinDown == other.spinDown; }
        /* Operations */
        [[nodiscard]] SpinElectron hopUp(unsigned char from, unsigned char to) const;
        [[nodiscard]] SpinElectron hopDown(unsigned char from, unsigned char to) const;
        void swap(SpinElectron& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] SpinlessElectron getSpinUp() const noexcept { return spinUp; }
        [[nodiscard]] SpinlessElectron getSpinDown() const noexcept { return spinDown; }
        [[nodiscard]] bool isVacuum() const noexcept { return spinUp.isVacuum() && spinDown.isVacuum(); }
        [[nodiscard]] bool isUpOccupy(unsigned char site) const noexcept { return spinUp.isOccupy(site); }
        [[nodiscard]] bool isDownOccupy(unsigned char site) const noexcept { return spinDown.isOccupy(site); }
        [[nodiscard]] unsigned int getNumSpinUpElectron() const noexcept { return spinUp.getNumElectron(); }
        [[nodiscard]] unsigned int getNumSpinDownElectron() const noexcept { return spinDown.getNumElectron(); }
    };

    SpinElectron::SpinElectron(SpinlessElectron spinUp_, SpinlessElectron spinDown_)
            : spinUp(spinUp_), spinDown(spinDown_) {}

    SpinElectron SpinElectron::hopUp(unsigned char from, unsigned char to) const {
        auto newSpinUp = spinUp.hop(from, to);
        const bool hopFailed = newSpinUp.isVacuum() && !spinUp.isVacuum();
        if (hopFailed)
            return SpinElectron();
        return SpinElectron(std::move(newSpinUp), spinDown);
    }

    SpinElectron SpinElectron::hopDown(unsigned char from, unsigned char to) const {
        auto newSpinDown = spinDown.hop(from, to);
        const bool hopFailed = newSpinDown.isVacuum() && !spinDown.isVacuum();
        if (hopFailed)
            return SpinElectron();
        return SpinElectron(spinUp, std::move(newSpinDown));
    }

    void SpinElectron::swap(SpinElectron& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spinUp.swap(obj.spinUp);
        spinDown.swap(obj.spinDown);
    }
}
