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
#include "Physica/Core/IO/Xdatcar.h"

namespace Physica::Core {
    Xdatcar::Xdatcar(std::ifstream fin_) : data()
                                        , fin(std::move(fin_))
                                        , stepNum(0)
                                        , init(false) {}
    Xdatcar& Xdatcar::operator=(Xdatcar xdatcar) noexcept {
        swap(xdatcar);
        return *this;
    }
    /**
     * \returns true if read successfully
     */
    bool Xdatcar::step() {
        if (!init) {
            fin >> data;
            init = true;
        }
        else {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            data.readAtomPos(fin);
        }
        ++stepNum;
        return bool(fin);
    }

    void Xdatcar::swap(Xdatcar& xdatcar) noexcept {
        data.swap(xdatcar.data);
    }
}

namespace std {
    template<>
    inline void swap<Physica::Core::Xdatcar>(
            Physica::Core::Xdatcar& car1, Physica::Core::Xdatcar& car2) noexcept {
        car1.swap(car2);
    }
}
