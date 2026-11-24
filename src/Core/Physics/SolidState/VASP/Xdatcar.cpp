/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/SolidState/VASP/Xdatcar.h"

using namespace Physica;

Xdatcar::Xdatcar(std::ifstream fin_)
        : data()
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

void Xdatcar::swap(Xdatcar& __restrict xdatcar) noexcept {
    assert(this != &xdatcar && "[Error]: Self swap is likely a bug");
    data.swap(xdatcar.data);
    fin.swap(xdatcar.fin);
    std::swap(stepNum, xdatcar.stepNum);
    std::swap(init, xdatcar.init);
}
