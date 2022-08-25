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
#include <iostream>
#include <cassert>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <fcntl.h>
#include "Physica/Core/Interface/VaspWarpper.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Utils/Unix/UnixHelper.h"

namespace Physica::Core {
    VaspWarpper::VaspWarpper()
            : Base([this]() { this->run(); })
            , vaspWorkingDir()
            , logFilePath()
            , core() {}

    VaspWarpper::VaspWarpper(size_t core_,
                             std::string pathToVasp_,
                             std::string workingDir,
                             std::string logFilePath_,
                             Poscar poscar_)
            : Base([this]() { this->run(); })
            , pathToVasp(std::move(pathToVasp_))
            , vaspWorkingDir(std::move(workingDir))
            , logFilePath(std::move(logFilePath_))
            , core(core_)
            , poscar(std::move(poscar_)) {
        std::ofstream fout((vaspWorkingDir + std::string("/POSCAR")).c_str(), std::ios_base::out | std::ios_base::trunc);
        fout << poscar;
    }

    VaspWarpper::VaspWarpper(VaspWarpper&& vasp) noexcept
            : SubProcess(std::move(vasp))
            , pathToVasp(std::move(vasp.pathToVasp))
            , vaspWorkingDir(std::move(vasp.vaspWorkingDir))
            , logFilePath(std::move(vasp.logFilePath))
            , core(vasp.core) {}

    VaspWarpper& VaspWarpper::operator=(VaspWarpper vasp) noexcept {
        swap(vasp);
        return *this;
    }

    void VaspWarpper::run() const {
        int standardErr = dup(2);
        if (!logFilePath.empty()) {
            int log_fd = open(logFilePath.c_str()
                    , O_WRONLY | O_TRUNC | O_CREAT
                    , S_IRUSR | S_IWUSR);
            dup2(log_fd, 1);
            dup2(log_fd, 2);
        }
        [[maybe_unused]] int err = chdir(vaspWorkingDir.c_str());
        /* Execute */ {
            constexpr size_t bufferLength = 21;
            char coreStr[bufferLength];
            [[maybe_unused]] int count = sprintf(coreStr, "%ld", core);
            assert(0 <= count && static_cast<size_t>(count) < bufferLength);
            execlp("mpirun", "mpirun", "-n", coreStr, pathToVasp.c_str(), nullptr);
        }
        dup2(standardErr, 2);
        perror("[Error]: Failed to execute VASP");
        _exit(EXIT_FAILURE);
    }

    typename VaspWarpper::ScalarType VaspWarpper::getEnergy() const {
        const char* tempfile = tmpnam(nullptr);
        const std::string command = std::string("grep energy ") +
                                    vaspWorkingDir +
                                    std::string("/OUTCAR | tail -n 1 | tr -s ' ' | cut -d ' ' -f 5 >") +
                                    std::string(tempfile);
        [[maybe_unused]] int err = system(command.c_str());
        ScalarType result;
        /* Read data */ {
            std::ifstream fin(tempfile);
            fin >> result;
        }
        unlink(tempfile);
        return result;
    }

    typename VaspWarpper::ScalarType VaspWarpper::getPress() const {
        const char* tempfile = tmpnam(nullptr);
        const std::string command = std::string("grep 'in kB' ") +
                                    vaspWorkingDir +
                                    std::string("/OUTCAR | tr -s ' ' | cut -d ' ' -f 4,5,6 >") +
                                    std::string(tempfile);
        [[maybe_unused]] int err = system(command.c_str());
        ScalarType press_x, press_y, press_z;
        /* Read data */ {
            std::ifstream fin(tempfile);
            fin >> press_x >> press_y >> press_z;
        }
        unlink(tempfile);
        return (press_x + press_y + press_z) / 3.0f;
    }

    Vector<typename VaspWarpper::ScalarType> VaspWarpper::getForce() const {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dynamic, 6>;

        std::ifstream fin((vaspWorkingDir + std::string("/OUTCAR")).c_str());
        if (fin) {
            fin.seekg(0, std::ios::end);
            const auto size = fin.tellg();
            fin.seekg(0, std::ios::beg);
            Utils::Array<char> buffer(size);
            std::string str{};
            Vector<ScalarType> force(3 * poscar.getAtomCount());
            do {
                fin.getline(buffer.data(), size);
                str = buffer.data();
                const bool success = str.find("TOTAL-FORCE") != std::string::npos;
                if (success) {
                    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    MatrixType pos_force(poscar.getAtomCount(), 6);
                    fin >> pos_force;
                    size_t index = 0;
                    for (size_t r = 0; r < pos_force.getRow(); ++r) {
                        for (size_t c = 3; c < pos_force.getColumn(); ++c) {
                            force[index] = pos_force(r, c);
                            ++index;
                        }
                    }
                    break;
                }
            } while(bool(fin));

            if (!fin)
                throw BadFileFormatException();
            return force;
        }
        else
            throw IOException();
    }

    void VaspWarpper::swap(VaspWarpper& vasp) noexcept {
        pathToVasp.swap(vasp.pathToVasp);
        vaspWorkingDir.swap(vasp.vaspWorkingDir);
        logFilePath.swap(vasp.logFilePath);
        std::swap(core, vasp.core);
        poscar.swap(vasp.poscar);
    }
}
