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

#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/IO/QE_scf.h"
#include "Physica/Core/Exception/SyscallException.h"
#include "Physica/Core/Parallel/Executor/ProcessExecutor.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Utils/Unix/TempFile.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType>
    class QEModel {
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using ElementTypeArray = typename Poscar::ElementTypeArray;

        std::string pathToPW;
        Utils::Array<char> input;
        ElementTypeArray elementTypes;
        unsigned int numMPIProcess;
    public:
        QEModel(const char* pathToPW_, const char* pathToInput, ElementTypeArray elementTypes_, unsigned int numMPIProcess_);
        QEModel(const QEModel&) = default;
        QEModel(QEModel&&) noexcept = default;
        ~QEModel() = default;
        /* Operators */
        QEModel& operator=(QEModel obj) noexcept;
        /* Operations */
        template<class Executor, bool IsSmallCell = true>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        void swap(QEModel& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumAtom() const noexcept { return elementTypes.getLength(); }
        [[nodiscard]] unsigned int getNumMPIProcess() const noexcept { return numMPIProcess; }
    private:
        template<size_t N> ProcessFuture run_qe(Utils::TempFile<N>& input, Utils::TempFile<N>& output) const;
    };

    template<class ScalarType, class PosScalarType>
    QEModel<ScalarType, PosScalarType>::QEModel(const char* pathToPW_, const char* pathToInput, ElementTypeArray elementTypes_, unsigned int numMPIProcess_)
            : pathToPW(pathToPW_)
            , elementTypes(std::move(elementTypes_))
            , numMPIProcess(numMPIProcess_) {
        std::ifstream fin(pathToInput);
        if (!fin)
            throw IOException("[Error]: No QE input file found");
        fin.seekg(0, std::ios::end);
        const auto size = fin.tellg();
        fin.seekg(0, std::ios::beg);
        input.resize(size + 1);
        fin.read(input.data(), size);
        input[size] = '\0';
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor, bool IsSmallCell>
    Vector<ScalarType> QEModel<ScalarType, PosScalarType>::force(const MDCellType& cell) const {
        assert(cell.getNumParticle() == getNumAtom());
        auto inputTmp = Utils::TempFile("/tmp/tmpXXXXXX");
        /* Make input */ {
            std::ofstream os(inputTmp.getName());
            os << input.data() << '\n';
            os << "CELL_PARAMETERS bohr\n";
            os << cell.getLattice() << '\n';
            os << "ATOMIC_POSITIONS bohr\n";
            for (size_t i = 0; i < getNumAtom(); ++i) {
                const auto row = cell.getPos().row(i);
                os << PhyConst<SI>::elementSymbol[elementTypes[i]] << ' '
                   << row[0] << ' '
                   << row[1] << ' '
                   << row[2] << '\n';
            }
        }
        auto outputTmp = Utils::TempFile("/tmp/tmpXXXXXX");
        run_qe(inputTmp, outputTmp).wait("[Error]: QE finished with non zero exit code");
        QE_scf out_scf(outputTmp.getName(), getNumAtom());
        return out_scf.getForce();
    }

    template<class ScalarType, class PosScalarType>
    QEModel<ScalarType, PosScalarType>& QEModel<ScalarType, PosScalarType>::operator=(QEModel obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    void QEModel<ScalarType, PosScalarType>::swap(QEModel& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pathToPW.swap(obj.pathToPW);
        input.swap(obj.input);
        elementTypes.swap(obj.elementTypes);
        std::swap(numMPIProcess, obj.numMPIProcess);
    }

    template<class ScalarType, class PosScalarType>
    template<size_t N>
    ProcessFuture QEModel<ScalarType, PosScalarType>::run_qe(Utils::TempFile<N>& input, Utils::TempFile<N>& output) const {
        int fd[2];
        if (pipe(fd) == -1)
            throw SyscallException();

        ProcessFuture future;
        future = ProcessExecutor::schedule([this, fd, &output]() {
            const int standardErr = dup(STDERR_FILENO);
            if (dup2(output.getFd(), STDOUT_FILENO) == -1) [[unlikely]]
                throw SyscallException();
            if (dup2(output.getFd(), STDERR_FILENO) == -1) [[unlikely]]
                throw SyscallException();
            close(fd[1]);
            if (dup2(fd[0], STDIN_FILENO) == -1) [[unlikely]]
                throw SyscallException();
            close(fd[0]);
            /* Execute */ {
                char numProcess[16]; // 16 is enough for unsigned int
                sprintf(numProcess, "%d", getNumMPIProcess());
                execlp("mpirun", "mpirun", "-n", numProcess, pathToPW.c_str(), "-n", numProcess, "-nk", numProcess, nullptr);
            }
            dup2(standardErr, 2);
            perror("[Error]: Failed to execute QE");
            _exit(EXIT_FAILURE);
        });

        close(fd[0]);
        std::ifstream fin(input.getName());
        if (!fin) [[unlikely]]
            throw IOException("[Error]: No QE input found, this is impossible");
        fin.seekg(0, std::ios::end);
        const auto size = fin.tellg();
        fin.seekg(0, std::ios::beg);
        Utils::Array<char> buffer(size);
        fin.read(buffer.data(), size);
        if (write(fd[1], buffer.data(), size) == -1)
            throw SyscallException();
        close(fd[1]);
        return future;
    }
}
