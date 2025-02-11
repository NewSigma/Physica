/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/IO/VASP/Poscar.h"
#include "Physica/Core/IO/QE/PWscfOut.h"
#include "Physica/Core/Exception/SystemException.h"
#include "Physica/Core/Parallel/Executor/ProcessExecutor.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Utils/Unix/TempFile.h"

namespace Physica {
    template<Scalar T>
    class QEModel {
        using MDCellType = MDCell<T>;
        using ElementTypeArray = Poscar<T>::ElementTypeArray;

        std::string pathToPW;
        Array<char> input;
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
        template<class Executor>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const;
        template<Vector V, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<V>& result) const noexcept;
        template<class Executor>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const { return VectorND<T>(cell.getDOF(), 0); }
        void swap(QEModel& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return elementTypes.getLength(); }
        [[nodiscard]] unsigned int getNumMPIProcess() const noexcept { return numMPIProcess; }
    private:
        template<size_t N> ProcessFuture run_qe(TempFile<N>& __restrict input, TempFile<N>& __restrict output) const;
    };

    template<Scalar T>
    QEModel<T>::QEModel(const char* pathToPW_, const char* pathToInput, ElementTypeArray elementTypes_, unsigned int numMPIProcess_)
            : pathToPW(pathToPW_)
            , elementTypes(std::move(elementTypes_))
            , numMPIProcess(numMPIProcess_) {
        std::ifstream fin(pathToInput);
        if (!fin)
            throw IOException("[Error]: No QE input file found");
        fin.seekg(0, std::ios::end);
        const std::streamsize size = fin.tellg();
        fin.seekg(0, std::ios::beg);
        input.resize(size + 1);
        fin.read(input.data(), size);
        input[size] = '\0';
    }

    template<Scalar T>
    template<class Executor>
    VectorND<T> QEModel<T>::force(const MDCellType& cell) const {
        VectorND<T> result(cell.getNumParticle());
        forceAsync<VectorND<T>, Executor>(cell, result);
        return result;
    }

    template<Scalar T>
    template<Vector V, class Executor>
    void QEModel<T>::forceAsync(const MDCellType& cell, ContinuousVector<V>& result) const noexcept {
        assert(cell.getNumParticle() == getNumParticle());
        try {
            auto inputTmp = TempFile("/tmp/tmpXXXXXX");
            /* Make input */ {
                std::ofstream os(inputTmp.getName());
                os << input.data() << '\n';
                os << "CELL_PARAMETERS bohr\n";
                os << cell.getLattice().format() << '\n';
                os << "ATOMIC_POSITIONS bohr\n";
                for (size_t i = 0; i < getNumParticle(); ++i) {
                    const auto row = cell.getPos().row(i);
                    os << PhyConst<SI>::elementSymbol[elementTypes[i]] << ' '
                    << row[0] << ' '
                    << row[1] << ' '
                    << row[2] << '\n';
                }
            }
            auto outputTmp = TempFile("/tmp/tmpXXXXXX");
            const int err = run_qe(inputTmp, outputTmp).wait();
            if (err != 0) [[unlikely]] {
                inputTmp.release();
                outputTmp.release();
                throw std::runtime_error("[Error]: QE finished with non zero exit code");
            }
            PWscfOut out_scf(outputTmp.getName(), getNumParticle());
            result.getDerived() = out_scf.getForce();
        }
        catch (std::exception& e) {
            fprintf(stderr, "%s\n", e.what());
            exit(EXIT_FAILURE);
        }
    }

    template<Scalar T>
    QEModel<T>& QEModel<T>::operator=(QEModel obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void QEModel<T>::swap(QEModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pathToPW.swap(obj.pathToPW);
        input.swap(obj.input);
        elementTypes.swap(obj.elementTypes);
        std::swap(numMPIProcess, obj.numMPIProcess);
    }

    template<Scalar T>
    template<size_t N>
    ProcessFuture QEModel<T>::run_qe(
            TempFile<N>& __restrict input, TempFile<N>& __restrict output) const {
        int fd[2];
        if (pipe(fd) == -1)
            throw SystemException();

        ProcessFuture future;
        future = ProcessExecutor::schedule([this, fd, &output]() {
            const int standardErr = dup(STDERR_FILENO);
            if (dup2(output.getFd(), STDOUT_FILENO) == -1) [[unlikely]]
                throw SystemException();
            if (dup2(output.getFd(), STDERR_FILENO) == -1) [[unlikely]]
                throw SystemException();
            close(fd[1]);
            if (dup2(fd[0], STDIN_FILENO) == -1) [[unlikely]]
                throw SystemException();
            close(fd[0]);
            /* Execute */ {
                constexpr int bufferLength = 16; // 16 is enough for unsigned int
                char numProcess[bufferLength];
                [[maybe_unused]] const int count = sprintf(numProcess, "%d", getNumMPIProcess());
                assert(0 <= count && count < bufferLength && "[Error]: Unexpected bad printf");
                execlp("mpirun", "mpirun", "-n", numProcess, pathToPW.c_str(), "-n", numProcess, "-nk", numProcess, nullptr);
            }
            dup2(standardErr, STDERR_FILENO);
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
        Array<char> buffer(size);
        fin.read(buffer.data(), size);
        if (write(fd[1], buffer.data(), size) == -1)
            throw SystemException();
        close(fd[1]);
        return future;
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<QEModel<T>> {
    public:
        constexpr static bool IsPeriodBoundary = true;
        constexpr static bool IsContractable = false;
    };
}
