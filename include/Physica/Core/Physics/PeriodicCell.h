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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"
#include "ReciprocalCell.h"

namespace Physica::Core {
    namespace Internal {
        class PeriodicCellImpl {
        public:
            enum class Type : bool {
                Direct,
                Cartesian
            };
        };
    }

    template<class ScalarType, unsigned int Dim>
    class PeriodicCell : public Internal::PeriodicCellImpl {
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid ScalarType");
        static_assert(Dim == 1 || Dim == 2 || Dim == 3, "[Error]: Unsupported dimention");
    public:
        using LatticeMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dim, Dim>;
        using PositionMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dynamic, Dim>;
        using MomentumMatrix = PositionMatrix;
        using SearchRangeType = Utils::Array<ssize_t, Dim>;
    protected:
        using InvLatticeMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Element, Dim, Dim>;
        using VectorType = Vector<ScalarType, Dim>;

        LatticeMatrix lattice;
        PositionMatrix pos;
        Type type;
    public:
        PeriodicCell();
        PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_, Type type_);
        PeriodicCell(const PeriodicCell&) = default;
        PeriodicCell(PeriodicCell&&) noexcept = default;
        ~PeriodicCell() = default;
        /* Operators */
        PeriodicCell& operator=(PeriodicCell cell) noexcept;
        /* Operations */
        [[nodiscard]] VectorType minDistVector(size_t id_from, size_t id_to) const;
        [[nodiscard]] VectorType minDistVector(VectorType from, size_t id_to) const;
        void normalize();
        void scale(ScalarType factor);
        void toDirect() { toDirect(makeInvLattice()); }
        void toCartesian();
        void unitToSuper(unsigned int x, unsigned int y, unsigned int z);
        void superToUnit(unsigned int x, unsigned int y, unsigned int z);
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] const PositionMatrix& getPos() const noexcept { return pos; }
        [[nodiscard]] Type getType() const noexcept { return type; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return pos.getRow(); }
        [[nodiscard]] ReciprocalCell<ScalarType> reciprocal() const { return ReciprocalCell(lattice); }
        [[nodiscard]] ScalarType getVolume() const noexcept { return getVolume(lattice); }
        /* Setters */
        void setLattice(LatticeMatrix new_lattice) { lattice = new_lattice; }
        void setPos(PositionMatrix new_pos) noexcept { pos = std::move(new_pos); }
        /* Helper */
        void swap(PeriodicCell& cell) noexcept;
        /* Static members */
        [[nodiscard]] static ScalarType getVolume(const LatticeMatrix& lattice);
        static void toDirect(PositionMatrix& target, const LatticeMatrix& lattice);
        static void toCartesian(PositionMatrix& target, const LatticeMatrix& lattice);
        [[nodiscard]] static SearchRangeType estimateRange(const LatticeMatrix& cell, ScalarType cutoff);
        template<class Functor> static void forCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func);
        template<class Functor> static void forReducedCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func);
    protected:
        [[nodiscard]] InvLatticeMatrix makeInvLattice() const noexcept { return lattice.inverse(); }
        void toDirect(const InvLatticeMatrix& invLattice);
        static void toDirect(PositionMatrix& target, const InvLatticeMatrix& invLattice);
        void normalize_direct();
        void scale_direct(ScalarType factor);
        void scale_cartesian(ScalarType factor);
        static void unitToSuper_direct(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z);
        static void superToUnit_direct(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z);
        void unitToSuper_direct(unsigned int x, unsigned int y, unsigned int z);
        void superToUnit_direct(unsigned int x, unsigned int y, unsigned int z);
    };

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>::PeriodicCell()
            : lattice(LatticeMatrix::unitMatrix(Dim))
            , pos()
            , type(Type::Direct) {}

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>::PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_, Type type_)
            : lattice(std::move(lattice_))
            , pos(std::move(pos_))
            , type(type_) {}

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>& PeriodicCell<ScalarType, Dim>::operator=(PeriodicCell cell) noexcept {
        swap(cell);
        return *this;
    }

    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::VectorType
    PeriodicCell<ScalarType, Dim>::minDistVector(size_t id_from, size_t id_to) const {
        ScalarType record_dist = std::numeric_limits<ScalarType>::max();
        VectorType result{};
        VectorType delta = pos.row(id_to).asVector() - pos.row(id_from).asVector();
        if (type == Type::Direct)
            delta = lattice.transpose() * delta;
        if constexpr (Dim == 1) {
            VectorType v1;
            for (int x = -1; x <= 1; ++x) {
                const bool isSelf = id_from == id_to && x == 0;
                if (isSelf)
                    continue; // We are not interested to distance between particle and itself
                v1 = delta + lattice.row(0).asVector() * ScalarType(x);
                const ScalarType squared_norm = v1.squaredNorm();
                if (squared_norm < record_dist) {
                    record_dist = squared_norm;
                    result = v1;
                }
            }
        }
        else if constexpr (Dim == 2) {
            VectorType v1, v2;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0).asVector() * ScalarType(x);
                for (int y = -1; y <= 1; ++y) {
                    const bool isSelf = id_from == id_to && x == 0 && y == 0;
                    if (isSelf)
                        continue; // We are not interested to distance between particle and itself
                    v2 = v1 + lattice.row(1).asVector() * ScalarType(y);
                    const ScalarType squared_norm = v2.squaredNorm();
                    if (squared_norm < record_dist) {
                        record_dist = squared_norm;
                        result = v2;
                    }
                }
            }
        }
        else if constexpr (Dim == 3) {
            VectorType v1, v2, v3;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0).asVector() * ScalarType(x);
                for (int y = -1; y <= 1; ++y) {
                    v2 = v1 + lattice.row(1).asVector() * ScalarType(y);
                    for (int z = -1; z <= 1; ++z) {
                        const bool isSelf = id_from == id_to && x == 0 && y == 0 && z == 0;
                        if (isSelf)
                            continue; // We are not interested to distance between particle and itself
                        v3 = v2 + lattice.row(2).asVector() * ScalarType(z);
                        const ScalarType squared_norm = v3.squaredNorm();
                        if (squared_norm < record_dist) {
                            record_dist = squared_norm;
                            result = v3;
                        }
                    }
                }
            }
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::VectorType
    PeriodicCell<ScalarType, Dim>::minDistVector(VectorType from, size_t id_to) const {
        ScalarType record_dist = std::numeric_limits<ScalarType>::max();
        VectorType result{};
        VectorType delta;
        if (type == Type::Direct)
            delta = lattice.transpose() * pos.row(id_to).asVector();
        else
            delta = pos.row(id_to).asVector();
        delta -= from;

        if constexpr (Dim == 1) {
            VectorType v1;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0).asVector() * ScalarType(x);
                const ScalarType squared_norm = v1.squaredNorm();
                if (squared_norm < record_dist) {
                    record_dist = squared_norm;
                    result = v1;
                }
            }
        }
        else if constexpr (Dim == 2) {
            VectorType v1, v2;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0).asVector() * ScalarType(x);
                for (int y = -1; y <= 1; ++y) {
                    v2 = v1 + lattice.row(1).asVector() * ScalarType(y);
                    const ScalarType squared_norm = v2.squaredNorm();
                    if (squared_norm < record_dist) {
                        record_dist = squared_norm;
                        result = v2;
                    }
                }
            }
        }
        else if constexpr (Dim == 3) {
            VectorType v1, v2, v3;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0).asVector() * ScalarType(x);
                for (int y = -1; y <= 1; ++y) {
                    v2 = v1 + lattice.row(1).asVector() * ScalarType(y);
                    for (int z = -1; z <= 1; ++z) {
                        v3 = v2 + lattice.row(2).asVector() * ScalarType(z);
                        const ScalarType squared_norm = v3.squaredNorm();
                        if (squared_norm < record_dist) {
                            record_dist = squared_norm;
                            result = v3;
                        }
                    }
                }
            }
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::normalize() {
        if (type == Type::Cartesian) {
            toDirect(makeInvLattice());
            normalize_direct();
            toCartesian();
        }
        else
            normalize_direct();
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::scale(ScalarType factor) {
        if (type == Type::Direct)
            scale_direct(factor);
        else
            scale_cartesian(factor);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::unitToSuper(unsigned int x, unsigned int y, unsigned int z) {
        if (type == Type::Cartesian) {
            toDirect();
            unitToSuper_direct(x, y, z);
            toCartesian();
        }
        else
            unitToSuper_direct(x, y, z);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::superToUnit(unsigned int x, unsigned int y, unsigned int z) {
        if (type == Type::Cartesian) {
            toDirect();
            superToUnit_direct(x, y, z);
            toCartesian();
        }
        else
            superToUnit_direct(x, y, z);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toCartesian() {
        if (type == Type::Direct) {
            pos *= lattice;
            type = Type::Cartesian;
        }
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::swap(PeriodicCell& cell) noexcept {
        lattice.swap(cell.lattice);
        pos.swap(cell.pos);
        std::swap(type, cell.type);
    }

    template<class ScalarType, unsigned int Dim>
    ScalarType PeriodicCell<ScalarType, Dim>::getVolume(const LatticeMatrix& lattice) {
        if constexpr (Dim == 1)
            return abs(lattice(0, 0));
        else if constexpr (Dim == 2)
            return (lattice.row(0).crossProduct(lattice.row(1))).compute().norm();
        else
            return abs((lattice.row(0).crossProduct(lattice.row(1))).compute() * lattice.row(2).asVector());
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toDirect(PositionMatrix& target, const LatticeMatrix& lattice) {
        const InvLatticeMatrix inv = lattice.inverse();
        toDirect(target, inv);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toCartesian(PositionMatrix& target, const LatticeMatrix& lattice) {
        target *= lattice;
    }

    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::SearchRangeType
    PeriodicCell<ScalarType, Dim>::estimateRange(const LatticeMatrix& lattice, ScalarType cutoff) {
        const ReciprocalCell repCell(lattice);
        const auto& repLatt = repCell.getLattice();
        const ScalarType factor = cutoff * (1 / (2 * M_PI));
        SearchRangeType range{};
        if constexpr (Dim == 1) {
            range[0] = static_cast<ssize_t>(double(factor * repLatt.row(0).norm()) + 1);
        }
        else if constexpr (Dim == 2) {
            range[0] = static_cast<ssize_t>(double(factor * repLatt.row(0).norm()) + 1);
            range[1] = static_cast<ssize_t>(double(factor * repLatt.row(1).norm()) + 1);
        }
        else if (Dim == 3) {
            range[0] = static_cast<ssize_t>(double(factor * repLatt.row(0).norm()) + 1);
            range[1] = static_cast<ssize_t>(double(factor * repLatt.row(1).norm()) + 1);
            range[2] = static_cast<ssize_t>(double(factor * repLatt.row(2).norm()) + 1);
        }
        return range;
    }

    template<class ScalarType, unsigned int Dim>
    template<class Functor>
    void PeriodicCell<ScalarType, Dim>::forCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func) {
        if constexpr (Dim == 1) {
            auto a1 = lattice.row(0);

            VectorType v1;
            for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                v1 = ScalarType(x) * a1.asVector();
                func(v1);
            }
        }
        else if constexpr (Dim == 2) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);

            VectorType v1, v2;
            for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                v1 = ScalarType(x) * a1.asVector();
                for (ssize_t y = -range[1]; y <= range[1]; ++y) {
                    v2 = v1 + ScalarType(y) * a2.asVector();
                    func(v2);
                }
            }
        }
        else if constexpr (Dim == 3) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);
            auto a3 = lattice.row(2);

            VectorType v1, v2, v3;
            for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                v1 = ScalarType(x) * a1.asVector();
                for (ssize_t y = -range[1]; y <= range[1]; ++y) {
                    v2 = v1 + ScalarType(y) * a2.asVector();
                    for (ssize_t z = -range[2]; z <= range[2]; ++z) {
                        v3 = v2 + ScalarType(z) * a3.asVector();
                        func(v3);
                    }
                }
            }
        }
    }

    template<class ScalarType, unsigned int Dim>
    template<class Functor>
    void PeriodicCell<ScalarType, Dim>::forReducedCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func) {
        if constexpr (Dim == 1) {
            auto a1 = lattice.row(0);

            VectorType v1;
            for (ssize_t x = 0; x <= range[0]; ++x) {
                v1 = ScalarType(x) * a1.asVector();
                func(v1);
            }
        }
        else if constexpr (Dim == 2) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);

            VectorType v1, v2;
            for (ssize_t x = 0; x <= range[0]; ++x) {
                v1 = ScalarType(x) * a1.asVector();
                const ssize_t minY = x == 0 ? 0 : -range[1];
                for (ssize_t y = minY; y <= range[1]; ++y) {
                    v2 = v1 + ScalarType(y) * a2.asVector();
                    func(v2);
                }
            }
        }
        else if constexpr (Dim == 3) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);
            auto a3 = lattice.row(2);

            VectorType v1, v2, v3;
            for (ssize_t x = 0; x <= range[0]; ++x) {
                v1 = ScalarType(x) * a1.asVector();
                const ssize_t minY = x == 0 ? 0 : -range[1];
                for (ssize_t y = minY; y <= range[1]; ++y) {
                    v2 = v1 + ScalarType(y) * a2.asVector();
                    const ssize_t minZ = (x == 0 && y == 0) ? 0 : -range[2];
                    for (ssize_t z = minZ; z <= range[2]; ++z) {
                        v3 = v2 + ScalarType(z) * a3.asVector();
                        func(v3);
                    }
                }
            }
        }
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toDirect(const InvLatticeMatrix& invLattice) {
        if (type == Type::Cartesian) {
            toDirect(pos, invLattice);
            type = Type::Direct;
        }
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toDirect(PositionMatrix& target, const InvLatticeMatrix& invLattice) {
        target *= invLattice;
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::normalize_direct() {
        for (auto& elem : pos) {
            elem -= floor(elem);
            assert(ScalarType::Zero() <= elem && elem <= ScalarType::One());
        }
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::scale_direct(ScalarType factor) {
        assert(factor.isPositive());
        assert(type == Type::Direct);
        lattice *= factor;
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::scale_cartesian(ScalarType factor) {
        assert(factor.isPositive());
        assert(type == Type::Cartesian);
        lattice *= factor;
        pos *= factor;
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::unitToSuper_direct(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z) {
        assert(Dim == 3);
        assert(x > 0 && y > 0 && z > 0);
        const size_t numParticle = target.getRow();
        const size_t newNumParticle = x * y * z * target.getRow();
        PositionMatrix new_pos(newNumParticle, Dim);
        size_t index = 0;
        for (size_t i = 0; i < numParticle; ++i) {
            for (unsigned int x_ = 0; x_ < x; ++x_) {
                for (unsigned int y_ = 0; y_ < y; ++y_) {
                    for (unsigned int z_ = 0; z_ < z; ++z_) {
                        new_pos(index, 0) = (target(i, 0) + ScalarType(x_)) / ScalarType(x);
                        new_pos(index, 1) = (target(i, 1) + ScalarType(y_)) / ScalarType(y);
                        new_pos(index, 2) = (target(i, 2) + ScalarType(z_)) / ScalarType(z);
                        ++index;
                    }
                }
            }
        }
        target.swap(new_pos);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::superToUnit_direct(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z) {
        assert(Dim == 3);
        assert(x > 0 && y > 0 && z > 0);
        assert(target.getRow() % (x * y * z) == 0);

        auto colX = target.col(0);
        colX *= ScalarType(x);
        auto colY = target.col(1);
        colY *= ScalarType(y);
        auto colZ = target.col(2);
        colZ *= ScalarType(z);

        const size_t numParticle = target.getRow();
        const size_t newNumParticle = numParticle / (x * y * z);
        PositionMatrix new_pos(newNumParticle, Dim);
        size_t toFill = 0;
        size_t toCheck = 0;
        const ScalarType one = ScalarType::One();
        for (; toFill < newNumParticle; ++toFill) {
            for (; toCheck < numParticle; ++toCheck) {
                auto rowToCheck = target.row(toCheck);
                if (rowToCheck[0] <= one && rowToCheck[1] <= one && rowToCheck[2] <= one) {
                    auto rowToFill = new_pos.row(toFill);
                    rowToFill = rowToCheck.asVector();
                    ++toCheck;
                    break;
                }
            }
        }
        target.swap(new_pos);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::unitToSuper_direct(unsigned int x, unsigned int y, unsigned int z) {
        unitToSuper_direct(pos, x, y, z);
        auto rowX = lattice.row(0);
        rowX *= ScalarType(x);
        auto rowY = lattice.row(1);
        rowY *= ScalarType(y);
        auto rowZ = lattice.row(2);
        rowZ *= ScalarType(z);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::superToUnit_direct(unsigned int x, unsigned int y, unsigned int z) {
        superToUnit_direct(pos, x, y, z);
        const ScalarType inv_x = Core::reciprocal(ScalarType(x));
        const ScalarType inv_y = Core::reciprocal(ScalarType(y));
        const ScalarType inv_z = Core::reciprocal(ScalarType(z));

        auto rowX = lattice.row(0);
        rowX *= inv_x;
        auto rowY = lattice.row(1);
        rowY *= inv_y;
        auto rowZ = lattice.row(2);
        rowZ *= inv_z;
    }
}
