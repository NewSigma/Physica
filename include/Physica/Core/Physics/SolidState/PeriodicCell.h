/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Exception/BadConvergenceException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"

namespace Physica {
    namespace Internal {
        class PeriodicCellImpl {
        public:
            enum class Type : bool {
                Direct,
                Cartesian
            };
        };
    }

    enum class ExtendCellOption : char {
        AtomMajor,
        CellMajor
    };

    template<Scalar T, unsigned int Dim>
    class PeriodicCell : public Internal::PeriodicCellImpl {
        static_assert(Dim == 1 || Dim == 2 || Dim == 3, "[Error]: Unsupported dimention");

        using This = PeriodicCell<T, Dim>;
        using Base = Internal::PeriodicCellImpl;
        using Tv = T::ValueType;
        static_assert(!T::isComplex(), "[Error]: Complex position is not allowed");
    public:
        using LatticeMatrix = DenseMatrix<T, MatrixMajor::Row, Dim, Dim>;
        using InvLatticeMatrix = DenseMatrix<T, MatrixMajor::Col, Dim, Dim>;
        using PositionMatrix = DenseMatrix<T, MatrixMajor::Row, Dynamic, Dim>;
        using SearchRangeType = Array<ssize_t, Dim>;
    protected:
        using VectorType = DenseVector<T, Dim>;

        LatticeMatrix lattice;
        PositionMatrix pos;
        Type type;
    public:
        PeriodicCell();
        PeriodicCell(size_t numParticle, Type type);
        PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_, Type type_);
        template<Scalar U>
        PeriodicCell(const PeriodicCell<U, Dim>& cell);
        PeriodicCell(const This&) = default;
        PeriodicCell(This&&) noexcept = default;
        ~PeriodicCell() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] VectorType minDistVector(size_t id_from, size_t id_to) const;
        [[nodiscard]] VectorType minDistVector(VectorType from, size_t id_to) const;
        void normalize();
        void standardizeLattice();
        void scale(const T& factor);
        void niggliReduce(double precision, unsigned int maxIteration);
        void niggliReduce2D(unsigned int maxIteration);
        [[nodiscard]] LatticeMatrix makeRepLattice() const;

        void toDirect() { toDirect(makeInvLattice()); }
        void toCartesian();
        void toCartesian(PositionMatrix& target) const { toCartesian(target, lattice); }
        template<ExtendCellOption Option>
        void toSuperCell(unsigned int x, unsigned int y, unsigned int z);
        void toUnitCell(unsigned int x, unsigned int y, unsigned int z);
        template<ExtendCellOption Option>
        [[nodiscard]] PeriodicCell makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const;
        [[nodiscard]] PeriodicCell makeUnitCell(unsigned int x, unsigned int y, unsigned int z) const;

        H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static unsigned int getDim() { return Dim; }
        [[nodiscard]] auto&& getLattice(this auto&& self) noexcept { return self.lattice; }
        [[nodiscard]] auto&& getPos(this auto&& self) noexcept { return self.pos; }
        [[nodiscard]] Type getType() const noexcept { return type; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return pos.getRow(); }
        [[nodiscard]] T getVolume() const noexcept { return getVolume(lattice); }
        /* Setters */
        void setLattice(LatticeMatrix new_lattice) { lattice = new_lattice; }
        void setPos(PositionMatrix new_pos) noexcept { pos = std::move(new_pos); }
        void swapPos(PositionMatrix& __restrict new_pos) noexcept { pos.swap(new_pos); }
        /* Static members */
        [[nodiscard]] static LatticeMatrix makeLattice(
                T normA,
                T normB,
                T normC,
                T alpha,
                T beta,
                T gamma);
        [[nodiscard]] static LatticeMatrix makeLattice2D(
                T normA,
                T normB,
                T normC,
                T gamma);
        [[nodiscard]] static LatticeMatrix makeNiggliLattice(const LatticeMatrix& lattice, double precision, unsigned int maxIteration);
        [[nodiscard]] static LatticeMatrix makeNiggliLattice2D(const LatticeMatrix& lattice, unsigned int maxIteration);
        [[nodiscard]] static auto makeRepLattice(const LatticeMatrix& lattice);
        [[nodiscard]] static CoDiff<T> getVolume(const LatticeMatrix& lattice);
        static void toDirect(PositionMatrix& target, const LatticeMatrix& lattice);
        static void toCartesian(PositionMatrix& target, const LatticeMatrix& lattice);
        [[nodiscard]] static SearchRangeType estimateRange(const LatticeMatrix& lattice, Tv cutoff);
        static void forCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, std::invocable<VectorType> auto fn);
        static void forReducedCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, std::invocable<VectorType> auto fn);

        template<bool IsUnitLattice>
        static void forLatticeCloud(std::invocable<Vector3D<T>, Index3D> auto fn, const LatticeMatrix& lattice, Index3D dim);
    protected:
        [[nodiscard]] InvLatticeMatrix makeInvLattice() const noexcept { return lattice.inv(); }
        void toDirect(const InvLatticeMatrix& invLattice);
        void normalize_direct();
        void scale_direct(const T& factor);
        void scale_cartesian(const T& factor);
        template<ExtendCellOption Option>
        void toSuperCellDirect(unsigned int x, unsigned int y, unsigned int z);
        void toUnitCellDirect(unsigned int x, unsigned int y, unsigned int z);
        /* Static members */
        static void toDirect(PositionMatrix& target, const InvLatticeMatrix& invLattice);
        template<ExtendCellOption Option>
        static void toSuperPosDirect(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z);
        static void toUnitPosDirect(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z);

        friend class device_obj<This>;
    };

    template<Scalar T, unsigned int Dim>
    PeriodicCell<T, Dim>::PeriodicCell()
            : lattice(LatticeMatrix::identity(Dim))
            , pos()
            , type(Type::Direct) {}

    template<Scalar T, unsigned int Dim>
    PeriodicCell<T, Dim>::PeriodicCell(size_t numParticle, Type type)
            : lattice(LatticeMatrix::identity(Dim))
            , pos(numParticle, Dim)
            , type(type) {}

    template<Scalar T, unsigned int Dim>
    PeriodicCell<T, Dim>::PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_, Type type_)
            : lattice(std::move(lattice_))
            , pos(std::move(pos_))
            , type(type_) {}

    template<Scalar T, unsigned int Dim>
    template<Scalar U>
    PeriodicCell<T, Dim>::PeriodicCell(const PeriodicCell<U, Dim>& cell)
            : lattice(cell.getLattice())
            , pos(cell.getPos())
            , type(cell.getType()) {}

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::minDistVector(size_t id_from, size_t id_to) const -> VectorType {
        T record_dist = std::numeric_limits<T>::max();
        VectorType result{};
        VectorType delta = pos.row(id_to) - pos.row(id_from);
        if (type == Type::Direct)
            delta = lattice.transpose() * delta;
        if constexpr (Dim == 1) {
            VectorType v1;
            for (int x = -1; x <= 1; ++x) {
                const bool isSelf = id_from == id_to && x == 0;
                if (isSelf)
                    continue; // We are not interested to distance between particle and itself
                v1 = delta + lattice.row(0) * T(x);
                const T squared_norm = v1.squaredNorm();
                if (squared_norm < record_dist) {
                    record_dist = squared_norm;
                    result = v1;
                }
            }
        }
        else if constexpr (Dim == 2) {
            VectorType v1, v2;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0) * T(x);
                for (int y = -1; y <= 1; ++y) {
                    const bool isSelf = id_from == id_to && x == 0 && y == 0;
                    if (isSelf)
                        continue; // We are not interested to distance between particle and itself
                    v2 = v1 + lattice.row(1) * T(y);
                    const T squared_norm = v2.squaredNorm();
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
                v1 = delta + lattice.row(0) * T(x);
                for (int y = -1; y <= 1; ++y) {
                    v2 = v1 + lattice.row(1) * T(y);
                    for (int z = -1; z <= 1; ++z) {
                        const bool isSelf = id_from == id_to && x == 0 && y == 0 && z == 0;
                        if (isSelf)
                            continue; // We are not interested to distance between particle and itself
                        v3 = v2 + lattice.row(2) * T(z);
                        const T squared_norm = v3.squaredNorm();
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

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::minDistVector(VectorType from, size_t id_to) const -> VectorType {
        T record_dist = std::numeric_limits<T>::max();
        VectorType result{};
        VectorType delta;
        if (type == Type::Direct)
            delta = lattice.transpose() * pos.row(id_to);
        else
            delta = pos.row(id_to);
        delta -= from;

        if constexpr (Dim == 1) {
            VectorType v1;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0) * T(x);
                const T squared_norm = v1.squaredNorm();
                if (squared_norm < record_dist) {
                    record_dist = squared_norm;
                    result = v1;
                }
            }
        }
        else if constexpr (Dim == 2) {
            VectorType v1, v2;
            for (int x = -1; x <= 1; ++x) {
                v1 = delta + lattice.row(0) * T(x);
                for (int y = -1; y <= 1; ++y) {
                    v2 = v1 + lattice.row(1) * T(y);
                    const T squared_norm = v2.squaredNorm();
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
                v1 = delta + lattice.row(0) * T(x);
                for (int y = -1; y <= 1; ++y) {
                    v2 = v1 + lattice.row(1) * T(y);
                    for (int z = -1; z <= 1; ++z) {
                        v3 = v2 + lattice.row(2) * T(z);
                        const T squared_norm = v3.squaredNorm();
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

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::normalize() {
        if (type == Type::Cartesian) {
            toDirect(makeInvLattice());
            normalize_direct();
            toCartesian();
        }
        else
            normalize_direct();
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::standardizeLattice() {
        static_assert(Dim == 3, "[Error]: standardizeLattice only implemented for 3D");
        using MatrixType = LatticeMatrix::ColMatrix;
        MatrixType temp = getLattice().transpose();
        Vector3D<T> buffer{};
        temp.col(0).householder(buffer);
        applyHouseholder(buffer, temp);

        auto buffer1 = buffer.template head<2>();
        temp.col(1).tail(1).householder(buffer1);
        auto corner = temp.bottomRightCorner(1);
        applyHouseholder(buffer1, corner);
        for (int i = 0; i < 3; ++i) {
            if (temp[i, i].isNegative()) {
                auto row = temp.row(i);
                row = -row;
            }
        }
        temp[1, 0] = temp[2, 0] = temp[2, 1] = T(0);
        getLattice() = temp.transpose();
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::scale(const T& factor) {
        if (type == Type::Direct)
            scale_direct(factor);
        else
            scale_cartesian(factor);
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::niggliReduce(double precision, unsigned int maxIteration) {
        if (type == Type::Direct)
            toCartesian();
        setLattice(makeNiggliLattice(lattice, precision, maxIteration));
        toDirect();
        normalize_direct();
        if (type == Type::Cartesian)
            toCartesian();
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::niggliReduce2D(unsigned int maxIteration) {
        if (type == Type::Direct)
            toCartesian();
        setLattice(makeNiggliLattice2D(lattice, maxIteration));
        toDirect();
        normalize_direct();
        if (type == Type::Cartesian)
            toCartesian();
    }

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::makeRepLattice() const -> LatticeMatrix {
        return makeRepLattice(lattice);
    }

    template<Scalar T, unsigned int Dim>
    template<ExtendCellOption Option>
    void PeriodicCell<T, Dim>::toSuperCell(unsigned int x, unsigned int y, unsigned int z) {
        if (type == Type::Cartesian) {
            toDirect();
            toSuperCellDirect<Option>(x, y, z);
            toCartesian();
        }
        else
            toSuperCellDirect<Option>(x, y, z);
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toUnitCell(unsigned int x, unsigned int y, unsigned int z) {
        if (type == Type::Cartesian) {
            toDirect();
            toUnitCellDirect(x, y, z);
            toCartesian();
        }
        else
            toUnitCellDirect(x, y, z);
    }

    template<Scalar T, unsigned int Dim>
    template<ExtendCellOption Option>
    auto PeriodicCell<T, Dim>::makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const -> This {
        This result = *this;
        result.toSuperCell<Option>(x, y, z);
        return result;
    }

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::makeUnitCell(unsigned int x, unsigned int y, unsigned int z) const -> This {
        This result = *this;
        result.toUnitCell(x, y, z);
        return result;
    }
#ifdef PHYSICA_HDF5
    template<Scalar T, unsigned int Dim>
    H5Group PeriodicCell<T, Dim>::read(const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        lattice.read(group, "lattice");
        const auto posDataset = pos.read(group, "pos");
        const auto typeAttr = posDataset.openAttribute("Type");
        typeAttr.read(H5::PredType::NATIVE_INT8, &type);
        return H5Group(group);
    }

    template<Scalar T, unsigned int Dim>
    H5Group PeriodicCell<T, Dim>::write(H5Loc& loc, const char* name) const {
        auto group = loc.openGroup(name);
        lattice.write(group, "lattice");
        auto posDataset = pos.write(group, "pos");
        auto typeAttr = posDataset.createAttribute("Type", H5::PredType::NATIVE_INT8, H5DataSpace<1>(1));
        typeAttr.write(H5::PredType::NATIVE_INT8, &type);
        return H5Group(group);
    }
#endif
    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lattice.swap(obj.lattice);
        pos.swap(obj.pos);
        std::swap(type, obj.type);
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toCartesian() {
        if (type == Type::Direct) {
            pos *= lattice;
            type = Type::Cartesian;
        }
    }

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::makeLattice(T normA, T normB, T normC, T alpha, T beta, T gamma) -> LatticeMatrix {
        LatticeMatrix result{};
        result[0, 0] = normA;
        result[0, 1] = T(0);
        result[0, 2] = T(0);
        result[1, 0] = normB * cos(gamma);
        result[1, 1] = normB * sin(gamma);
        result[1, 2] = T(0);
        result[2, 0] = normC * cos(beta);
        result[2, 1] = normC * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
        result[2, 2] = sqrt(square(normC) - square(result[2, 0]) - square(result[2, 1]));
        return result;
    }

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::makeLattice2D(T normA, T normB, T normC, T gamma) -> LatticeMatrix {
        LatticeMatrix result{};
        result[0, 0] = normA;
        result[0, 1] = T(0);
        result[0, 2] = T(0);
        result[1, 0] = normB * cos(gamma);
        result[1, 1] = normB * sin(gamma);
        result[1, 2] = T(0);
        result[2, 0] = T(0);
        result[2, 1] = T(0);
        result[2, 2] = normC;
        return result;
    }
    /**
     * Step 1-8 referenced from [1]
     *
     * Reference:
     * [1] Acta Cryst. A32, 297 (1976); https://doi.org/10.1107/S0567739476000636
     */
    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::makeNiggliLattice(const LatticeMatrix& lattice, double precision, unsigned int maxIteration) -> LatticeMatrix{
        assert(precision > 0 && "[Error]: Precision should be positive");
        assert(maxIteration > 0 && "[Error]: Set maxIteration = 0 does nothing");
        T squaredNormA = lattice.row(0).squaredNorm();
        T squaredNormB = lattice.row(1).squaredNorm();
        T squaredNormC = lattice.row(2).squaredNorm();
        T dot1 = T(2) * (lattice.row(1) * lattice.row(2));
        T dot2 = T(2) * (lattice.row(0) * lattice.row(2));
        T dot3 = T(2) * (lattice.row(0) * lattice.row(1));
        unsigned int iteration = 0;
        while (true) {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Cannot reduce to niggli cell within required iterations");
            iteration += 1;
            /* Step 1 */ {
                const bool changeAB1 = squaredNormA > squaredNormB;
                const bool changeAB2 = scalarNear(squaredNormA, squaredNormB, precision) && (abs(dot1) > abs(dot2));
                if (changeAB1 || changeAB2) {
                    squaredNormA.swap(squaredNormB);
                    dot1.swap(dot2);
                }
            }
            /* Step 2 */ {
                const bool changeBC1 = squaredNormB > squaredNormC;
                const bool changeBC2 = scalarNear(squaredNormB, squaredNormC, precision) && (abs(dot2) > abs(dot3));
                if (changeBC1 || changeBC2) {
                    squaredNormB.swap(squaredNormA);
                    dot2.swap(dot3);
                    continue;
                }
            }
            /* Step 3 and 4 */ {
                const T temp = dot1 * dot2 * dot3;
                if (temp.isPositive()) {
                    dot1 = abs(dot1);
                    dot2 = abs(dot2);
                    dot3 = abs(dot3);
                }
                else {
                    dot1 = -abs(dot1);
                    dot2 = -abs(dot2);
                    dot3 = -abs(dot3);
                }
            }
            /* Step 5 */ {
                const bool cond1 = abs(dot1) > squaredNormB;
                const bool cond2 = scalarNear(dot1, squaredNormB, precision) && (T(2) * dot2 < dot3);
                const bool cond3 = scalarNear(dot1, -squaredNormB, precision) && dot3.isNegative();
                if (cond1 || cond2 || cond3) {
                    assert(!dot1.isZero() && "[Error]: Condition 1~3 shall ensure dot1 is not zero");
                    const T sign = unit(dot1);
                    squaredNormC += squaredNormB - abs(dot1);
                    dot2 -= dot3 * sign;
                    dot1 -= T(2) * squaredNormB * sign;
                    continue;
                }
            }
            /* Step 6 */ {
                const bool cond1 = abs(dot2) > squaredNormA;
                const bool cond2 = scalarNear(dot2, squaredNormA, precision) && (T(2) * dot1 < dot3);
                const bool cond3 = scalarNear(dot2, -squaredNormA, precision) && dot3.isNegative();
                if (cond1 || cond2 || cond3) {
                    assert(!dot2.isZero() && "[Error]: Condition 1~3 shall ensure dot2 is not zero");
                    const T sign = unit(dot2);
                    squaredNormC += squaredNormA - abs(dot2);
                    dot1 -= dot3 * sign;
                    dot2 -= T(2) * squaredNormA * sign;
                    continue;
                }
            }
            /* Step 7 */ {
                const bool cond1 = abs(dot3) > squaredNormA;
                const bool cond2 = scalarNear(dot3, squaredNormA, precision) && (T(2) * dot1 < dot2);
                const bool cond3 = scalarNear(dot3, -squaredNormA, precision) && dot2.isNegative();
                if (cond1 || cond2 || cond3) {
                    assert(!dot3.isZero() && "[Error]: Condition 1~3 shall ensure dot3 is not zero");
                    const T sign = unit(dot3);
                    squaredNormB += squaredNormA - abs(dot3);
                    dot1 -= dot2 * sign;
                    dot3 -= T(2) * squaredNormA * sign;
                    continue;
                }
            }
            /* Step 8 */ {
                const T temp = dot1 + dot2 + dot3 + squaredNormA + squaredNormB;
                const bool cond1 = temp.isNegative();
                const bool cond2 = scalarNear(temp, T(0), precision) && (T(2) * (squaredNormA + dot2) + dot3).isPositive();
                if (cond1 || cond2) {
                    squaredNormC += temp;
                    dot1 += T(2) * squaredNormB + dot3;
                    dot2 += T(2) * squaredNormA + dot3;
                    continue;
                }
            }
            break;
        }
        const T normA = sqrt(squaredNormA);
        const T normB = sqrt(squaredNormB);
        const T normC = sqrt(squaredNormC);
        const T alpha = arccos(dot1 / (T(2) * normB * normC));
        const T beta = arccos(dot2 / (T(2) * normA * normC));
        const T gamma = arccos(dot3 / (T(2) * normA * normB));
        return makeLattice(normA, normB, normC, alpha, beta, gamma);
    }
    /**
     * A directly simplified version of above function
     */
    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::makeNiggliLattice2D(const LatticeMatrix& lattice, unsigned int maxIteration) -> LatticeMatrix {
        assert(maxIteration > 0 && "[Error]: Set maxIteration = 0 does nothing");
        T squaredNormA = lattice.row(0).squaredNorm();
        T squaredNormB = lattice.row(1).squaredNorm();
        T dot = T(2) * (lattice.row(0) * lattice.row(1));
        unsigned int iteration = 0;
        while (true) {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Cannot reduce to niggli cell within required iterations");
            iteration += 1;

            if (squaredNormA > squaredNormB)
                squaredNormA.swap(squaredNormB);
            dot = -abs(dot);

            if (abs(dot) > squaredNormA) {
                assert(!dot.isZero() && "[Error]: Condition 1~3 shall ensure dot is not zero");
                squaredNormB += squaredNormA - abs(dot);
                dot -= T(2) * squaredNormA * unit(dot);
                continue;
            }
            break;
        }
        const T normA = sqrt(squaredNormA);
        const T normB = sqrt(squaredNormB);
        const T gamma = arccos(dot / (T(2) * normA * normB));
        return makeLattice2D(normA, normB, lattice[2, 2], gamma);
    }

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::makeRepLattice(const LatticeMatrix& lattice) {
        using RtnTy = LatticeMatrix::template rebind_scalar<Tv>;
        const auto& latt = lattice.values();
        RtnTy result{};
        result.row(0) = cross(latt.row(1), latt.row(2));
        result.row(1) = cross(latt.row(2), latt.row(0));
        result.row(2) = cross(latt.row(0), latt.row(1));
        const auto factor = Tv(2 * M_PI) / (latt.row(0) * result.row(0));
        result *= factor;
        return result;
    }

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::getVolume(const LatticeMatrix& lattice) -> CoDiff<T> {
        return abs(lattice.det());
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toDirect(PositionMatrix& target, const LatticeMatrix& lattice) {
        const InvLatticeMatrix inv = lattice.inv();
        toDirect(target, inv);
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toCartesian(PositionMatrix& target, const LatticeMatrix& lattice) {
        target *= lattice;
    }

    template<Scalar T, unsigned int Dim>
    auto PeriodicCell<T, Dim>::estimateRange(const LatticeMatrix& lattice, Tv cutoff) -> SearchRangeType {
        const auto repLatt = makeRepLattice(lattice);
        const Tv factor = cutoff * Tv(1 / (2 * M_PI));
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

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::forCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, std::invocable<VectorType> auto fn) {
        if constexpr (Dim == 1) {
            auto a1 = lattice.row(0);

            VectorType v1;
            for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                v1 = T(x) * a1;
                fn(v1);
            }
        }
        else if constexpr (Dim == 2) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);

            VectorType v1, v2;
            for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                v1 = T(x) * a1;
                for (ssize_t y = -range[1]; y <= range[1]; ++y) {
                    v2 = v1 + T(y) * a2;
                    fn(v2);
                }
            }
        }
        else if constexpr (Dim == 3) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);
            auto a3 = lattice.row(2);

            VectorType v1, v2, v3;
            for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                v1 = T(x) * a1;
                for (ssize_t y = -range[1]; y <= range[1]; ++y) {
                    v2 = v1 + T(y) * a2;
                    for (ssize_t z = -range[2]; z <= range[2]; ++z) {
                        v3 = v2 + T(z) * a3;
                        fn(v3);
                    }
                }
            }
        }
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::forReducedCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, std::invocable<VectorType> auto fn) {
        if constexpr (Dim == 1) {
            auto a1 = lattice.row(0);

            VectorType v1;
            for (ssize_t x = 0; x <= range[0]; ++x) {
                v1 = T(x) * a1;
                fn(v1);
            }
        }
        else if constexpr (Dim == 2) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);

            VectorType v1, v2;
            for (ssize_t x = 0; x <= range[0]; ++x) {
                v1 = T(x) * a1;
                const ssize_t minY = x == 0 ? 0 : -range[1];
                for (ssize_t y = minY; y <= range[1]; ++y) {
                    v2 = v1 + T(y) * a2;
                    fn(v2);
                }
            }
        }
        else if constexpr (Dim == 3) {
            auto a1 = lattice.row(0);
            auto a2 = lattice.row(1);
            auto a3 = lattice.row(2);

            VectorType v1, v2, v3;
            for (ssize_t x = 0; x <= range[0]; ++x) {
                v1 = T(x) * a1;
                const ssize_t minY = x == 0 ? 0 : -range[1];
                for (ssize_t y = minY; y <= range[1]; ++y) {
                    v2 = v1 + T(y) * a2;
                    const ssize_t minZ = (x == 0 && y == 0) ? 0 : -range[2];
                    for (ssize_t z = minZ; z <= range[2]; ++z) {
                        v3 = v2 + T(z) * a3;
                        fn(v3);
                    }
                }
            }
        }
    }

    template<Scalar T, unsigned int Dim>
    template<bool IsUnitLattice>
    void PeriodicCell<T, Dim>::forLatticeCloud(std::invocable<Vector3D<T>, Index3D> auto fn, const LatticeMatrix& lattice, Index3D dim) {
        static_assert(Dim == 3, "[Error]: Not implemented");
        LatticeMatrix sub_lattice{};
        if constexpr (IsUnitLattice)
            sub_lattice = lattice;
        else {
            for (int i = 0; i < Dim; ++i)
                sub_lattice.row(i) = lattice.row(i) * reciprocal(T(dim[i]));
        }

        Vector3D<T> v1, v2, v3;
        auto a1 = sub_lattice.row(0);
        auto a2 = sub_lattice.row(1);
        auto a3 = sub_lattice.row(2);
        for (size_t x = 0; x < dim[0]; ++x) {
            v1 = T(x) * a1;
            for (size_t y = 0; y < dim[1]; ++y) {
                v2 = v1 + T(y) * a2;
                for (size_t z = 0; z < dim[2]; ++z) {
                    v3 = v2 + T(z) * a3;
                    fn(v3, Index3D{x, y, z});
                }
            }
        }
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toDirect(const InvLatticeMatrix& invLattice) {
        if (type == Type::Cartesian) {
            toDirect(pos, invLattice);
            type = Type::Direct;
        }
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::normalize_direct() {
        for (auto& elem : pos.asArray()) {
            elem -= floor(elem);
            assert(T(0) <= elem && elem <= T(1));
        }
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::scale_direct(const T& factor) {
        assert(factor.isPositive());
        assert(type == Type::Direct);
        lattice *= factor;
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::scale_cartesian(const T& factor) {
        assert(factor.isPositive());
        assert(type == Type::Cartesian);
        lattice *= factor;
        pos *= factor;
    }

    template<Scalar T, unsigned int Dim>
    template<ExtendCellOption Option>
    void PeriodicCell<T, Dim>::toSuperCellDirect(unsigned int x, unsigned int y, unsigned int z) {
        assert(type == Type::Direct);
        toSuperPosDirect<Option>(pos, x, y, z);
        auto rowX = lattice.row(0);
        rowX *= T(x);
        auto rowY = lattice.row(1);
        rowY *= T(y);
        auto rowZ = lattice.row(2);
        rowZ *= T(z);
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toUnitCellDirect(unsigned int x, unsigned int y, unsigned int z) {
        assert(type == Type::Direct);
        toUnitPosDirect(pos, x, y, z);
        const T inv_x = reciprocal(T(x));
        const T inv_y = reciprocal(T(y));
        const T inv_z = reciprocal(T(z));

        auto rowX = lattice.row(0);
        rowX *= inv_x;
        auto rowY = lattice.row(1);
        rowY *= inv_y;
        auto rowZ = lattice.row(2);
        rowZ *= inv_z;
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toDirect(PositionMatrix& target, const InvLatticeMatrix& invLattice) {
        target *= invLattice;
    }

    template<Scalar T, unsigned int Dim>
    template<ExtendCellOption Option>
    void PeriodicCell<T, Dim>::toSuperPosDirect(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z) {
        assert(Dim == 3);
        assert(x > 0 && y > 0 && z > 0);
        const size_t numParticle = target.getRow();
        const size_t newNumParticle = x * y * z * numParticle;
        PositionMatrix new_pos(newNumParticle, Dim);
        size_t index = 0;
        if constexpr (Option == ExtendCellOption::AtomMajor) {
            for (size_t i = 0; i < numParticle; ++i) {
                for (unsigned int x_ = 0; x_ < x; ++x_) {
                    for (unsigned int y_ = 0; y_ < y; ++y_) {
                        for (unsigned int z_ = 0; z_ < z; ++z_) {
                            new_pos[index, 0] = (target[i, 0] + T(x_)) / T(x);
                            new_pos[index, 1] = (target[i, 1] + T(y_)) / T(y);
                            new_pos[index, 2] = (target[i, 2] + T(z_)) / T(z);
                            ++index;
                        }
                    }
                }
            }
        }
        else if constexpr (Option == ExtendCellOption::CellMajor) {
            for (unsigned int x_ = 0; x_ < x; ++x_) {
                for (unsigned int y_ = 0; y_ < y; ++y_) {
                    for (unsigned int z_ = 0; z_ < z; ++z_) {
                        for (size_t i = 0; i < numParticle; ++i) {
                            new_pos[index, 0] = (target[i, 0] + T(x_)) / T(x);
                            new_pos[index, 1] = (target[i, 1] + T(y_)) / T(y);
                            new_pos[index, 2] = (target[i, 2] + T(z_)) / T(z);
                            ++index;
                        }
                    }
                }
            }
        }
        else
            static_assert(Option == ExtendCellOption::AtomMajor, "[Error]: Invalid option");
        target.swap(new_pos);
    }

    template<Scalar T, unsigned int Dim>
    void PeriodicCell<T, Dim>::toUnitPosDirect(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z) {
        assert(Dim == 3);
        assert(x > 0 && y > 0 && z > 0);
        assert(target.getRow() % (x * y * z) == 0);

        auto colX = target.col(0);
        colX *= T(x);
        auto colY = target.col(1);
        colY *= T(y);
        auto colZ = target.col(2);
        colZ *= T(z);

        const size_t numParticle = target.getRow();
        const size_t newNumParticle = numParticle / (x * y * z);
        PositionMatrix new_pos(newNumParticle, Dim);
        size_t toFill = 0;
        size_t toCheck = 0;
        const T one = T(1);
        for (; toFill < newNumParticle; ++toFill) {
            for (; toCheck < numParticle; ++toCheck) {
                auto rowToCheck = target.row(toCheck);
                if (rowToCheck[0] <= one && rowToCheck[1] <= one && rowToCheck[2] <= one) {
                    auto rowToFill = new_pos.row(toFill);
                    rowToFill = rowToCheck;
                    ++toCheck;
                    break;
                }
            }
        }
        target.swap(new_pos);
    }
}
