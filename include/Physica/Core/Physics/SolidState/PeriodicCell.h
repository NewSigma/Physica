/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Exception/BadConvergenceException.h>
#include <Physica/Core/MultiPrecision/Real.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>

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

    enum class ExtendCellOption {
        AtomMajor,
        CellMajor
    };

    template<class ScalarType, unsigned int Dim>
    class PeriodicCell : public Internal::PeriodicCellImpl {
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid ScalarType");
        static_assert(Dim == 1 || Dim == 2 || Dim == 3, "[Error]: Unsupported dimention");

        using This = PeriodicCell<ScalarType, Dim>;
        using Base = Internal::PeriodicCellImpl;
        using ValueType = typename ScalarType::ValueType;
    public:
        using LatticeMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dim, Dim>;
        using InvLatticeMatrix = DenseMatrix<ScalarType, MatrixOption::Col | MatrixOption::Element, Dim, Dim>;
        using PositionMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dynamic, Dim>;
        using MomentumMatrix = PositionMatrix;
        using SearchRangeType = Array<ssize_t, Dim>;
    protected:
        using VectorType = Vector<ScalarType, Dim>;

        LatticeMatrix lattice;
        PositionMatrix pos;
        Type type;
    public:
        PeriodicCell();
        PeriodicCell(size_t numParticle, Type type_);
        PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_, Type type_);
        template<class OtherScalar>  PeriodicCell(const PeriodicCell<OtherScalar, Dim>& cell);
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
        void niggliReduce(double precision, unsigned int maxIteration);
        void niggliReduce2D(unsigned int maxIteration);
        [[nodiscard]] inline LatticeMatrix makeRepLattice() const;

        void toDirect() { toDirect(makeInvLattice()); }
        void toCartesian();
        void toCartesian(PositionMatrix& target) const { toCartesian(target, lattice); }
        template<ExtendCellOption Option>
        void toSuperCell(unsigned int x, unsigned int y, unsigned int z);
        void toUnitCell(unsigned int x, unsigned int y, unsigned int z);
        template<ExtendCellOption Option>
        [[nodiscard]] PeriodicCell makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const;
        [[nodiscard]] PeriodicCell makeUnitCell(unsigned int x, unsigned int y, unsigned int z) const;

        H5Group read(const H5Location& loc, const char* name);
        H5Group write(H5Location& loc, const char* name) const;
        void swap(PeriodicCell& __restrict cell) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static unsigned int getDim() { return Dim; }
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] const PositionMatrix& getPos() const noexcept { return pos; }
        [[nodiscard]] Type getType() const noexcept { return type; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return pos.getRow(); }
        [[nodiscard]] ScalarType getVolume() const noexcept { return getVolume(lattice); }
        /* Setters */
        void setLattice(LatticeMatrix new_lattice) { lattice = new_lattice; }
        void setPos(PositionMatrix new_pos) noexcept { pos = std::move(new_pos); }
        void swapPos(PositionMatrix& __restrict new_pos) noexcept { pos.swap(new_pos); }
        /* Static members */
        [[nodiscard]] static LatticeMatrix makeLattice(
                ScalarType normA,
                ScalarType normB,
                ScalarType normC,
                ScalarType alpha,
                ScalarType beta,
                ScalarType gamma);
        [[nodiscard]] static LatticeMatrix makeLattice2D(
                ScalarType normA,
                ScalarType normB,
                ScalarType normC,
                ScalarType gamma);
        [[nodiscard]] static LatticeMatrix makeNiggliLattice(const LatticeMatrix& lattice, double precision, unsigned int maxIteration);
        [[nodiscard]] static LatticeMatrix makeNiggliLattice2D(const LatticeMatrix& lattice, unsigned int maxIteration);
        [[nodiscard]] static LatticeMatrix makeRepLattice(const LatticeMatrix& lattice);
        [[nodiscard]] static ScalarType getVolume(const LatticeMatrix& lattice);
        static void toDirect(PositionMatrix& target, const LatticeMatrix& lattice);
        static void toCartesian(PositionMatrix& target, const LatticeMatrix& lattice);
        [[nodiscard]] static SearchRangeType estimateRange(const LatticeMatrix& cell, ValueType cutoff);
        template<class Functor> static void forCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func);
        template<class Functor> static void forReducedCellInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func);
    protected:
        [[nodiscard]] InvLatticeMatrix makeInvLattice() const noexcept { return lattice.inverse(); }
        void toDirect(const InvLatticeMatrix& invLattice);
        void normalize_direct();
        void scale_direct(ScalarType factor);
        void scale_cartesian(ScalarType factor);
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

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>::PeriodicCell()
            : lattice(LatticeMatrix::unitMatrix(Dim))
            , pos()
            , type(Type::Direct) {}

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>::PeriodicCell(size_t numParticle, Type type_)
            : PeriodicCell() {
        pos.resize(numParticle, Dim);
        type = type_;
    }

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>::PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_, Type type_)
            : lattice(std::move(lattice_))
            , pos(std::move(pos_))
            , type(type_) {}

    template<class ScalarType, unsigned int Dim>
    template<class OtherScalar>
    PeriodicCell<ScalarType, Dim>::PeriodicCell(const PeriodicCell<OtherScalar, Dim>& cell)
            : lattice(cell.getLattice())
            , pos(cell.getPos())
            , type(cell.getType()) {}

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
    void PeriodicCell<ScalarType, Dim>::niggliReduce(double precision, unsigned int maxIteration) {
        if (type == Type::Direct)
            toCartesian();
        setLattice(makeNiggliLattice(lattice, precision, maxIteration));
        toDirect();
        normalize_direct();
        if (type == Type::Cartesian)
            toCartesian();
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::niggliReduce2D(unsigned int maxIteration) {
        if (type == Type::Direct)
            toCartesian();
        setLattice(makeNiggliLattice2D(lattice, maxIteration));
        toDirect();
        normalize_direct();
        if (type == Type::Cartesian)
            toCartesian();
    }

    template<class ScalarType, unsigned int Dim>
    inline typename PeriodicCell<ScalarType, Dim>::LatticeMatrix PeriodicCell<ScalarType, Dim>::makeRepLattice() const {
        return makeRepLattice(lattice);
    }

    template<class ScalarType, unsigned int Dim>
    template<ExtendCellOption Option>
    void PeriodicCell<ScalarType, Dim>::toSuperCell(unsigned int x, unsigned int y, unsigned int z) {
        if (type == Type::Cartesian) {
            toDirect();
            toSuperCellDirect<Option>(x, y, z);
            toCartesian();
        }
        else
            toSuperCellDirect<Option>(x, y, z);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toUnitCell(unsigned int x, unsigned int y, unsigned int z) {
        if (type == Type::Cartesian) {
            toDirect();
            toUnitCellDirect(x, y, z);
            toCartesian();
        }
        else
            toUnitCellDirect(x, y, z);
    }

    template<class ScalarType, unsigned int Dim>
    template<ExtendCellOption Option>
    PeriodicCell<ScalarType, Dim> PeriodicCell<ScalarType, Dim>::makeSuperCell(unsigned int x, unsigned int y, unsigned int z) const {
        PeriodicCell<ScalarType, Dim> result = *this;
        result.toSuperCell<Option>(x, y, z);
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim> PeriodicCell<ScalarType, Dim>::makeUnitCell(unsigned int x, unsigned int y, unsigned int z) const {
        PeriodicCell<ScalarType, Dim> result = *this;
        result.toUnitCell(x, y, z);
        return result;
    }
#ifdef PHYSICA_HDF5
    template<class ScalarType, unsigned int Dim>
    H5Group PeriodicCell<ScalarType, Dim>::read(const H5Location& loc, const char* name) {
        const auto group = loc.openGroup(name);
        lattice.read(group, "lattice");
        const auto posDataset = pos.read(group, "pos");
        const auto typeAttr = posDataset.openAttribute("Type");
        typeAttr.read(H5::PredType::NATIVE_INT8, &type);
        return H5Group(group);
    }

    template<class ScalarType, unsigned int Dim>
    H5Group PeriodicCell<ScalarType, Dim>::write(H5Location& loc, const char* name) const {
        auto group = loc.openGroup(name);
        lattice.write(group, "lattice");
        auto posDataset = pos.write(group, "pos");
        auto typeAttr = posDataset.createAttribute("Type", H5::PredType::NATIVE_INT8, H5DataSpace<1>(1));
        typeAttr.write(H5::PredType::NATIVE_INT8, &type);
        return H5Group(group);
    }
#endif
    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::swap(PeriodicCell& __restrict cell) noexcept {
        assert(this != &cell && "[Error]: Self swap is likely a bug");
        lattice.swap(cell.lattice);
        pos.swap(cell.pos);
        std::swap(type, cell.type);
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toCartesian() {
        if (type == Type::Direct) {
            pos *= lattice;
            type = Type::Cartesian;
        }
    }

    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::LatticeMatrix PeriodicCell<ScalarType, Dim>::makeLattice(
            ScalarType normA,
            ScalarType normB,
            ScalarType normC,
            ScalarType alpha,
            ScalarType beta,
            ScalarType gamma) {
        LatticeMatrix result{};
		result(0, 0) = normA;
		result(0, 1) = ScalarType(0);
		result(0, 2) = ScalarType(0);
		result(1, 0) = normB * cos(gamma);
		result(1, 1) = normB * sin(gamma);
		result(1, 2) = ScalarType(0);
		result(2, 0) = normC * cos(beta);
		result(2, 1) = normC * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
		result(2, 2) = sqrt(square(normC) - square(result(2, 0)) - square(result(2, 1)));
		return result;
    }

    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::LatticeMatrix PeriodicCell<ScalarType, Dim>::makeLattice2D(
            ScalarType normA,
            ScalarType normB,
            ScalarType normC,
            ScalarType gamma) {
        LatticeMatrix result{};
		result(0, 0) = normA;
		result(0, 1) = ScalarType(0);
		result(0, 2) = ScalarType(0);
		result(1, 0) = normB * cos(gamma);
		result(1, 1) = normB * sin(gamma);
		result(1, 2) = ScalarType(0);
		result(2, 0) = ScalarType(0);
		result(2, 1) = ScalarType(0);
		result(2, 2) = normC;
		return result;
    }
    /**
     * Step 1-8 referenced from [1]
     * 
     * Reference:
     * [1] Acta Cryst. A32, 297 (1976); https://doi.org/10.1107/S0567739476000636
     */
    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::LatticeMatrix PeriodicCell<ScalarType, Dim>::makeNiggliLattice(
            const LatticeMatrix& lattice, double precision, unsigned int maxIteration) {
        assert(precision > 0 && "[Error]: Precision should be positive");
        assert(maxIteration > 0 && "[Error]: Set maxIteration = 0 does nothing");
        ScalarType squaredNormA = lattice.row(0).squaredNorm();
        ScalarType squaredNormB = lattice.row(1).squaredNorm();
        ScalarType squaredNormC = lattice.row(2).squaredNorm();
        ScalarType dot1 = ScalarType(2) * (lattice.row(1).asVector() * lattice.row(2).asVector());
        ScalarType dot2 = ScalarType(2) * (lattice.row(0).asVector() * lattice.row(2).asVector());
        ScalarType dot3 = ScalarType(2) * (lattice.row(0).asVector() * lattice.row(1).asVector());
        unsigned int iteration = 0;
        do {
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
                const ScalarType temp = dot1 * dot2 * dot3;
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
                const bool cond2 = scalarNear(dot1, squaredNormB, precision) && (ScalarType(2) * dot2 < dot3);
                const bool cond3 = scalarNear(dot1, -squaredNormB, precision) && dot3.isNegative();
                if (cond1 || cond2 || cond3) {
                    assert(!dot1.isZero() && "[Error]: Condition 1~3 shall ensure dot1 is not zero");
                    const ScalarType unit = dot1.unit();
                    squaredNormC += squaredNormB - abs(dot1);
                    dot2 -= dot3 * unit;
                    dot1 -= ScalarType(2) * squaredNormB * unit;
                    continue;
                }
            }
            /* Step 6 */ {
                const bool cond1 = abs(dot2) > squaredNormA;
                const bool cond2 = scalarNear(dot2, squaredNormA, precision) && (ScalarType(2) * dot1 < dot3);
                const bool cond3 = scalarNear(dot2, -squaredNormA, precision) && dot3.isNegative();
                if (cond1 || cond2 || cond3) {
                    assert(!dot2.isZero() && "[Error]: Condition 1~3 shall ensure dot2 is not zero");
                    const ScalarType unit = dot2.unit();
                    squaredNormC += squaredNormA - abs(dot2);
                    dot1 -= dot3 * unit;
                    dot2 -= ScalarType(2) * squaredNormA * unit;
                    continue;
                }
            }
            /* Step 7 */ {
                const bool cond1 = abs(dot3) > squaredNormA;
                const bool cond2 = scalarNear(dot3, squaredNormA, precision) && (ScalarType(2) * dot1 < dot2);
                const bool cond3 = scalarNear(dot3, -squaredNormA, precision) && dot2.isNegative();
                if (cond1 || cond2 || cond3) {
                    assert(!dot3.isZero() && "[Error]: Condition 1~3 shall ensure dot3 is not zero");
                    const ScalarType unit = dot3.unit();
                    squaredNormB += squaredNormA - abs(dot3);
                    dot1 -= dot2 * unit;
                    dot3 -= ScalarType(2) * squaredNormA * unit;
                    continue;
                }
            }
            /* Step 8 */ {
                const ScalarType temp = dot1 + dot2 + dot3 + squaredNormA + squaredNormB;
                const bool cond1 = temp.isNegative();
                const bool cond2 = scalarNear(temp, ScalarType(0), precision) && (ScalarType(2) * (squaredNormA + dot2) + dot3).isPositive();
                if (cond1 || cond2) {
                    squaredNormC += temp;
                    dot1 += ScalarType(2) * squaredNormB + dot3;
                    dot2 += ScalarType(2) * squaredNormA + dot3;
                    continue;
                }
            }
        } while(false);
        const ScalarType normA = sqrt(squaredNormA);
        const ScalarType normB = sqrt(squaredNormB);
        const ScalarType normC = sqrt(squaredNormC);
        const ScalarType alpha = arccos(dot1 / (ScalarType(2) * normB * normC));
        const ScalarType beta = arccos(dot2 / (ScalarType(2) * normA * normC));
        const ScalarType gamma = arccos(dot3 / (ScalarType(2) * normA * normB));
        return makeLattice(normA, normB, normC, alpha, beta, gamma);
    }
    /**
     * A directly simplified version of above function
     */
    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::LatticeMatrix PeriodicCell<ScalarType, Dim>::makeNiggliLattice2D(
            const LatticeMatrix& lattice, unsigned int maxIteration) {
        assert(maxIteration > 0 && "[Error]: Set maxIteration = 0 does nothing");
        ScalarType squaredNormA = lattice.row(0).squaredNorm();
        ScalarType squaredNormB = lattice.row(1).squaredNorm();
        ScalarType dot = ScalarType(2) * (lattice.row(0).asVector() * lattice.row(1).asVector());
        unsigned int iteration = 0;
        do {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Cannot reduce to niggli cell within required iterations");
            iteration += 1;

            if (squaredNormA > squaredNormB)
                squaredNormA.swap(squaredNormB);
            dot = -abs(dot);

            if (abs(dot) > squaredNormA) {
                assert(!dot.isZero() && "[Error]: Condition 1~3 shall ensure dot is not zero");
                const ScalarType unit = dot.unit();
                squaredNormB += squaredNormA - abs(dot);
                dot -= ScalarType(2) * squaredNormA * unit;
                continue;
            }
        } while(false);
        const ScalarType normA = sqrt(squaredNormA);
        const ScalarType normB = sqrt(squaredNormB);
        const ScalarType gamma = arccos(dot / (ScalarType(2) * normA * normB));
        return makeLattice2D(normA, normB, lattice(2, 2), gamma);
    }

    template<class ScalarType, unsigned int Dim>
    typename PeriodicCell<ScalarType, Dim>::LatticeMatrix PeriodicCell<ScalarType, Dim>::makeRepLattice(const LatticeMatrix& lattice) {
        LatticeMatrix result{};
        result.row(0) = lattice.row(1).crossProduct(lattice.row(2));
        result.row(1) = lattice.row(2).crossProduct(lattice.row(0));
        result.row(2) = lattice.row(0).crossProduct(lattice.row(1));
        const ScalarType factor = ScalarType(2 * M_PI) / (lattice.row(0) * result.row(0).asVector());
        result *= factor;
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    ScalarType PeriodicCell<ScalarType, Dim>::getVolume(const LatticeMatrix& lattice) {
        if constexpr (Dim == 1)
            return abs(lattice(0, 0));
        else if constexpr (Dim == 2)
            return (lattice.row(0).crossProduct(lattice.row(1))).compute().norm();
        else
            return abs(VectorType(lattice.row(0).crossProduct(lattice.row(1))) * lattice.row(2).asVector());
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
    PeriodicCell<ScalarType, Dim>::estimateRange(const LatticeMatrix& lattice, ValueType cutoff) {
        const auto repLatt = makeRepLattice(lattice);
        const ValueType factor = cutoff * ValueType(1 / (2 * M_PI));
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
    void PeriodicCell<ScalarType, Dim>::normalize_direct() {
        for (auto& elem : pos.asArray()) {
            elem -= floor(elem);
            assert(ScalarType(0) <= elem && elem <= ScalarType(1));
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
    template<ExtendCellOption Option>
    void PeriodicCell<ScalarType, Dim>::toSuperCellDirect(unsigned int x, unsigned int y, unsigned int z) {
        assert(type == Type::Direct);
        toSuperPosDirect<Option>(pos, x, y, z);
        auto rowX = lattice.row(0);
        rowX.asVector() *= ScalarType(x);
        auto rowY = lattice.row(1);
        rowY.asVector() *= ScalarType(y);
        auto rowZ = lattice.row(2);
        rowZ.asVector() *= ScalarType(z);
        if constexpr (ScalarType::isReverseDiff)
            lattice.makeContinuous();
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toUnitCellDirect(unsigned int x, unsigned int y, unsigned int z) {
        assert(type == Type::Direct);
        toUnitPosDirect(pos, x, y, z);
        const ScalarType inv_x = Core::reciprocal(ScalarType(x));
        const ScalarType inv_y = Core::reciprocal(ScalarType(y));
        const ScalarType inv_z = Core::reciprocal(ScalarType(z));

        auto rowX = lattice.row(0);
        rowX.asVector() *= inv_x;
        auto rowY = lattice.row(1);
        rowY.asVector() *= inv_y;
        auto rowZ = lattice.row(2);
        rowZ.asVector() *= inv_z;
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toDirect(PositionMatrix& target, const InvLatticeMatrix& invLattice) {
        target *= invLattice;
    }

    template<class ScalarType, unsigned int Dim>
    template<ExtendCellOption Option>
    void PeriodicCell<ScalarType, Dim>::toSuperPosDirect(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z) {
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
                            new_pos(index, 0) = (target(i, 0) + ScalarType(x_)) / ScalarType(x);
                            new_pos(index, 1) = (target(i, 1) + ScalarType(y_)) / ScalarType(y);
                            new_pos(index, 2) = (target(i, 2) + ScalarType(z_)) / ScalarType(z);
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
                            new_pos(index, 0) = (target(i, 0) + ScalarType(x_)) / ScalarType(x);
                            new_pos(index, 1) = (target(i, 1) + ScalarType(y_)) / ScalarType(y);
                            new_pos(index, 2) = (target(i, 2) + ScalarType(z_)) / ScalarType(z);
                            ++index;
                        }
                    }
                }
            }
        }
        else
            static_assert(Option == ExtendCellOption::AtomMajor, "[Error]: Invalid option");
        target.swap(new_pos);
        if constexpr (ScalarType::isReverseDiff)
            target.makeContinuous();
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::toUnitPosDirect(PositionMatrix& target, unsigned int x, unsigned int y, unsigned int z) {
        assert(Dim == 3);
        assert(x > 0 && y > 0 && z > 0);
        assert(target.getRow() % (x * y * z) == 0);

        auto colX = target.col(0);
        colX.asVector() *= ScalarType(x);
        auto colY = target.col(1);
        colY.asVector() *= ScalarType(y);
        auto colZ = target.col(2);
        colZ.asVector() *= ScalarType(z);

        const size_t numParticle = target.getRow();
        const size_t newNumParticle = numParticle / (x * y * z);
        PositionMatrix new_pos(newNumParticle, Dim);
        size_t toFill = 0;
        size_t toCheck = 0;
        const ScalarType one = ScalarType(1);
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
}
