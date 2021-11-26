/*
 * Copyright 2021 WeiBo He.
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

#include <QtCharts/QLineSeries>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

namespace Physica::Gui {
    template<class MatrixType>
    class ContourSeries : public QObject {
        using ScalarType = typename MatrixType::ScalarType;
        using ContourLine = std::pair<Utils::Array<double>, Utils::Array<double>>;

        struct Quad {
            using Vertex = std::pair<size_t, size_t>;

            size_t row;
            size_t column;

            Quad(size_t row_, size_t column_) : row(row_), column(column_) {}
            /* Operators */
            [[nodiscard]] bool operator==(Quad quad) const noexcept;
            /* Getters */
            [[nodiscard]] Vertex getTopLeftVertex() const;
            [[nodiscard]] Vertex getTopRightVertex() const;
            [[nodiscard]] Vertex getBottomLeftVertex() const;
            [[nodiscard]] Vertex getBottomRightVertex() const;
            [[nodiscard]] Quad topNeigh() const;
            [[nodiscard]] Quad bottomNeigh() const;
            [[nodiscard]] Quad leftNeigh() const;
            [[nodiscard]] Quad rightNeigh() const;
        };

        struct Edge {
            using Vertex = typename Quad::Vertex;

            enum Direction {
                Top,
                Bottom,
                Left,
                Right
            };

            Quad quad;
            Direction dir;

            Edge(Quad quad_, Direction dir_) : quad(quad_), dir(dir_) {}
            /* Operations */
            void moveToNextEdge();
            /* Getters */
            [[nodiscard]] Vertex getVertex1() const;
            [[nodiscard]] Vertex getVertex2() const;
        };

        class Grid {
            using FlagMatrix = Utils::Array<Utils::Array<bool>>;

            const Core::LValueMatrix<MatrixType>& x;
            const Core::LValueMatrix<MatrixType>& y;
            const Core::LValueMatrix<MatrixType>& z;
            FlagMatrix flags;
        public:            
            Grid(const Core::LValueMatrix<MatrixType>& x_,
                 const Core::LValueMatrix<MatrixType>& y_,
                 const Core::LValueMatrix<MatrixType>& z_);
            Grid(const Grid&) = delete;
            Grid(Grid&&) noexcept = delete;
            ~Grid() = default;
            /* Operations */
            Grid& operator=(const Grid&) = delete;
            Grid& operator=(Grid&&) noexcept = delete;
            ContourLine interpolateFromEdge(Edge edge, double level);
            /* Getters */
            [[nodiscard]] size_t getRow() const noexcept { assert(x.getRow() > 0); return x.getRow() - 1; }
            [[nodiscard]] size_t getColumn() const noexcept { assert(x.getColumn() > 0); return x.getColumn() - 1; }
            [[nodiscard]] bool canInterpolate(Edge edge, double level) const noexcept;
            [[nodiscard]] bool isBoundary(Edge edge) const noexcept;
            [[nodiscard]] bool haveVisited(Quad quad) const { return flags[quad.row][quad.column]; }
            /* Setters */
            void setVisited(Quad quad) { flags[quad.row][quad.column] = true; }
        private:
            void interpolateEdge(Edge edge, double level, ContourLine& line) const;
        };
        
        Utils::Array<ContourLine> contourLines;
        Utils::Array<QLineSeries*> splines;
    public:
        ContourSeries(const Core::LValueMatrix<MatrixType>& x,
                      const Core::LValueMatrix<MatrixType>& y,
                      const Core::LValueMatrix<MatrixType>& z,
                      Utils::Array<double> level,
                      QObject* parent = nullptr);
        ContourSeries(const ContourSeries&) = delete;
        ContourSeries(ContourSeries&&) noexcept = delete;
        ~ContourSeries() = default;
        /* Operators */
        ContourSeries& operator=(const ContourSeries&) = delete;
        ContourSeries& operator=(ContourSeries&&) noexcept = delete;
        /* Operations */
        void attachTo(QChart& chart);
    private:
        void initContourLine(Grid& grid, double level);
        void tryInterpolate(Grid& grid, double level, Edge edge);
    };
}

#include "ContourSeriesImpl.h"
