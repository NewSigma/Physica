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

#include <iostream>

namespace Physica::Gui {
    template<class MatrixType>
    bool ContourSeries<MatrixType>::Quad::operator==(Quad quad) const noexcept {
        return row == quad.row && column == quad.column;
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad::Vertex ContourSeries<MatrixType>::Quad::getTopLeftVertex() const {
        return {row, column};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad::Vertex ContourSeries<MatrixType>::Quad::getTopRightVertex() const {
        return {row, column + 1};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad::Vertex ContourSeries<MatrixType>::Quad::getBottomLeftVertex() const {
        return {row + 1, column};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad::Vertex ContourSeries<MatrixType>::Quad::getBottomRightVertex() const {
        return {row + 1, column + 1};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad ContourSeries<MatrixType>::Quad::topNeigh() const {
        assert(row > 0);
        return {row - 1, column};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad ContourSeries<MatrixType>::Quad::bottomNeigh() const {
        return {row + 1, column};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad ContourSeries<MatrixType>::Quad::leftNeigh() const {
        assert(column > 0);
        return {row, column - 1};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Quad ContourSeries<MatrixType>::Quad::rightNeigh() const {
        return {row, column + 1};
    }

    template<class MatrixType>
    void ContourSeries<MatrixType>::Edge::moveToNextEdge() {
        switch (dir) {
            case Edge::Top:
                quad = quad.topNeigh();
                dir = Edge::Bottom;
                break;
            case Edge::Bottom:
                quad = quad.bottomNeigh();
                dir = Edge::Top;
                break;
            case Edge::Left:
                quad = quad.leftNeigh();
                dir = Edge::Right;
                break;
            case Edge::Right:
                quad = quad.rightNeigh();
                dir = Edge::Left;
        }
    }
    /**
     * Label 1 and 2 is by clockwise
     */
    template<class MatrixType>
    typename ContourSeries<MatrixType>::Edge::Vertex ContourSeries<MatrixType>::Edge::getVertex1() const {
        switch (dir) {
            case Edge::Top:
                return quad.getTopLeftVertex();
            case Edge::Bottom:
                return quad.getBottomRightVertex();
            case Edge::Left:
                return quad.getBottomLeftVertex();
            case Edge::Right:
                return quad.getTopRightVertex();
        }
        assert(false);
        return {0, 0};
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::Edge::Vertex ContourSeries<MatrixType>::Edge::getVertex2() const {
        switch (dir) {
            case Edge::Top:
                return quad.getTopRightVertex();
            case Edge::Bottom:
                return quad.getBottomLeftVertex();
            case Edge::Left:
                return quad.getTopLeftVertex();
            case Edge::Right:
                return quad.getBottomRightVertex();
        }
        assert(false);
        return {0, 0};
    }

    template<class MatrixType>
    bool ContourSeries<MatrixType>::Grid::canInterpolate(Edge edge, double level) const noexcept {
        size_t vertex1_row = 0, vertex1_column = 0, vertex2_row = 0, vertex2_column = 0;
        const size_t quad_row = edge.quad.row;
        const size_t quad_column = edge.quad.column;
        switch (edge.dir) {
            case Edge::Top:
                vertex1_row = vertex2_row = quad_row;
                vertex1_column = quad_column;
                vertex2_column = quad_column + 1;
                break;
            case Edge::Bottom:
                vertex1_row = vertex2_row = quad_row + 1;
                vertex1_column = quad_column;
                vertex2_column = quad_column + 1;
                break;
            case Edge::Left:
                vertex1_column = vertex2_column = quad_column;
                vertex1_row = quad_row;
                vertex2_row = quad_row + 1;
                break;
            case Edge::Right:
                vertex1_column = vertex2_column = quad_column + 1;
                vertex1_row = quad_row;
                vertex2_row = quad_row + 1;
        }
        const auto value1 = double(z(vertex1_row, vertex1_column));
        const auto value2 = double(z(vertex2_row, vertex2_column));
        return (value1 >= level && level >= value2) || (value1 <= level && level <= value2);
    }

    template<class MatrixType>
    ContourSeries<MatrixType>::Grid::Grid(const Core::LValueMatrix<MatrixType>& x_,
                                          const Core::LValueMatrix<MatrixType>& y_,
                                          const Core::LValueMatrix<MatrixType>& z_) : x(x_), y(y_), z(z_), flags() {
        flags.resize(x.getRow() - 1);
        for (size_t i = 0; i < x.getRow() - 1; ++i)
            flags[i].resize(x.getColumn() - 1, false);
    }

    template<class MatrixType>
    typename ContourSeries<MatrixType>::ContourLine ContourSeries<MatrixType>::Grid::interpolateFromEdge(Edge edge, double level) {
        assert(!haveVisited(edge.quad));
        ContourLine line{};
        /* Initialize */ {
            interpolateEdge(edge, level, line);
        }

        const Quad initialQuad = edge.quad;
        while (true) {
            const Quad quad = edge.quad;
            if (haveVisited(quad))
                break;
            setVisited(quad);

            if (edge.dir != Edge::Top && canInterpolate(Edge(quad, Edge::Top), level)) {
                interpolateEdge(edge, level, line);
                edge.dir = Edge::Top;
            }
            else if (edge.dir != Edge::Right && canInterpolate(Edge(quad, Edge::Right), level)) {
                interpolateEdge(edge, level, line);
                edge.dir = Edge::Right;
            }
            else if (edge.dir != Edge::Bottom && canInterpolate(Edge(quad, Edge::Bottom), level)) {
                interpolateEdge(edge, level, line);
                edge.dir = Edge::Bottom;
            }
            else if (edge.dir != Edge::Left && canInterpolate(Edge(quad, Edge::Left), level)) {
                interpolateEdge(edge, level, line);
                edge.dir = Edge::Left;
            }
            else
                break;

            if (isBoundary(edge))
                break;
            edge.moveToNextEdge();
        }
        
        const bool backToStartPoint = initialQuad == edge.quad;
        if (backToStartPoint) {
            const double factor = 1 + std::numeric_limits<double>::epsilon(); //Multiply with a factor because start point and end point cannot be same in Qt
            line.first.append(line.first[0] * factor);
            line.second.append(line.second[0] * factor);
        }
        return ContourLine(std::move(line));
    }

    template<class MatrixType>
    bool ContourSeries<MatrixType>::Grid::isBoundary(Edge edge) const noexcept {
        return (edge.quad.row == 0 && edge.dir == Edge::Top)
            || (edge.quad.column == 0 && edge.dir == Edge::Left)
            || (edge.quad.row == (x.getRow() - 2) && edge.dir == Edge::Bottom)
            || (edge.quad.column == (x.getColumn() - 2) && edge.dir == Edge::Right);
    }

    template<class MatrixType>
    void ContourSeries<MatrixType>::Grid::interpolateEdge(Edge edge, double level, ContourLine& line) const {
        auto vertex1 = edge.getVertex1();
        const double z1 = z(vertex1.first, vertex1.second);

        auto vertex2 = edge.getVertex2();
        const double z2 = z(vertex2.first, vertex2.second);

        const double factor = (level - z1) / (z2 - z1);
        const double x1 = x(vertex1.first, vertex1.second);
        const double x2 = x(vertex2.first, vertex2.second);
        const double point_x = x1 + (x2 - x1) * factor;
        line.first.append(point_x);

        const double y1 = y(vertex1.first, vertex1.second);
        const double y2 = y(vertex2.first, vertex2.second);
        const double point_y = y1 + (y2 - y1) * factor;
        line.second.append(point_y);
    }

    template<class MatrixType>
    ContourSeries<MatrixType>::ContourSeries(const Core::LValueMatrix<MatrixType>& x,
                                             const Core::LValueMatrix<MatrixType>& y,
                                             const Core::LValueMatrix<MatrixType>& z,
                                             Utils::Array<double> levels,
                                             QObject* parent) : QObject(parent) {
        for (double level : levels) {
            Grid grid = Grid(x, y, z);
            initContourLine(grid, level);
        }
    }

    template<class MatrixType>
    void ContourSeries<MatrixType>::attachTo(QChart& chart) {
        QObject::setParent(&chart);
        for (const auto& line : contourLines) {
            const auto& x = line.first;
            const auto& y = line.second;
            QLineSeries* series = new QLineSeries();
            for (size_t i = 0; i < x.getLength(); ++i)
                *series << QPointF(x[i], y[i]);
            splines.append(series);
            chart.addSeries(series);
        }
    }

    template<class MatrixType>
    void ContourSeries<MatrixType>::initContourLine(Grid& grid, double level) {
        for (size_t i = 0; i < grid.getRow(); ++i)
            tryInterpolate(grid, level, Edge(Quad(i, 0), Edge::Left));

        for (size_t i = 0; i < grid.getRow(); ++i)
            tryInterpolate(grid, level, Edge(Quad(i, grid.getColumn() - 1), Edge::Right));

        for (size_t i = 0; i < grid.getColumn(); ++i)
            tryInterpolate(grid, level, Edge(Quad(0, i), Edge::Top));

        for (size_t i = 0; i < grid.getColumn(); ++i)
            tryInterpolate(grid, level, Edge(Quad(grid.getRow() - 1, i), Edge::Bottom));
        
        for (size_t i = 1; i < grid.getColumn() - 1; ++i) {
            for (size_t j = 1; j < grid.getRow() - 1; ++j) {
                const Quad quad = Quad(j, i);
                tryInterpolate(grid, level, Edge(quad, Edge::Top));
                tryInterpolate(grid, level, Edge(quad, Edge::Bottom));
                tryInterpolate(grid, level, Edge(quad, Edge::Left));
                tryInterpolate(grid, level, Edge(quad, Edge::Right));
            }
        }
    }

    template<class MatrixType>
    void ContourSeries<MatrixType>::tryInterpolate(Grid& grid, double level, Edge edge) {
        if (!grid.haveVisited(edge.quad) && grid.canInterpolate(edge, level)) {
            auto line = grid.interpolateFromEdge(edge, level);
            if (line.first.getLength() > 1)
                contourLines.append(std::move(line));
        }
    }
}
