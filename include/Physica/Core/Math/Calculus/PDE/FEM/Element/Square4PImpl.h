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

namespace Physica::Core {
    template<class Scalar>
    Scalar Square4P<Scalar>::shapePartialS1(ShapeIndex shapeIndex, const Poing<2>& p) {
        Scalar temp;
        switch (shapeIndex) {
            case ShapeIndex::Shape1:
                temp = -(1 - p.y());
            case ShapeIndex::Shape2:
                temp = -(1 + p.y());
            case ShapeIndex::Shape3:
                temp = 1 + p.y();
            case ShapeIndex::Shape4:
                temp = 1 - p.y(); 
        }
        return temp / 4;
    }

    template<class Scalar>
    Scalar Square4P<Scalar>::shapePartialS2(ShapeIndex shapeIndex, const Poing<2>& p) {
        Scalar temp;
        switch (shapeIndex) {
            case ShapeIndex::Shape1:
                temp = -(1 - p.x());
            case ShapeIndex::Shape2:
                temp = 1 - p.x();
            case ShapeIndex::Shape3:
                temp = 1 + p.x();
            case ShapeIndex::Shape4:
                temp = -(1 + p.x()); 
        }
        return temp / 4;
    }

    template<class Scalar>
    Matrix Square4P<Scalar>::jacobi(const Poing<2>& p) {
        return Matrix({{shapePartialS1(Shape1, p) * (*nodes[0]).x()
                        + shapePartialS1(Shape2, p) * (*nodes[1]).x()
                        + shapePartialS1(Shape3, p) * (*nodes[2]).x()
                        + shapePartialS1(Shape4, p) * (*nodes[3]).x()
                        , shapePartialS1(Shape1, p) * (*nodes[0]).y()
                        + shapePartialS1(Shape2, p) * (*nodes[1]).y()
                        + shapePartialS1(Shape3, p) * (*nodes[2]).y()
                        + shapePartialS1(Shape4, p) * (*nodes[3]).y()}
                            , {shapePartialS2(Shape1, p) * (*nodes[0]).x()
                            + shapePartialS2(Shape2, p) * (*nodes[1]).x()
                            + shapePartialS2(Shape3, p) * (*nodes[2]).x()
                            + shapePartialS2(Shape4, p) * (*nodes[3]).x()
                            , shapePartialS2(Shape1, p) * (*nodes[0]).y()
                            + shapePartialS2(Shape2, p) * (*nodes[1]).y()
                            + shapePartialS2(Shape3, p) * (*nodes[2]).y()
                            + shapePartialS2(Shape4, p) * (*nodes[3]).y()}});
    }
}