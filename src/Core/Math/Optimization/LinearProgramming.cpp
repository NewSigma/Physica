/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#include "Physica/Core/Math/Optimization/LinearProgramming.h"
#include <iostream>
namespace Physica::Core {
    LinearProgramming::LinearProgramming() : order(nullptr), state(Normal) {}

    LinearProgramming::LinearProgramming(LinearProgramming&& lp) noexcept
            : target(std::move(lp.target)), data(std::move(lp.data)), order(lp.order), state(lp.state) {
        lp.order = nullptr;
    }

    LinearProgramming& LinearProgramming::operator=(const LinearProgramming& lp) {
        if(this != &lp) {
            this->~LinearProgramming();
            target = lp.target;
            data = lp.data;
            const auto length = data.getRow() + data.getColumn() + 1;
            order = new size_t[length];
            for(size_t i = 0; i < length; ++i)
                order[i] = lp.order[i];
            state = lp.state;
        }
        return *this;
    }

    LinearProgramming& LinearProgramming::operator=(LinearProgramming&& lp) noexcept {
        this->~LinearProgramming();
        target = std::move(lp.target);
        data = std::move(lp.data);
        order = lp.order;
        lp.order = nullptr;
        state = lp.state;
        return *this;
    }

    LinearProgramming::~LinearProgramming() {
        delete[] order;
    }
    /*!
     * Add a restrain to the LinearProgramming, return true if successful, false if failed.
     * Fails when sizes are not compatible.
     */
    bool LinearProgramming::addRestrain(Vector<> v, RestrainType type) {
        const auto length = target.getLength();
        if(length == 0 || length != v.getLength())
            return false;
        switch(type) {
            case Equal:
                data.appendRow((-v).calc());
                data.appendRow(std::move(v));
                break;
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wimplicit-fallthrough="
            case LessEqual:
                v.toOpposite();
            case GreaterEqual:
                data.appendRow(std::move(v));
    #pragma GCC diagnostic push
        }
        return true;
    }
    /*!
     * Variable are non-negative by default, if a variable is able to be negative
     * , call this function.
     */
    void LinearProgramming::forceNonNegative(size_t index) {
        if(index == 0 || index >= target.getLength())
            return;
        MultiScalar copy(target[index]);
        target << std::move(copy);
        const auto l = data.getRow();
        for(size_t i = 0; i < l; ++i) {
            auto& temp = data[i];
            copy = temp[index];
            temp << std::move(copy);
        }
    }

    void LinearProgramming::solve() {
        initialize();
        if(state == Normal)
            solveImpl();
    }
    /*!
     * This function successes when the target is null or the new target has the same amount
     * of elements.
     * Return true if successful, false if failed.
     */
    bool LinearProgramming::setTarget(const Vector<>& v) {
        const auto length = target.getLength();
        if(length == 0 || length == v.getLength()) {
            target = v;
            return true;
        }
        return false;
    }
    /*!
     * This function successes when the target is null or the new target has the same amount
     * of elements.
     * Return true if successful, false if failed.
     */
    bool LinearProgramming::setTarget(Vector<>&& v) {
        const auto length = target.getLength();
        if(length == 0 || length == v.getLength()) {
            target = std::move(v);
            return true;
        }
        return false;
    }
    /*!
     * initialize() generates the variable order and forces the initial solution valid.
     */
    void LinearProgramming::initialize() {
        const auto row = data.getRow();
        const auto column = data.getColumn();
        /*
         * Allocated 1 more element, which comes from constant term.
         * This extra space will be used for a temp variable.
         */
        const auto orderLength = row + column;

        //Return directly if the initial basic solution is feasiable.
        size_t minimumIndex = findMinConst();
        if(!data(minimumIndex, 0).isNegative()) {
            order = new size_t[orderLength];
            for(size_t i = 0; i < orderLength; ++i)
                order[i] = i + 1;
            return;
        }

        /* Generate Order */ {
            order = new size_t[orderLength];
            for(size_t i = 0; i < column - 1; ++i)
                order[i] = i + 1;
            order[column - 1] = 0;
            for(size_t i = column; i < orderLength; ++i)
                order[i] = i;
        }
        //Save the old target.
        Vector<> copyTarget(std::move(target));
        //Construct the auxiliary problem.
        target = Vector<>::zeroVector(column + 1);
        target[column] = MultiScalar((SignedMPUnit)-1);
        for(size_t i = 0; i < row; ++i)
            data[i] << MultiScalar::getOne();
        pivot(minimumIndex, column - 1);
        //Solve the auxiliary problem.
        solveImpl();
        if(target[0].getPower() >= 0) { //Scalar whose power < 0 is a very small value, approximately equal to zero.
            state = Infeasiable;
            return;
        }
        //Find the position of temp term.
        size_t tempTermOrder;
        for(tempTermOrder = 0; tempTermOrder < orderLength; ++tempTermOrder)
            if(order[tempTermOrder] == 0)
                break;
        Q_ASSERT(tempTermOrder != orderLength);
        //If temp term is a basic variable, exchange it with a non-basic variable.
        if(tempTermOrder >= column) {
            size_t i;
            for(i = 1; i <= column; ++i)
                if(!data(tempTermOrder - column, i).isZero())
                    break;
            --i;
            pivot(tempTermOrder - column, i);
            tempTermOrder = i;
        }
        data.removeColumnAt(tempTermOrder + 1);
        //Remove 0 from order
        memmove(order + tempTermOrder, order + tempTermOrder + 1
                , (orderLength - tempTermOrder - 1) * sizeof(size_t));
        //Restore the original target.
        target = Vector<>::zeroVector(column);
        target[0] = copyTarget[0];
        for(size_t i = 1; i < column; ++i) {
            //Is variable i basic?
            bool isNonBasic = false;
            for(size_t j = 0; j < column - 1; ++j) {
                if(i == order[j]) {
                    target[j + 1] += copyTarget[i];
                    isNonBasic = true;
                    break;
                }
            }

            if(!isNonBasic) {
                for(size_t j = column - 1; j < orderLength - 1; ++j) {
                    if(i == order[j]) {
                        isNonBasic = true; //Set to true for the convenience of debug.
                        target += data[j - column + 1] * copyTarget[i];
                        break;
                    }
                }
            }
            Q_ASSERT(isNonBasic);
        }
    }
    /*!
     * Exchange the basic variable at index \from and non-basic variable at index \to.
     *
     * Optimize: Unnecessary multiplies on column \to.
     */
    void LinearProgramming::pivot(size_t basic, size_t nonBasic) {
        const auto dataSize = data.getRow();
        Q_ASSERT(basic < dataSize && nonBasic < target.getLength());
        /* Handle order */
        std::swap(order[data.getColumn() - 1 + basic], order[nonBasic]);
        /* Handle row basic */
        //nonBasic starts from 0, while column 0 is consisted by constants.
        ++nonBasic;
        auto& fromRow = data[basic];
        auto temp = reciprocal(data(basic, nonBasic));
        fromRow *= -temp;
        auto& ele = data(basic, nonBasic);
        ele = std::move(temp);

        MultiScalar copy;
        for(size_t i = 0; i < dataSize; ++i) {
            if(i == basic)
                continue;
            copy = data(i, nonBasic);
            data[i] += fromRow * data(i, nonBasic);
            data(i, nonBasic) = copy * ele;
        }
        copy = target[nonBasic];
        target += fromRow * target[nonBasic];
        target[nonBasic] = copy * ele;
    }

    void LinearProgramming::solveImpl() {
        const auto row = data.getRow();
        const auto length = target.getLength();
        size_t candidate = findPositiveVariable();
        while(candidate < length) {
            size_t exchange = 0;
            //Select a acceptable start.
            while(exchange < row && !data(exchange, candidate).isNegative())
                ++exchange;
            //Select the most strict restrain.
            for(size_t i = exchange + 1; i < row; ++i) {
                auto& temp = data(i, candidate);
                if(temp.isNegative()) {
                    bool b = (data(exchange, 0) / data(exchange, candidate)) < (data(i, 0) / temp);
                    exchange = b ? i : exchange;
                }
            }
            //Return if no restrain is selected.
            if(exchange == row || !data(exchange, candidate).isNegative()) {
                state = Infinity;
                return;
            }
            pivot(exchange, candidate - 1);
            //Update condition.
            candidate = findPositiveVariable();
        }
    }

    size_t LinearProgramming::findMinConst() const {
        const auto row = data.getRow();
        Q_ASSERT(row > 0);
        size_t minimumIndex = 0; //Assume element at 0 is minimum
        for(size_t i = 1; i < row; ++i)
            minimumIndex = data(minimumIndex, 0) < data(i, 0) ? minimumIndex : i;
        return minimumIndex;
    }

    size_t LinearProgramming::findPositiveVariable() const {
        const auto length = target.getLength();
        for(size_t i = 1; i < length; ++i)
            if(target[i].isPositive())
                return i;
        return length;
    }
}