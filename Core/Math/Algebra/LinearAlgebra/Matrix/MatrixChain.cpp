/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/MatrixChain.h"
#include "Core/Header/Matrix.h"

namespace Physica::Core {
    /*!
     * \class MatrixChain provides a effective method to calculate the production of a chain of matrices.
     *
     * Note: matrices should be initialized and pass into the MatrixChain before calling solve(),
     * matrices must not be deleted before the calculation finished.
     *
     * Warning:
     * 1.If the chain is too long, size_t or stack may overflow and lead to wrong results.
     * 2.length should be at least 2 or solve() will simply return the matrix passed in
     * , which may result in double free.
     *
     * Reference: 算法导论 第三版 Page: 210-215
     */
    MatrixChain::MatrixChain(size_t length)
            : chain(new Matrix*[length]), length(length), price(new size_t*[length]), point(new size_t*[length - 1]) {
        const auto length_1 = length - 1;
        for(size_t i = 0; i < length_1; ++i) {
            price[i] = new size_t[length];
            point[i] = new size_t[length_1];
        }
        //price is longer than point by 1.
        price[length_1] = new size_t[length];
    }

    MatrixChain::~MatrixChain() {
        delete[] chain;
        const auto length_1 = length - 1;
        for(size_t i = 0; i < length_1; ++i) {
            delete[] price[i];
            delete[] point[i];
        }
        //price is longer than point by 1.
        delete[] price[length_1];
        delete[] price;
        delete[] point;
    }
    //!Optimize: Only half of the space of price and point is used. Maybe change them into a 1D array.
    std::unique_ptr<Matrix> MatrixChain::solve() {
        for(size_t i = 0; i < length; ++i)
            price[i][i] = 0;

        for(size_t l = 1; l < length; ++l) { //(l + 1) is the length of sub-chain.
            const auto m_max = length - l;
            for(size_t m = 0; m < m_max; ++m) { //m is the start index of each sub-chain.
                const auto m_end = m + l; //End index of each sub-chain.
                auto n = m; //Cut (m ... ) into (m ... n) and (n + 1 ... m_end).
                auto chain_n = chain[n];
                /* Handle n = m */ {
                    price[m][m_end] = price[n + 1][m_end] + chain_n->row() * chain_n->column() * chain[n + 1]->row();
                    point[m][m_end] = n;
                    ++n;
                }
                for(; n < m_end; ++n) {
                    chain_n = chain[n];
                    size_t temp = price[m][n] + price[n + 1][m_end]
                                                + chain_n->row() * chain_n->column() * chain[n + 1]->row();
                    price[m][m_end] = temp < price[m][m_end] ? temp : price[m][m_end];
                    point[m][m_end] = temp < price[m][m_end] ? n : point[m][m_end];
                }
                /* Handle n = m_end */ {
                    chain_n = chain[n];
                    size_t temp = price[m][n] + chain_n->row() * chain_n->column() * chain[n + 1]->row();
                    price[m][m_end] = temp < price[m][m_end] ? temp : price[m][m_end];
                    point[m][m_end] = temp < price[m][m_end] ? n : point[m][m_end];
                }
            }
        }
        return multiply(0, length - 1);
    }
    //!Both \from and \to are included.
    std::unique_ptr<Matrix> MatrixChain::multiply(size_t from, size_t to) {
        if(from == to)
            return std::unique_ptr<Matrix>(chain[from]);
        const auto cut_at = point[from][to];
        auto first = multiply(from, cut_at);
        auto second = multiply(cut_at + 1, to);
        auto first_p = from == cut_at ? first.release() : first.get();
        auto second_p = to == cut_at ? second.release() : second.get();
        return (*first_p * *second_p);
    }
}