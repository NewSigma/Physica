/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "CalcBasic.h"
#include "ArraySupport.h"
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    //Reference: GMP Doc BaseCase Multiplication
    Numerical square(const Numerical& n) {
        if(n == BasicConst::getInstance().get_1())
            return Numerical(n);
        else {
            auto n_size = n.getSize();
            //Estimate the ed of result first. we will calculate it accurately later.
            auto length = 2 * n_size;
            auto byte = reinterpret_cast<NumericalUnit*>(calloc(length, sizeof(NumericalUnit)));

            for(int i = 0; i < n_size - 1; ++i)
                byte[i + n_size] = mulAddArrByWord(byte + i + i + 1, n.byte + i + 1, n_size - i - 1, n.byte[i]);
            //Shift count is known, possible to optimize the performance.
            byteLeftShift(byte, length, 1);

            NumericalUnit high, low;
            for(int i = 0; i < n_size; ++i) {
                mulWordByWord(high, low, n.byte[i], n.byte[i]);
                byte[2 * i] += low;
                byte[2 * i + 1] += high;
            }
            auto power = n.power * 2 + 1;
            if(high == 0) {
                --power;
                byte = reinterpret_cast<NumericalUnit*>(realloc(byte, (length - 1) * sizeof(NumericalUnit)));
            }
            return Numerical(byte, length, power);
        }
    }
}