/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PHYSICATESTS_H
#define PHYSICA_PHYSICATESTS_H

class Numerical;

void checkGPU();
void const_test();
void elementary_function_test();
void numerical_test();
void printElements(const Numerical& n);
void simpleNet();

#endif