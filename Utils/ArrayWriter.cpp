/*
 * Write several elements into a array, which makes it easy to construct vertex attribute of OpenGL.
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <GL/gl.h>
#include <iostream>
#include "Header/ArrayWriter.h"

ArrayWriter::ArrayWriter(GLfloat* arr, const int len) {
    array = arr;
    length = len;
    index = 0;
}

void ArrayWriter::writeNext(float a, float b, float c) {
    if(length - 1 >= index) {
        std::cout << "[ArrayWriter]: Out of index!" << std::endl;
        return;
    }
    array[index] = a;
    array[index + 1] = b;
    array[index + 2] = c;
    index += 3;
}