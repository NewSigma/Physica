#ifndef PHYSICA_ARRAYWRITER_H
#define PHYSICA_ARRAYWRITER_H

#include <GL/gl.h>

class ArrayWriter {
public:
    ArrayWriter(GLfloat* array, int length);
    void writeNext(float a, float b, float c);
private:
    float* array;
    int length;
    int index;
};

#endif
