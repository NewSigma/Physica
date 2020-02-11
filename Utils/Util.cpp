#include <cstring>
#include <iostream>
#include "Header/Util.h"

void error(const char* file, const unsigned long line, const char* msg) {
    std::cout << '[';
    unsigned long begin, end;
    for(end = strlen(file) - 1; end  >= 0; --end)
        if(file[end] == '.')
            break;

    for(begin = end - 1; begin >= 0; --begin)
        if(file[begin] == '/')
            break;

    for(++begin; begin < end; ++begin)
        std::cout << file[begin];
    std::cout << ": " << line << "] " << msg << '\n';
}