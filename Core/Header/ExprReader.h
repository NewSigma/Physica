#ifndef PHYSICA_EXPRREADER_H
#define PHYSICA_EXPRREADER_H

#include <list>
#include "RealNumber.h"

class ExprReader {
private:
    std::list<std::string> anti_poland;
public:
    ExprReader(const std::string& s);

    RealNumber* calc();
private:
    void clearString(std::string& s);
};

#endif
