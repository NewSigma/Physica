#ifndef PHYSICA_EXPRREADER_H
#define PHYSICA_EXPRREADER_H

#include <list>
#include "RealNumber.h"

class ExprReader {
private:
    std::list<std::wstring> anti_poland;
public:
    ExprReader(const std::wstring& s);

    RealNumber* calc();
private:
    void clearString(std::wstring& s);
};

#endif
