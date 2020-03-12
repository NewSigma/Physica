#ifndef PHYSICA_EXPRREADER_H
#define PHYSICA_EXPRREADER_H

#include <list>
#include "Numerical.h"

class ExprReader {
private:
    std::list<std::wstring> anti_poland;
public:
    ExprReader(const std::wstring& str);
    Numerical* calc();
private:
    static bool isSign(wchar_t c);
};

#endif
