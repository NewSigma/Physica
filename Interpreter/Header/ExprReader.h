#ifndef PHYSICA_EXPRREADER_H
#define PHYSICA_EXPRREADER_H

#include <list>
#include <iosfwd>

namespace Physica::Core {
    class Scalar;
}

using Physica::Core::Scalar;

namespace Physica::Interpreter {
    class ExprReader {
    private:
        std::list<std::wstring> anti_poland;
    public:
        ExprReader(const std::wstring& str);
        Scalar calc();
    private:
        static bool isSign(wchar_t c);
    };
}

#endif
