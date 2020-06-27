#ifndef PHYSICA_EXPRREADER_H
#define PHYSICA_EXPRREADER_H

#include <list>
#include <iosfwd>
#include <Physica/Core/MultiPrecition/Scalar.h>

using Physica::Core::MultiScalar;

namespace Physica::Interpreter {
    class ExprReader {
    private:
        std::list<std::wstring> anti_poland;
    public:
        explicit ExprReader(const std::wstring& str);
        MultiScalar calc();
    private:
        static bool isSign(wchar_t c);
    };
}

#endif
